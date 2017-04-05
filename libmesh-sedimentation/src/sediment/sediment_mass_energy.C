/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/gmv_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/vtk_io.h"
#include "libmesh/getpot.h"
#include "libmesh/perf_log.h"

// This example will solve a linear transient system,
// so we need to include the \p TransientLinearImplicitSystem definition.
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/vector_value.h"


#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"


#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/analytic_function.h"
#include "libmesh/zero_function.h"

#include "stab_helper.h"

#include "define.h"

// The definition of a geometric element
#include "libmesh/elem.h"


#include "sedimentation_transport.h"

//#define MAX(a,b) a>b?a:b
//#define MIN(a,b) a<b?a:b

// Bring in everything from the libMesh namespace
using namespace libMesh;
using namespace std;

void SedimentationTransport::PrintMass(const char *fname) {

    double local[2] = {0.0, 0.0};
    double global[2] = {0.0, 0.0};

    const MeshBase & mesh = es.get_mesh();

    const unsigned int dim = mesh.mesh_dimension();
    TransientLinearImplicitSystem & system =
            es.get_system<TransientLinearImplicitSystem> ("transport");

    const unsigned int s_var = system.variable_number("s");

    const DofMap & dof_map = system.get_dof_map();

    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_face;

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = system.variable_type(s_var);

    // Build a Finite Element object of the specified type.  Since the
    // \p FEBase::build() member dynamically creates memory we will
    // store the object as an \p UniquePtr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
    UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_type));

    // A Gauss quadrature rule for numerical integration.
    // Let the \p FEType object decide what order rule is appropriate.
    QGauss qrule(dim, fe_type.default_quadrature_order());
    QGauss qface(dim - 1, fe_type.default_quadrature_order());

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);
    fe_face->attach_quadrature_rule(&qface);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.  We will start
    // with the element Jacobian * quadrature weight at each integration point.
    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<Real>& JxW_face = fe_face->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const std::vector<std::vector<Real> >& phi_face = fe_face->get_phi();

    const Real Us = es.parameters.get<Real>("Us");
    const Real dt = es.parameters.get<Real>("dt");


    MeshBase::const_element_iterator it = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_local_elements_end();

    for (; it != end; ++it) {
        
        const Elem *elem = *it;

        dof_map.dof_indices(elem, dof_indices, s_var);

        fe->reinit(elem);


        for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
            // Values to hold the old solution & its gradient.
            Number s = 0.0;

            // Compute the old solution & its gradient.
            for (unsigned int l = 0; l < phi.size(); l++) {
                s += phi[l][qp] * system.current_solution(dof_indices[l]);
                //s +=  system.current_solution(dof_indices[l]);
            }
            
            local[0] += JxW[qp] * s;
            

        }

        {


            for (unsigned int s = 0; s < elem->n_sides(); s++)
                if (elem->neighbor(s) == NULL) {

                    fe_face->reinit(elem, s);

                    // Applying flux advective boundary condition
                    if (mesh.boundary_info->boundary_id(elem, s) == this->deposition_id) {


                        for (unsigned int qp = 0; qp < qface.n_points(); qp++) {
                            Number s = 0.0;
                            for (unsigned int i = 0; i < phi_face.size(); i++) {
                                s += phi_face[i][qp] * system.old_solution(dof_indices[i]);
                            }
                           
                            local[1] += JxW_face[qp]*s*dt*Us;

                        }


                    }
                }

        }

    }

    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    this->mass_dep  += global[1];
    this->total_mass = global[0] + this->mass_dep;

    if (libMesh::global_processor_id() == 0) {

        std::ofstream fout1;

        if (system.time == 0) {
            this->mass_dep = 0.0;
            this->init_mass = global[0] + global[1];
            fout1.open(fname);
            if (system.time == 0) fout1 << "#time        M_susp       M_dep        M_tot        M_tot/M_init" << endl;
        } else
            fout1.open(fname, ios_base::app);

        fout1 << std::left << std::scientific <<
                std::setw(8) << system.time << " " <<
                std::setw(8) << global[0] << " " <<
                std::setw(8) << global[1] << " " <<
                std::setw(8) << this->total_mass << " " <<
                std::setw(8) << this->total_mass / this->init_mass << endl;
        
        fout1.close();

    }

    return;
}