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

void SedimentationTransport::PrintMass(std::ofstream& fmass, unsigned int t_step ) {

    double local[4] = {0.0, 0.0, 0.0, 0.0};
    double global[4] = {0.0, 0.0, 0.0, 0.0};

    const MeshBase & mesh = es.get_mesh();

    const unsigned int dim = mesh.mesh_dimension();
    TransientLinearImplicitSystem & transport_system =
            es.get_system<TransientLinearImplicitSystem> ("transport");   

    const unsigned int s_var = transport_system.variable_number("s");
    std::vector<dof_id_type> dof_indices_s;
    const DofMap & dof_map_transp = transport_system.get_dof_map();
    
    // Get a reference to the flow system object.
    TransientLinearImplicitSystem & flow_system =
            es.get_system<TransientLinearImplicitSystem> ("flow");
    
    // Numeric ids corresponding to each variable in the system
    const unsigned int u_var = flow_system.variable_number("u");
    const unsigned int v_var = flow_system.variable_number("v");
    unsigned int w_var;
    if (this->dim==3)
        w_var = flow_system.variable_number("w");    

    std::vector<dof_id_type> dof_indices_u;
    std::vector<dof_id_type> dof_indices_v;
    std::vector<dof_id_type> dof_indices_w;
    const DofMap & dof_map_flow = flow_system.get_dof_map();    
    
    std::vector<dof_id_type> dof_indices_face;

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = transport_system.variable_type(s_var);

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
    
    RealVectorValue inlet_velocity, outlet_velocity;

    for (; it != end; ++it) {
        
        const Elem *elem = *it;

        dof_map_transp.dof_indices(elem, dof_indices_s, s_var);
        dof_map_flow.dof_indices(elem, dof_indices_u, u_var);
        dof_map_flow.dof_indices(elem, dof_indices_v, v_var);
        if (this->dim==3)
            dof_map_flow.dof_indices(elem, dof_indices_w, w_var);        

        fe->reinit(elem);

        // Computing total suspended mass
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
            // Values to hold the current solution
            Number s = 0.0;
            // Compute the current solution
            for (unsigned int l = 0; l < phi.size(); l++) {
                s += phi[l][qp] * transport_system.current_solution(dof_indices_s[l]);
            }
            // total suspended mass
            local[0] += JxW[qp] * s;
        }

        // computing mass at boundaries: inlet and deposition for while and also discountint the "suspended" at deposition wall
        for (unsigned int sd = 0; sd < elem->n_sides(); sd++) {
            if (elem->neighbor(sd) == NULL) {

                fe_face->reinit(elem, sd);
                UniquePtr<Elem> side (elem->build_side(sd));

                // Applying flux advective boundary condition at deposition wall
                if (mesh.boundary_info->boundary_id(elem, sd) == this->deposition_id) {

                    for (unsigned int qp = 0; qp < qface.n_points(); qp++) {
                        Number s = 0.0;
                        for (unsigned int i = 0; i < phi_face.size(); i++) {
                            s += phi_face[i][qp] * transport_system.current_solution(dof_indices_s[i]);
                        }
                           
                        // mass leaving the domain through the deposition wall with velocity component Us (vertical direction)
                        local[1] += JxW_face[qp] * s * dt * Us;
                    }
                    
                    // Computing mass at bottom wall
                    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
                        // Values to hold the current solution
                        Number s = 0.0;                    
                        for (unsigned int sideNode=0; sideNode < side->n_nodes(); sideNode++) {                
                            for (unsigned int elemNode=0; elemNode<elem->n_nodes(); elemNode++) {
                                if (elem->node(elemNode) == side->node(sideNode)) {
                                    s += phi[elemNode][qp] * transport_system.current_solution(dof_indices_s[elemNode]);
                                }
                            }
                        }
                        // Discounting mass at deposition wall from total suspended mass
                        local[0] -= JxW[qp] * s;
                    }
                } // end if(deposition_id)
                    
                // Computing mass entering the domain through the "inlet" wall
                if (mesh.boundary_info->boundary_id(elem, sd) == this->inlet_bc_id) {
                    //double elem_v_x = 0.0, elem_v_y = 0.0;
                    // normal to the element face
                    const std::vector<Point> normal = fe_face->get_normals();
                    Number s = 0.0, vel_x = 0.0, vel_y = 0.0, vel_z = 0.0, vel_normal = 0.0;
                    for (unsigned int qp = 0; qp < qface.n_points(); qp++) {
                        s = 0.0, vel_x = 0.0, vel_y = 0.0, vel_z = 0.0, vel_normal = 0.0;
                        
                        for (unsigned int i = 0; i < phi_face.size(); i++) {
                            s += phi_face[i][qp] * transport_system.current_solution(dof_indices_s[i]);
                            vel_x += phi_face[i][qp] * flow_system.current_solution(dof_indices_u[i]);
                            vel_y += phi_face[i][qp] * flow_system.current_solution(dof_indices_v[i]);
                            if (this->dim==3)
                                vel_z += phi_face[i][qp] * flow_system.current_solution(dof_indices_w[i]);
                        }
                        //elem_v_x+=vel_x; elem_v_y+= vel_y;
                        inlet_velocity(0) = vel_x;
                        inlet_velocity(1) = vel_y;
                        if (this->dim==3)
                            inlet_velocity(2) = vel_z;
                        vel_normal = (inlet_velocity * normal[qp]);                        
                        // mass entering the domain through the inlet wall                     
                        local[2] -= JxW_face[qp] * s * vel_normal * dt; // to consider inlet, negative sign must be
                    }
                    //cout<<" iel # = "<<elem->id()<<" vx = "<< elem_v_x/qface.n_points()<< "  vy = "<<elem_v_y/qface.n_points()<<" s = "<< s<<endl;
                } //end if (inlet_BC)
                
                // Computing volume of fluid leaving the domain through the "outlet" wall
                if (mesh.boundary_info->boundary_id(elem, sd) == es.parameters.set<int> ("outflow_id")) {
                     
                    double elem_v_x = 0.0, elem_v_y = 0.0;
                    // normal to the element face
                    const std::vector<Point> normal = fe_face->get_normals();
                    Number vel_x = 0.0, vel_y = 0.0, vel_z = 0.0, vel_normal = 0.0;
                    for (unsigned int qp = 0; qp < qface.n_points(); qp++) {
                        vel_x = 0.0, vel_y = 0.0, vel_z = 0.0, vel_normal = 0.0;
                        for (unsigned int i = 0; i < phi_face.size(); i++) {
                            vel_x += phi_face[i][qp] * flow_system.current_solution(dof_indices_u[i]);
                            vel_y += phi_face[i][qp] * flow_system.current_solution(dof_indices_v[i]);
                            if (this->dim==3)
                                vel_z += phi_face[i][qp] * flow_system.current_solution(dof_indices_w[i]);
                        }
                        elem_v_x+=vel_x; elem_v_y+= vel_y;
                        outlet_velocity(0) = vel_x;
                        outlet_velocity(1) = vel_y;
                        if (this->dim==3)
                            outlet_velocity(2) = vel_z;
                        vel_normal = (outlet_velocity * normal[qp]);
                        // Volume of fluid leaving the domain through the outlet wall
                        local[3] += JxW_face[qp] * vel_normal * dt;
                    }
                    //cout<<" iel # = "<<elem->id()<<" vx = "<< elem_v_x/qface.n_points()<< "  vy = "<<elem_v_y/qface.n_points()<<endl;
                } //end if (outlet_BC)
            } //end if (elem->neighbor(sd) == NULL)
        } // end loop over element's sides
    } // end element iterator

    MPI_Allreduce(local, global, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    this->mass_dep  += global[1];
    this->total_mass = global[0] + this->mass_dep;
    this->inlet_mass += global[2];
    this->outlet_volume += global[3];

    if (libMesh::global_processor_id() == 0) {

        std::ofstream fout1;

        if (transport_system.time == 0) {
            this->mass_dep = 0.0;
            this->total_mass = this->init_mass = global[0];
            if (this->init_mass>0.0) {
                fmass << " initial Mass = "<< std::setw(8) << this->init_mass<<endl;
                fmass << "    time        M_dep      M_dep_tot      M_susp       M_tot    M_tot/M_init" << endl;
            } else {
                fmass << "    time       M_inlet    M_inlet_tot     M_susp       M_dep      M_dep_tot      M_tot    M_tot/M_in_tot   V_out      V_tot_out" << endl;
            }
        } else {
            if ((t_step) % es.parameters.get<unsigned int>("write_interval") == 0 || abs( transport_system.time - es.parameters.get<Real>("tmax") )<1.0e-8 ) {
                if (this->init_mass>0.0) {
                    fmass << std::left << std::scientific <<
                        std::setw(8) << transport_system.time << " " <<
                        std::setw(8) << global[1] << " " <<
                        std::setw(8) << this->mass_dep << " " <<
                        std::setw(8) << global[0] << " " <<
                        std::setw(8) << this->total_mass << " " <<
                        std::setw(8) << this->total_mass / this->init_mass << endl;
                } else {
                        fmass << std::left << std::scientific <<
                        std::setw(8) << transport_system.time << " " <<
                        std::setw(8) << global[2] << " " <<
                        std::setw(8) << this->inlet_mass << " " <<
                        std::setw(8) << global[0] << " " <<
                        std::setw(8) << global[1] << " " <<
                        std::setw(8) << this->mass_dep << " " <<
                        std::setw(8) << this->total_mass << " " <<
                        std::setw(8) << this->total_mass / this->inlet_mass << " " <<
                        std::setw(8) << global[3] << " " <<
                        std::setw(8) << this->outlet_volume << " " <<  endl;
                }
            }
        }
    }

    return;
}