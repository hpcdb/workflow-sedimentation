
// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>
using namespace std;

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/vtk_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/perf_log.h"
#include "libmesh/boundary_info.h"
#include "libmesh/utility.h"

#include "libmesh/getpot.h"
#include "libmesh/mesh_refinement.h"
// For systems of equations the \p DenseSubMatrix
// and \p DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/numeric_vector.h"


// To impose Dirichlet boundary conditions
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/analytic_function.h"

// The definition of a geometric element
#include "libmesh/elem.h"
#include "define.h"
using namespace libMesh;

#include "mesh_moviment.h"
#include "define.h"

#include "libmesh/zero_function.h"


void MeshMoviment::setup()
{

  LinearImplicitSystem & mesh_system = this->es.add_system<LinearImplicitSystem> ("mesh_moviment");
  mesh_system.add_vector("mesh-velocity");
  unsigned int dispz = mesh_system.add_variable ("disp-z");

  mesh_system.attach_assemble_object(*this);

}




void MeshMoviment::assemble()
{


  //std::vector<Real> tau_fsi;

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Convection-Diffusion system object.
  LinearImplicitSystem & system =
    es.get_system<LinearImplicitSystem> ("mesh_moviment");

  ExplicitSystem &      deposition_system = es.get_system<ExplicitSystem>("deposition");
  NumericVector<Number> &deposition_rate = deposition_system .get_vector("deposition_rate");


  // Numeric ids corresponding to each variable in the system
  const unsigned int dispz_var     = system.variable_number ("disp-z");
  const unsigned int d_var         = deposition_system.variable_number("d");

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = system.variable_type(dispz_var);

  // Build a Finite Element object of the specified type.  Since the
  // \p FEBase::build() member dynamically creates memory we will
  // store the object as an \p UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> fe      (FEBase::build(dim, fe_type));


  // A Gauss quadrature rule for numerical integration.
  // Let the \p FEType object decide what order rule is appropriate.
  QGauss qrule (dim,   fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule      (&qrule);

  UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
  QGauss qface(dim-1, FIFTH);
  fe_face->attach_quadrature_rule (&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.  We will start
  // with the element Jacobian * quadrature weight at each integration point.
  const std::vector<Real>& JxW      = fe->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi      = fe->get_phi();


  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> >& dphi      = fe->get_dphi();

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the \p DofMap
  // in future examples.
  const DofMap& dof_map      = system.get_dof_map();
  const DofMap& dof_map_dep  = deposition_system.get_dof_map();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_deposition;

  Number aux1, aux2;

  {

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    double he;
    double vmaxt = -1.0E8;
    double vmint = 1.0E8;
    double vmax, vmin;
    int iel = 0;

    for ( ; el != end_el; ++el)
    {
        // Store a pointer to the element we are currently
        // working on.  This allows for nicer syntax later.
        const Elem* elem = *el;
        he = elem->volume();
        vmaxt = (he>vmaxt)?he:vmaxt;
        vmint = (he<vmint)?he:vmint;
        iel++;
    }

    MPI_Allreduce(&vmaxt, &vmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&vmint, &vmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    cout << "||dV/dt|| = " << deposition_rate.l2_norm() << endl;
    cout << "Vol_{min}: " << vmin << endl;
    cout << "Vol_{max}: " << vmax << endl;

    aux1 = vmin/vmax;
    aux2 = 1.0/vmax;

  }

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);
      dof_map_dep.dof_indices (elem, dof_indices_deposition, d_var);

      // Compute the element-specific data for the current
      // element.
      fe->reinit (elem);

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Ke.resize (dof_indices.size(), dof_indices.size());
      Fe.resize (dof_indices.size());

      Number tau_fsi = elem->volume();
      Number aux3 = tau_fsi*aux2;
      tau_fsi  = (1.0 - aux1)/aux3;

      // loop over quadrature points
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {
            // Now compute the element matrix and RHS contributions.
            for (unsigned int i=0; i<phi.size(); i++) {
                for (unsigned int j=0; j<phi.size(); j++) {
                    // The Galerkin contribution
                    Ke(i,j) += JxW[qp]*(tau_fsi)*(dphi[i][qp]*dphi[j][qp]);
                }
            }
      }


      {

       // The penalty value.
       const Real penalty = 1.e10;

       for (unsigned int s=0; s<elem->n_sides(); s++)
       if (elem->neighbor(s) == NULL)
       {
              UniquePtr<Elem> side = elem->side(s);


              if( mesh.boundary_info->boundary_id(elem,s) == BOUNDARY_DEPOSITION)
              {

                const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
                const std::vector<Real>& JxW_face                = fe_face->get_JxW();
                fe_face->reinit(elem, s);

                for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {

                    const Real penalty = 1.e10;

                    const Real value = 0.001;
                  
                    for (unsigned int i=0; i<phi_face.size(); i++)
                      for (unsigned int j=0; j<phi_face.size(); j++)
                        Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];

                    for (unsigned int i=0; i<phi_face.size(); i++)
                      Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
                 }

                /*
                    // Loop over the nodes on the side.
                    for (unsigned int ns=0; ns<side->n_nodes(); ns++)
                    {
                        for (unsigned int n=0; n<elem->n_nodes(); n++)
                            if (elem->node(n) == side->node(ns))
                            {

                                int dof = elem->get_node(n)->dof_number(deposition_system.number(),0,0);
                                //Real value = deposition_rate(dof_indices_deposition[n]);
                                Real value = deposition_rate.el(dof);                                // Matrix contribution.
                                Ke(n,n) += penalty;
                                // Right-hand-side contribution.
                                Fe(n)   += penalty*value;
                            }

                    }
                 * */

               }


            }

         }

        // If this assembly program were to be used on an adaptive mesh,
        // we would have to apply any hanging node constraint equations
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The \p SparseMatrix::add_matrix()
        // and \p NumericVector::add_vector() members do this for us.
        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
    }


}

    // That concludes the system matrix assembly routine.

void MeshMoviment::updateMesh()
{
    // Get a constant reference to the mesh object.
  MeshBase& mesh =  es.get_mesh();
  const double dt = es.parameters.get<Real>("dt");

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Convection-Diffusion system object.
  LinearImplicitSystem & system =
    es.get_system<LinearImplicitSystem> ("mesh_moviment");

  NumericVector<Number> & mesh_velocity = system.get_vector("mesh-velocity");

  // Loop over all nodes and copy the location from the current system to
  // the auxiliary system.
  const MeshBase::const_node_iterator nd_end = mesh.local_nodes_end();
  for (MeshBase::const_node_iterator nd = mesh.local_nodes_begin(); nd != nd_end; ++nd)
  {
      Node *node = *nd;
      unsigned int dof      = node->dof_number(system.number(),0, 0);

      Number value = system.current_local_solution->el(dof);
      (*node)(2) += value;
      mesh_velocity.set(dof, value/dt);

  }


}
