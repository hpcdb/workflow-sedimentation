
// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>
using namespace std;

#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel_ghost_sync.h"

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/vtk_io.h"
#include "libmesh/gmsh_io.h"
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


#include "libmesh/fe_interface.h"
// To impose Dirichlet boundary conditions
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/analytic_function.h"

// The definition of a geometric element
#include "libmesh/elem.h"
#include "define.h"


#include "mesh_moviment.h"
#include "define.h"

#include "libmesh/zero_function.h"

using namespace libMesh;


void MeshMoviment::init()
{
    LinearImplicitSystem & mesh_system = this->es.add_system<LinearImplicitSystem> ("MeshMoving");
    unsigned int dispz = mesh_system.add_variable ("disp-z");
}


void MeshMoviment::setup(GetPot &infile)
{

  LinearImplicitSystem & mesh_system = this->es.get_system<LinearImplicitSystem> ("MeshMoving");
   unsigned int dispz = mesh_system.variable_number ("disp-z");
  mesh_system.attach_assemble_object(*this);

  this->deposition_id = infile("dirichlet/deposition", 0);
  this->fixedwall_id  = infile("dirichlet/fixedwall", 0);  
  
  std::set<boundary_id_type> fixed;
  std::vector<unsigned int> vars(1);
  vars[0] = dispz;
  ZeroFunction<Number> zero;
  
  if(fixedwall_id > 0) {
      fixed.insert(fixedwall_id);
      mesh_system.get_dof_map().add_dirichlet_boundary(DirichletBoundary(fixed,vars, &zero));
  }
  

}


void MeshMoviment::assemble()
{
   
  PerfLog* perf_log = es.parameters.get<PerfLog*>("PerfLog");
  perf_log->pause_event("Solver","Mesh Moviment");
  perf_log->start_event("Assembly","Mesh Moviment");
  
  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Convection-Diffusion system object.
  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem> ("MeshMoving");

  
  //ExplicitSystem        &  deposition_system = es.get_system<ExplicitSystem>("deposition");
  ExplicitSystem        &  deposition_rate = es.get_system<ExplicitSystem>("deposition rate");
  //NumericVector<Number> &  deposition_rate   = deposition_system .get_vector("deposition_rate");


  // Numeric ids corresponding to each variable in the system
  const unsigned int dispz_var     = system.variable_number ("disp-z");
  const unsigned int r_var         = deposition_rate.variable_number("r");

  
  
  const Real c_factor = es.parameters.get<Real>("c_factor");
  const Real Us       = es.parameters.get<Real>("Us");
  const Real dt       = es.parameters.get<Real>("dt");
  
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
  
  
  const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
  const std::vector<Real>& JxW_face                = fe_face->get_JxW();

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the \p DofMap
  // in future examples.
  const DofMap& dof_map      = system.get_dof_map();
  const DofMap& dof_map_dep  = deposition_rate.get_dof_map();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_dep;
  std::vector<dof_id_type> dof_indices_dep_face;
  

 
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

    //cout << "||dV/dt|| = " << deposition_rate.l2_norm() << endl;
    //cout << "Vol_{min}: " << vmin << endl;
    //cout << "Vol_{max}: " << vmax << endl;

    aux1 = vmin/vmax;
    aux2 = 1.0/vmax;

  }

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
   NumericVector<Number> & sys_soln(*deposition_rate.current_local_solution);
   
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

      //std::cout << " DoF Size = " << dof_indices.size() <<  std::endl;
      
      {

       // The penalty value.
       const Real penalty = 1.e10;
       std::vector<Number>   elem_soln;
       std::vector<Number>   nodal_soln;

       for (unsigned int s=0; s<elem->n_sides(); s++)
       if (elem->neighbor(s) == NULL)
       {
              UniquePtr<Elem> side = elem->side(s);

              if( mesh.boundary_info->boundary_id(elem,s) == deposition_id)
              {
                 elem_soln.resize(side->n_nodes());
                 
                 dof_map_dep.dof_indices(side.get(),dof_indices_dep_face);
                 
                 for(int i =0; i < dof_indices_dep_face.size(); i++)
                     elem_soln[i] = sys_soln(dof_indices_dep_face[i]);
                 
                 FEInterface::nodal_soln(dim,fe_type,side.get(),elem_soln, nodal_soln);
                 
                 libmesh_assert(nodal_soln.size() == elem_soln.size());
                  
                 for(int n = 0; n < elem->n_nodes(); n++)
                     for(int i =0; i < side->n_nodes(); i++)
                     {
                         if(side->node(i) == elem->node(n))
                         {
                             Ke(n,n) += penalty;
                             Fe(n)   += nodal_soln[i]*penalty*c_factor;
                            
                         }
                    
                     }


               }


            }

         }

        // If this assembly program were to be used on an adaptive mesh,
        // we would have to apply any hanging node constraint equations
         dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The \p SparseMatrix::add_matrix()
        // and \p NumericVector::add_vector() members do this for us.
        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
    }
  
    perf_log->stop_event("Assembly","Mesh Moviment");
    perf_log->restart_event("Solver","Mesh Moviment");

}

    // That concludes the system matrix assembly routine.

void MeshMoviment::updateMesh()
{
    // Get a constant reference to the mesh object.
  MeshBase& mesh =  es.get_mesh();
  const double dt = es.parameters.get<Real>("dt");

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();
  int d = dim -1;
  // Get a reference to the Convection-Diffusion system object.
  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem> ("MeshMoving");

  // Loop over all nodes and copy the location from the current system to
  // the auxiliary system.
  const MeshBase::const_node_iterator nd_end = mesh.local_nodes_end();
  for (MeshBase::const_node_iterator nd = mesh.local_nodes_begin(); nd != nd_end; ++nd)
  {
      Node *node = *nd;
      unsigned int dof      = node->dof_number(system.number(),0, 0);

      Number value = system.current_local_solution->el(dof);
      (*node)(d) += value;

  }

  SyncNodalPositions sync_object(mesh);
  Parallel::sync_dofobject_data_by_id(mesh.comm(), mesh.nodes_begin(), mesh.nodes_end(), sync_object);

  
}
