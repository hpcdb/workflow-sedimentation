
#include "libmesh/fe_interface.h"

#include "sedimentation_deposition.h"

#include "define.h"


void SedimentationDeposition::init()
{
    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & sediment_system =
     es.get_system<TransientLinearImplicitSystem> ("sediment");
    
    ExplicitSystem & deposition_system = es.add_system<ExplicitSystem>("deposition"); 
    deposition_system.add_variable("d");
    deposition_system.add_vector("deposition_rate");
    
    //ExplicitSystem & deposition_rate  = es.add_system<ExplicitSystem>("deposition_rate"); 
    //deposition_rate.add_variable("dVdt");
    
}


void SedimentationDeposition::setup(GetPot &infile)
{
    ExplicitSystem & deposition_system = es.get_system<ExplicitSystem>("deposition"); 
    deposition_system.add_vector("deposition_rate");
    this->deposition_id = infile("dirichlet/deposition", -1);
}


void SedimentationDeposition::print()
{
    
    TransientLinearImplicitSystem & sediment_system =
     es.get_system<TransientLinearImplicitSystem> ("sediment");
    
    ExplicitSystem & deposition_system = es.get_system<ExplicitSystem>("deposition"); 
    
    Real Max  = deposition_system.solution->max();
    Real Min  = deposition_system.solution->min();
    Real norm = deposition_system.solution->l2_norm(); 
    
    std:: cout << " Volume deposited statistics: " << std::endl;
    std:: cout << " Vol_{max}: " << Max << std::endl;
    std:: cout << " Vol_{min}: " << Min << std::endl;
    std:: cout << " Vol_{l2} : " << norm << std::endl;
          
}



void SedimentationDeposition::ComputeDeposition()
{
  // Get a constant reference to the mesh object.
  const MeshBase& mesh    = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Convection-Diffusion system object.
  TransientLinearImplicitSystem & sediment_system =
     es.get_system<TransientLinearImplicitSystem> ("sediment");
  
  ExplicitSystem & deposition_system = es.get_system<ExplicitSystem>("deposition"); 
  
  //ExplicitSystem &       deposition_rate = es.get_system<ExplicitSystem>("deposition_rate"); 
  NumericVector<Number> &deposition_rate = deposition_system .get_vector("deposition_rate");
  
  const DofMap& dof_map      = sediment_system.get_dof_map();
  std::vector<dof_id_type> dof_indices_face; 
  
  // Numeric ids corresponding to each variable in the system
  const unsigned int s_var  = sediment_system.variable_number ("s");
  const unsigned int d_var  = deposition_system.variable_number ("d");
  //const unsigned int dVdt_var  = deposition_rate.variable_number ("dVdt");
  FEType fe_type = sediment_system.variable_type(s_var);
  
  UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));

  // A Gauss quadrature rule for numerical integration.
  // Let the \p FEType object decide what order rule is appropriate.
  QGauss qface (dim-1, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe_face->attach_quadrature_rule (&qface);
 
  const Real c_factor = es.parameters.get<Real>("c_factor");
  const Real Us       = es.parameters.get<Real>("Us");
  const Real dt       = es.parameters.get<Real>("dt");
  
  
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  NumericVector<Number> & sys_soln(*sediment_system.current_local_solution);
  std::vector<Number>   elem_soln;
  std::vector<Number>   nodal_soln;
  for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem * elem = *el;
      
      for (unsigned int s=0; s<elem->n_sides(); s++)
          if (elem->neighbor(s) == libmesh_nullptr)
            {
              if(mesh.get_boundary_info().has_boundary_id(elem, s, deposition_id))
              {
                  
                 UniquePtr<const Elem> side (elem->build_side(s));
                 
                 elem_soln.resize(side->n_nodes());
                 
                 dof_map.dof_indices(side.get(),dof_indices_face);
                 
                 for(int i =0; i < dof_indices_face.size(); i++)
                     elem_soln[i] = sys_soln(dof_indices_face[i]);
                 
                 FEInterface::nodal_soln(dim,fe_type,side.get(),elem_soln, nodal_soln);
                 
                 libmesh_assert_equal_to (nodal_soln.size(), side->n_nodes());
                 libmesh_assert_equal_to (nodal_soln.size(), dof_indices_face.size());
                 
                 for(int i =0; i < side->n_nodes(); i++)
                 {
                     const Node *node = side->get_node(i);
                     unsigned int source_dof = node->dof_number(sediment_system.number(),s_var, 0);
                     unsigned int dep_dof    = node->dof_number(deposition_system.number(),d_var, 0);
                     Number value = nodal_soln[i]*Us*dt*c_factor;
                     deposition_rate.set(dep_dof, value);
                     
                 }
                 
              }
            }
    }
  

  
  deposition_rate.close();
  deposition_system.solution->add(deposition_rate);
  deposition_system.update();
  

  
}
