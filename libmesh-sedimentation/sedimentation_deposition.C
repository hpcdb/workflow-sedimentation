#include "sedimentation_deposition.h"

#include "define.h"


void SedimentationDeposition::setup()
{

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & sediment_system =
     es.get_system<TransientLinearImplicitSystem> ("sediment");
    
    ExplicitSystem & deposition_system = es.add_system<ExplicitSystem>("deposition"); 
    deposition_system.add_variable("d");
    deposition_system.add_vector("deposition_rate");
      
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
  
  NumericVector<Number> &deposition_rate = deposition_system .get_vector("deposition_rate");
  
 
  // Numeric ids corresponding to each variable in the system
  const unsigned int s_var  = sediment_system.variable_number ("s");
  const unsigned int d_var  = deposition_system.variable_number ("d");
 
  const Real c_factor = es.parameters.get<Real>("c_factor");
  const Real Us       = es.parameters.get<Real>("Us");
  const Real dt       = es.parameters.get<Real>("dt");
  
  
   
  
  const MeshBase::const_node_iterator nd_end = mesh.active_nodes_end();
  for (MeshBase::const_node_iterator nd      = mesh.active_nodes_begin(); nd != nd_end; ++nd) {
       const Node *node = *nd;
       std::vector<short int> boundary_ids = mesh.boundary_info->boundary_ids(node);
       if (boundary_ids.size() == 0 ) continue;
       if( std::find(boundary_ids.begin(), boundary_ids.end(),BOUNDARY_DEPOSITION) != boundary_ids.end()) 
       {
         
          unsigned int source_dof = node->dof_number(sediment_system.number(),s_var, 0);
          unsigned int dep_dof = node->dof_number(deposition_system.number(),d_var, 0);
          
          Number sed   = sediment_system.current_local_solution->el(source_dof);
          Number dep   = deposition_system.current_local_solution->el(dep_dof);
          Number value = sed*Us*dt*c_factor;
          
          deposition_rate.set(dep_dof, value);
          
  
       }
  }
  
  deposition_rate.close();
  deposition_system.solution->add(deposition_rate);
  //deposition_system.update();
  deposition_system.update();
  
}