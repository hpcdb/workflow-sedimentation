/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */




#include "viscosity_flow.h"

void ViscosityFlow::init() {

    LinearImplicitSystem &system = this->es.add_system<LinearImplicitSystem>("viscosity");
    system.add_variable("mu");

}


void ViscosityFlow::assemble() {

    const MeshBase& mesh = this->es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the LinearImplicitSystem we are solving
    LinearImplicitSystem & system = es.get_system<LinearImplicitSystem> ("viscosity");
    
    const unsigned int mu_var = system.variable_number("mu");
    
    // A reference to the  DofMap object for this system.  
    const DofMap & dof_map = system.get_dof_map();

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = dof_map.variable_type(mu_var);

    // Build a Finite Element object of the specified type.  
    UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));

    
    // A 5th order Gauss quadrature rule for numerical integration.
    QGauss qrule(dim, fe_type.default_quadrature_order());

    
    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);

    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> & JxW = fe->get_JxW();

    // The physical XY locations of the quadrature points on the element.
    // These might be useful for evaluating spatially varying material
    // properties at the quadrature points.
    const std::vector<Point> & q_point = fe->get_xyz();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real> > & phi = fe->get_phi();

   
    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  
    DenseMatrix<Number> Me;
    DenseVector<Number> Fe;

    
    // This vector will hold the degree of freedom indices for
    // the element.  
    std::vector<dof_id_type> dof_indices;

    // Now we will loop over all the elements in the mesh.
    // We will compute the element matrix and right-hand-side
    // contribution.
    //
  
    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    // Loop over the elements.  Note that  ++el is preferred to
    // el++ since the latter requires an unnecessary temporary
    // object.
    for (; el != end_el; ++el) {
        
        // Store a pointer to the element we are currently
        // working on.  This allows for nicer syntax later.
        const Elem * elem = *el;

        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        dof_map.dof_indices(elem, dof_indices);

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe->reinit(elem);

        // Zero the element matrix and right-hand side before
        // summing them.  

        Me.resize(dof_indices.size(), dof_indices.size());
        Fe.resize(dof_indices.size());

        // Now loop over the quadrature points.  This handles
        // the numeric integration.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

            // Now we will build the element matrix.  This involves
            // a double loop to integrate the test funcions (i) against
            // the trial functions (j).
            for (std::size_t i = 0; i < phi.size(); i++)
                for (std::size_t j = 0; j < phi.size(); j++) {
                    Me(i, j) += JxW[qp]*(phi[i][qp] * phi[j][qp]);
                }

            // This is the end of the matrix summation loop
            // Now we build the element right-hand-side contribution.
            // This involves a single loop in which we integrate the
            // "forcing function" in the PDE against the test functions.
            {
                const Real x = q_point[qp](0);
                const Real y = q_point[qp](1);
                const Real eps = 1.e-3;
                
                const Real mu_e = 1.0;

                for (std::size_t i = 0; i < phi.size(); i++)
                    Fe(i) += JxW[qp] * mu_e * phi[i][qp];
            }
        }
    }

}



