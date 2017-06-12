
#include "libmesh/fe_interface.h"

#include "sedimentation_deposition.h"

#include "define.h"

void SedimentationDeposition::init() {
    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & sediment_system =
            es.get_system<TransientLinearImplicitSystem> ("transport");

    ExplicitSystem & deposition_system = es.add_system<ExplicitSystem>("deposition");
    ExplicitSystem & deposition_rate   = es.add_system<ExplicitSystem>("deposition rate");
    //ExplicitSystem & bed_load = es.add_system<ExplicitSystem>("bed load");

    deposition_system.add_variable("d");
    deposition_rate.add_variable("r");
    //bed_load.add_variable("q-x");
    //bed_load.add_variable("q-y");
    //if (es.get_mesh().mesh_dimension() == 3)
    //    bed_load.add_variable("q-z");

}


void SedimentationDeposition::setup(GetPot &infile)
{
    ExplicitSystem & deposition_system = es.get_system<ExplicitSystem>("deposition"); 
    deposition_system.add_vector("deposition_rate");
    this->deposition_id = infile("transport/deposition", -1);
}


void SedimentationDeposition::print()
{
    
    TransientLinearImplicitSystem & sediment_system =
            es.get_system<TransientLinearImplicitSystem> ("transport");

    ExplicitSystem & deposition_system = es.get_system<ExplicitSystem>("deposition");

    Real Max = deposition_system.solution->max();
    Real Min = deposition_system.solution->min();
    Real norm = deposition_system.solution->l2_norm();

    std::cout << " Volume deposited statistics: " << std::endl;
    std::cout << " Vol_{max}: " << Max << std::endl;
    std::cout << " Vol_{min}: " << Min << std::endl;
    std::cout << " Vol_{l2} : " << norm << std::endl;

}

void SedimentationDeposition::ComputeDeposition() {
    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & sediment_system =
            es.get_system<TransientLinearImplicitSystem> ("transport");

    ExplicitSystem & deposition_system = es.get_system<ExplicitSystem>("deposition");

    ExplicitSystem & deposition_rate = es.get_system<ExplicitSystem>("deposition rate");


    const DofMap& dof_map = sediment_system.get_dof_map();
    std::vector<dof_id_type> dof_indices_face;

    // Numeric ids corresponding to each variable in the system
    const unsigned int s_var = sediment_system.variable_number("s");
    const unsigned int d_var = deposition_system.variable_number("d");
    const unsigned int r_var = deposition_rate.variable_number("r");
    FEType fe_type = sediment_system.variable_type(s_var);

    UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_type));

    // A Gauss quadrature rule for numerical integration.
    // Let the \p FEType object decide what order rule is appropriate.
    QGauss qface(dim - 1, fe_type.default_quadrature_order());

    // Tell the finite element object to use our quadrature rule.
    fe_face->attach_quadrature_rule(&qface);

    const Real c_factor = es.parameters.get<Real>("c_factor");
    const Real Us = es.parameters.get<Real>("Us");
    const Real dt = es.parameters.get<Real>("dt");


    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    NumericVector<Number> & sys_soln(*sediment_system.current_local_solution);
    std::vector<Number> elem_soln;
    std::vector<Number> nodal_soln;
    
    

    for (; el != end_el; ++el) {
        // Store a pointer to the element we are currently
        // working on.  This allows for nicer syntax later.
        const Elem * elem = *el;

        for (unsigned int s = 0; s < elem->n_sides(); s++)
            if (elem->neighbor(s) == libmesh_nullptr) {
                if (mesh.get_boundary_info().has_boundary_id(elem, s, deposition_id)) {

                    UniquePtr<const Elem> side(elem->build_side(s));

                    elem_soln.resize(side->n_nodes());

                    dof_map.dof_indices(side.get(), dof_indices_face);

                    for (int i = 0; i < dof_indices_face.size(); i++)
                        elem_soln[i] = sys_soln(dof_indices_face[i]);

                    FEInterface::nodal_soln(dim, fe_type, side.get(), elem_soln, nodal_soln);

                    libmesh_assert_equal_to(nodal_soln.size(), side->n_nodes());
                    libmesh_assert_equal_to(nodal_soln.size(), dof_indices_face.size());

                    for (int i = 0; i < side->n_nodes(); i++) {
                        const Node *node = side->get_node(i);
                        unsigned int source_dof = node->dof_number(sediment_system.number(), s_var, 0);
                        //unsigned int dep_dof = node->dof_number(deposition_system.number(), d_var, 0);
                        unsigned int r_dof = node->dof_number(deposition_rate.number(), r_var, 0);
                        Number value = nodal_soln[i] * Us * dt*c_factor;
                        deposition_rate.solution->set(r_dof, value);
                        //deposition_system.solution->add(dep_dof, value);


                    }

                }
            }
    }


    deposition_rate.solution->close();
    
    deposition_system.solution->add(*deposition_rate.solution);
    
    deposition_system.solution->close();

    deposition_rate.update();
    deposition_system.update();

}

#if 0 

void SedimentationDeposition::ComputeBedLoad() {
    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & flow_system =
            es.get_system<TransientLinearImplicitSystem> ("flow");

    ExplicitSystem & bedload_system = es.get_system<ExplicitSystem>("bed load");

    const DofMap& dof_map_flow = flow_system.get_dof_map();
    const DofMap& dof_map_bedload = bedload_system.get_dof_map();

    std::vector<dof_id_type> dof_indices_face;

    // Numeric ids corresponding to each variable in the system
    const unsigned int u_var = flow_system.variable_number("u");
    const unsigned int v_var = flow_system.variable_number("v");
    unsigned w_var = 0;
    if (dim > 2)
        w_var = flow_system.variable_number("w");

    // Numeric ids corresponding to each variable in the system
    const unsigned int qx_var = bedload_system.variable_number("q-x");
    const unsigned int qy_var = bedload_system.variable_number("q-y");
    unsigned qz_var = 0;
    if (dim > 2)
        qz_var = bedload_system.variable_number("qz");


    FEType fe_type = flow_system.variable_type(u_var);

    UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_type));

    // A Gauss quadrature rule for numerical integration.
    // Let the \p FEType object decide what order rule is appropriate.
    QGauss qface(dim - 1, fe_type.default_quadrature_order());



    // Tell the finite element object to use our quadrature rule.
    fe_face->attach_quadrature_rule(&qface);

    const std::vector<Real>& JxW_face = fe_face->get_JxW();

    const std::vector<std::vector<Real> >& phi_face = fe_face->get_phi();
    const std::vector<Point> & normals = fe_face->get_normals();

    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<RealGradient> >& dphi_face = fe_face->get_dphi();

    const Real c_factor = es.parameters.get<Real>("c_factor");
    const Real Us = es.parameters.get<Real>("Us");
    const Real dt = es.parameters.get<Real>("dt");

    std::vector<dof_id_type> dof_indices_u;
    std::vector<dof_id_type> dof_indices_v;
    std::vector<dof_id_type> dof_indices_w;

    std::vector<dof_id_type> dof_indices_qx;
    std::vector<dof_id_type> dof_indices_qy;
    std::vector<dof_id_type> dof_indices_qz;


    Real shields_c = 0.05;
    Real rhoL = es.parameters.get<Real>("rhoL");
    Real rhoS = es.parameters.get<Real>("rhoS");
    Real g = es.parameters.get<Real>("g");
    Real ds = es.parameters.get<Real>("ds");



    for (; el != end_el; ++el) {
        // Store a pointer to the element we are currently
        // working on.  This allows for nicer syntax later.
        const Elem * elem = *el;

        dof_map_flow.dof_indices(elem, dof_indices_u, u_var);
        dof_map_flow.dof_indices(elem, dof_indices_v, v_var);
        if (dim > 2) dof_map_flow.dof_indices(elem, dof_indices_w, w_var);
        
        
        for (unsigned int s = 0; s < elem->n_sides(); s++)
            if (elem->neighbor(s) == libmesh_nullptr) {
                if (mesh.get_boundary_info().has_boundary_id(elem, s, deposition_id)) {

                    fe_face->reinit(elem, s);

                    for (unsigned int qp = 0; qp < qface.n_points(); qp++) {

                        Gradient grad_u, grad_v, grad_w;

                        for (unsigned int i = 0; i < phi_face.size(); i++) {

                            grad_u.add_scaled(dphi_face[i][qp], flow_system.current_solution(dof_indices_u[i]));
                            grad_v.add_scaled(dphi_face[i][qp], flow_system.current_solution(dof_indices_v[i]));
                            if (dim > 2) grad_w.add_scaled(dphi_face[i][qp], flow_system.current_solution(dof_indices_w[i]));
                        }



                    }

                    // Bed shear stress (dimensional form (*ub/Lb))
                    RealTensor G;
                    G(0, 0) = visc * (2.0 * grad_u(0)) * ub / Lb;
                    G(0, 1) = G(1, 0) = visc * (grad_u(1) + grad_v(0)) * ub / Lb;
                    G(1, 1) = visc * 2.0 * grad_v(1) * ub / Lb;
                    RealVectorValue taub;
                    taub(0) = normals[qp](0) * G(0, 0) + normals[qp](1) * G(1, 0);
                    taub(1) = normals[qp](0) * G(0, 1) + normals[qp](1) * G(1, 1);
                    Number taub_norm = taub.norm();


                    // Shields number
                    Real R = (rhoS - rhoL) / rhoL;
                    Real shields = taub_norm / (rhoL * g * R * ds);

                    // Bed load rate
                    Real qstar = 0.0;
                    if (shields > shields_c) {
                        qstar = 18.74 * (shields - shields_c)*(std::sqrt(shields) - 0.7 * std::sqrt(shields_c));
                    }

                    Real q0 = qstar * std::sqrt(R * g * ds) * ds;

                    RealVectorValue qb;
                    qb(0) = q0 * taub(0) / taub_norm;
                    qb(1) = q0 * taub(1) / taub_norm;


                    UniquePtr<const Elem> side(elem->build_side(s));

                    dof_map_bedload.dof_indices(side.get(), dof_indices_qx, qx_var);
                    dof_map_bedload.dof_indices(side.get(), dof_indices_qy, qy_var);
                    if (dim > 2) dof_map_bedload.dof_indices(side.get(), dof_indices_qz, qz_var);

                        
                    int tmp = (dim > 2)? dof_indices_qz.size(): 0; 
                    elem_soln.resize(dof_indices_qx.size + dof_indices_qy.size + tmp);

                    for (int i = 0; i < dof_indices_qx.size(); i++)
                        elem_soln[i] = sys_soln(dof_indices_face[i]);

                    FEInterface::nodal_soln(dim, fe_type, side.get(), elem_soln, nodal_soln);

                    libmesh_assert_equal_to(nodal_soln.size(), side->n_nodes());
                    libmesh_assert_equal_to(nodal_soln.size(), dof_indices_face.size());

                    for (int i = 0; i < side->n_nodes(); i++) {
                        const Node *node = side->get_node(i);
                        unsigned int source_dof = node->dof_number(sediment_system.number(), s_var, 0);
                        unsigned int dep_dof = node->dof_number(deposition_system.number(), d_var, 0);
                        unsigned int r_dof = node->dof_number(deposition_rate.number(), r_var, 0);
                        Number value = nodal_soln[i] * Us * dt*c_factor;
                        deposition_rate.solution->set(dep_dof, value);
                        deposition_system.solution->add(r_dof, value);


                    }

                }
            }
    }



    deposition_rate.solution->close();
    deposition_system.solution->close();

    deposition_rate.update();
    deposition_system.update();
   
    


}

#endif