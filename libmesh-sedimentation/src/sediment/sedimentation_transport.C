

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
#include "libmesh/exodusII_io.h"

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


// Define a wrapper for exact_solution that will be needed below

Number init_value(const Point& p,
        const Parameters& parameters,
        const std::string&,
        const std::string&) {
    double xlock = parameters.get<Real> ("xlock");
    double hlock = parameters.get<Real> ("hlock");
    int d = parameters.get<int> ("dim");
    if (((p(0) - xlock) <= 1.0E-03) && ((p(d - 1) - hlock) <= 1.0E-03)) return 1.0;
    return 0.0;
}

void init_sedimentation(EquationSystems& es, const std::string& system_name) {

    libmesh_assert_equal_to(system_name, "transport");

    TransientLinearImplicitSystem & system =
            es.get_system<TransientLinearImplicitSystem>("transport");

    es.parameters.set<Real> ("time") = system.time = 0;

    system.project_solution(init_value, NULL, es.parameters);
}

void SedimentationTransport::init() {
    const MeshBase& mesh = es.get_mesh();

    // The dimension that we are running
    this->dim = mesh.mesh_dimension();

    TransientLinearImplicitSystem & transport_system = this->es.add_system<TransientLinearImplicitSystem> ("transport");

    unsigned int s_var = transport_system.add_variable("s");


}

void SedimentationTransport::setup(GetPot &infile) {
    
    const MeshBase& mesh = es.get_mesh();

    // The dimension that we are running
    this->dim = mesh.mesh_dimension();

    TransientLinearImplicitSystem & transport_system = this->es.get_system<TransientLinearImplicitSystem> ("transport");
    
    unsigned int s_var = transport_system.variable_number("s");

    transport_system.attach_assemble_object(*this);

    std::set<boundary_id_type> czero;
    std::set<boundary_id_type> cprescribed;

    std::vector<unsigned int> sed_vars(1);
    sed_vars[0] = s_var;

    int size = infile.vector_variable_size("transport/dirichlet/prescribed");
    for (int i = 0; i < size; i++) {
        int prescribed_id = infile("transport/dirichlet/prescribed", -1, i);
        if (prescribed_id != -1) {
            cprescribed.insert(prescribed_id);

        }
    }
    if (cprescribed.size() != 0) {
        Real prescribed_value = infile("dirichlet/sediment/prescribed/value", 1.0);
        ConstFunction<Number> pbc(prescribed_value);
        transport_system.get_dof_map().add_dirichlet_boundary(DirichletBoundary(cprescribed, sed_vars, &pbc));
    }


    this->noflux_bc_id  = infile("transport/neumann/noflux", -1);
    this->erosion_bc_id = infile("transport/neumann/erosion", -1);
    this->deposition_id = infile("transport/deposition", -1);

    bool init           = infile("sediment/init_condition", true);
    if (init) transport_system.attach_init_function(init_sedimentation);

}

void SedimentationTransport::assemble2D() {

    // It is a good idea to make sure we are assembling
    // the proper system.

    PerfLog* perf_log = es.parameters.get<PerfLog*>("PerfLog");
    perf_log->pause_event("Solver", "Transport");
    perf_log->start_event("Assembly", "Transport");

    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system =
            es.get_system<TransientLinearImplicitSystem> ("transport");

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & flow_system =
            es.get_system<TransientLinearImplicitSystem> ("flow");

    // Numeric ids corresponding to each variable in the system
    const unsigned int s_var = system.variable_number("s");

    // Numeric ids corresponding to each variable in the system
    const unsigned int u_var = flow_system.variable_number("u");
    const unsigned int v_var = flow_system.variable_number("v");
    const unsigned int p_var = flow_system.variable_number("p");

    std::vector<dof_id_type> dof_indices_u;
    std::vector<dof_id_type> dof_indices_v;
    std::vector<dof_id_type> dof_indices_p;

    std::vector<dof_id_type> face_dof_indices;


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
    const std::vector<Point> & normals = fe_face->get_normals();

    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
    const std::vector<std::vector<RealGradient> >& dphi_face = fe_face->get_dphi();

    // The XY locations of the quadrature points used for face integration
    const std::vector<Point>& qface_points = fe_face->get_xyz();

    // A reference to the \p DofMap object for this system.  The \p DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the \p DofMap
    // in future examples.
    const DofMap& dof_map = system.get_dof_map();
    const DofMap& dof_map_flow = flow_system.get_dof_map();

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_face;

    // Here we extract the velocity & parameters that we put in the
    // EquationSystems object.
    //const RealVectorValue velocity = es.parameters.get<RealVectorValue> ("velocity");
    const Real dt = es.parameters.get<Real> ("dt");

    const Real fopc = es.parameters.get<Real> ("fopc");
    const Real theta = es.parameters.get<Real> ("theta");
    const Real k = es.parameters.get<Real> ("Diffusivity");
    const Real Us = es.parameters.get<Real> ("Us");
    const Real ex = es.parameters.get<Real> ("ex");
    const Real ey = es.parameters.get<Real> ("ey");
    const Real ez = es.parameters.get<Real> ("ez");
    const Real Rp = es.parameters.get<Real> ("erosion/Rp");
    const Real dt_stab = es.parameters.get<Real> ("dt_stab");


    RealVectorValue e(ex, ey);

    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  Since the mesh
    // will be refined we want to only consider the ACTIVE elements,
    // hence we use a variant of the \p active_elem_iterator.
    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for (; el != end_el; ++el) {
        // Store a pointer to the element we are currently
        // working on.  This allows for nicer syntax later.
        const Elem* elem = *el;

        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        dof_map.dof_indices(elem, dof_indices);

        int n_dofs = dof_indices.size();

        dof_map_flow.dof_indices(elem, dof_indices_u, u_var);
        dof_map_flow.dof_indices(elem, dof_indices_v, v_var);
        dof_map_flow.dof_indices(elem, dof_indices_p, p_var);

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe->reinit(elem);

        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).
        Ke.resize(dof_indices.size(),
                dof_indices.size());

        Fe.resize(dof_indices.size());


        // Compute SUPG stabilization parameters: 
        
        // loop over quadrature points
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
            // Values to hold the old solution & its gradient.
            Number s_old = 0.0, s = 0.0;
            Gradient grad_s_old, grad_s, grad_u, grad_v, grad_p;
            Number u = 0.0, v = 0.0;
            Number u_old = 0.0, v_old = 0.0;

            // Compute the old solution & its gradient.
            for (unsigned int l = 0; l < phi.size(); l++) {

                s_old += phi[l][qp] * system.old_solution(dof_indices[l]);
                grad_s_old.add_scaled(dphi[l][qp], system.old_solution(dof_indices[l]));

                s += phi[l][qp] * system.current_solution(dof_indices[l]);
                grad_s.add_scaled(dphi[l][qp], system.current_solution(dof_indices[l]));

                u += phi[l][qp]*(flow_system.current_solution(dof_indices_u[l]));
                v += phi[l][qp]*(flow_system.current_solution(dof_indices_v[l]));

                grad_u.add_scaled(dphi[l][qp], flow_system.current_solution(dof_indices_u[l]));
                grad_v.add_scaled(dphi[l][qp], flow_system.current_solution(dof_indices_v[l]));

                u_old += phi[l][qp]*(flow_system.old_solution(dof_indices_u[l]));
                v_old += phi[l][qp]*(flow_system.old_solution(dof_indices_v[l]));

                grad_p.add_scaled(dphi[l][qp], flow_system.current_solution(dof_indices_p[l]));

            }


            RealGradient g = compute_g(fe.get(), dim, qp);
            RealTensor G = compute_G(fe.get(), dim, qp);


            Point f = e*s;
            
            RealVectorValue U(u, v);

            Real tau_vms = compute_tau_M(g, G, U, k, dt, dt_stab);

            // Now compute the element matrix and RHS contributions.
            Real rint_x = (u - u_old) / dt + (U * grad_u) + grad_p(0) - f(0);
            Real rint_y = (v - v_old) / dt + (U * grad_v) + grad_p(1) - f(1);

            RealVectorValue velocity(u - tau_vms*rint_x +e(0)*Us, v - tau_vms * rint_y + e(1)*Us);
            Real tau = compute_tau_M(g, G, velocity, k, dt, dt_stab);

            /*
            
            Real unorm = velocity.norm();
            Real aux2 = std::pow(2.0 * unorm / h_caract, 2.0) + aux1;
            Real tau = std::pow(aux2, -0.5);

            // CAU stabilization parameter
            Real res_mass = (s - s_old) / dt;
            Real res_adv = (velocity * grad_s);
            Real residuo = res_mass + res_adv;
            Real gcnorm = grad_s.norm();
            gcnorm = std::max(1.0E-10, gcnorm);
            Real ogcnorm = 1.0 / gcnorm;
            Real aux3 = res_adv / (ogcnorm * ogcnorm);
            RealVectorValue b(grad_s(0) * aux3, grad_s(1));
            Real bnorm = b.norm();
            bnorm = std::max(bnorm, 1.0E-10);
            Real bdb = k * bnorm*bnorm;
            bdb = std::max(bdb, 1.0E-10);
            Real Pe_p = h_caract * (bnorm * bnorm * bnorm) / bdb;
            Real alpha_c = std::min(0.25 * Pe_p, 0.70);
            Real delta_sco = 0.5 * h_caract * alpha_c * residuo * ogcnorm * fopc;
             * */
             
            // Now compute the element matrix and RHS contributions.
            for (unsigned int i = 0; i < phi.size(); i++) {
                // The RHS contribution
                // Galerkin term
                Fe(i) += JxW[qp]*(s_old * phi[i][qp] //mass term
                        -(1.0 - theta) * dt * (
                        // Convection term
                        (grad_s_old * velocity) * phi[i][qp] +
                        // Diffusion term
                        k * (grad_s_old * dphi[i][qp]))
                        );
                // SUPG term        

                Fe(i) += JxW[qp] * tau * (s_old * (velocity * dphi[i][qp])
                        -(1.0 - theta) * dt * (grad_s_old * velocity)*(velocity * dphi[i][qp])
                        );
                for (unsigned int j = 0; j < phi.size(); j++) {
                    // The Galerkin contribution
                    Ke(i, j) += JxW[qp]*(
                            phi[i][qp] * phi[j][qp] + // Mass-matrix
                            theta * dt * ((velocity * dphi[j][qp]) * phi[i][qp] + // Convection
                            k * (dphi[i][qp] * dphi[j][qp]) // Diffusion
                            )
                            );
                    // The SUPG contribution
                    Ke(i, j) += JxW[qp] * tau * (
                            phi[j][qp]*(velocity * dphi[i][qp]) +
                            theta * dt * (velocity * dphi[j][qp])*(velocity * dphi[i][qp])
                            );
                    // YZBetha
                    //Ke(i, j) += JxW[qp] * delta_sco * dt * (dphi[i][qp] * dphi[j][qp]);
                }
            }
        }

        // At this point the interior element integration has
        // been completed.  However, we have not yet addressed
        // boundary conditions.  For this example we will only
        // consider simple Dirichlet boundary conditions imposed
        // via the penalty method.
        //
        // The following loops over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.
        {

            for (unsigned int s = 0; s < elem->n_sides(); s++)
                if (elem->neighbor(s) == NULL) {

                    fe_face->reinit(elem, s);

                    // Applying flux advective boundary condition
                    if (mesh.boundary_info->boundary_id(elem, s) == this->deposition_id) {

                        for (unsigned int qp = 0; qp < qface.n_points(); qp++) {
                            Number s = 0.0, s_old = 0.0;
                            Gradient grad_u, grad_v;

                            //Gradient grad_s;

                            for (unsigned int i = 0; i < phi_face.size(); i++) {
                                s_old += phi_face[i][qp] * system.old_solution(dof_indices[i]);
                                grad_u.add_scaled(dphi_face[i][qp], flow_system.current_solution(dof_indices_u[i]));
                                grad_v.add_scaled(dphi_face[i][qp], flow_system.current_solution(dof_indices_v[i]));
                            }

                            /*
                            RealTensor G;
                            G(0, 0) = k * (2.0 * grad_u(0));
                            G(0, 1) = G(1, 0) = k * (grad_u(1) + grad_v(0));
                            G(1, 1) = k * 2.0 * grad_v(1);
                            RealVectorValue taub;
                            taub(0) = normals[qp](0) * G(0, 0) + normals[qp](1) * G(1, 0);
                            taub(1) = normals[qp](0) * G(0, 1) + normals[qp](1) * G(1, 1);

                            Number taub_norm = taub.norm();
                            const Real A = 1.3E-07;
                            const Real Z = std::sqrt(taub_norm) * std::pow(Rp, 0.6) / Us;
                            const Real Z5 = std::pow(Z, 5.0);

                            Number E = (A * Z5) / (1.0 + A * Z5 / 0.3);
                            */

                            Number tmp = (k/Us);
                            // RHS contribution
                            for (unsigned int i = 0; i < phi_face.size(); i++) {
                                Fe(i) += JxW_face[qp]*tmp*s_old*phi_face[i][qp];
                                for (unsigned int j = 0; j < phi_face.size(); j++)
                                    Ke(i, j) += JxW_face[qp]*tmp*phi_face[i][qp] * phi_face[j][qp];

                            }


                        }
                    }

                    if (mesh.boundary_info->boundary_id(elem, s) == this->noflux_bc_id) {
                        for (unsigned int qp = 0; qp < qface.n_points(); qp++) {

                            Number s_old = 0.0;

                            for (unsigned int i = 0; i < phi_face.size(); i++) {
                                s += phi_face[i][qp] * system.old_solution(dof_indices[i]);

                            }

                            Number tmp = Us * dt * JxW_face[qp];

                            // RHS contribution
                            for (unsigned int i = 0; i < phi_face.size(); i++)
                                Fe(i) += (1 - theta) * tmp * s_old * phi_face[i][qp];

                            //Matrix contribution
                            for (unsigned int i = 0; i < phi_face.size(); i++)
                                for (unsigned int j = 0; j < phi_face.size(); j++)
                                    Ke(i, j) += theta * tmp * phi_face[i][qp] * phi_face[j][qp];

                        }

                    }
                }

        }




        // If this assembly program were to be used on an adaptive mesh,
        // we would have to apply any hanging node constraint equations
        dof_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The \p SparseMatrix::add_matrix()
        // and \p NumericVector::add_vector() members do this for us.
        system.matrix->add_matrix(Ke, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);
    }
    // That concludes the system matrix assembly routine.

    perf_log->stop_event("Assembly", "Transport");
    perf_log->restart_event("Solver", "Transport");

}

void SedimentationTransport::assemble() {
    switch (this->dim) {
        case 2:
            this->assemble2D();
            break;
        default:
            this->assemble3D();
            break;
    }
}

void SedimentationTransport::assemble3D() {

    // It is a good idea to make sure we are assembling
    // the proper system.

    PerfLog* perf_log = es.parameters.get<PerfLog*>("PerfLog");
    perf_log->pause_event("Solver", "Transport");
    perf_log->start_event("Assembly", "Transport");


    
    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system =
            es.get_system<TransientLinearImplicitSystem> ("transport");

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & flow_system =
            es.get_system<TransientLinearImplicitSystem> ("flow");

    // Numeric ids corresponding to each variable in the system
    const unsigned int s_var = system.variable_number("s");

    // Numeric ids corresponding to each variable in the system
    const unsigned int u_var = flow_system.variable_number("u");
    const unsigned int v_var = flow_system.variable_number("v");
    const unsigned int w_var = flow_system.variable_number("w");
    const unsigned int p_var = flow_system.variable_number("p");

    std::vector<dof_id_type> dof_indices_u;
    std::vector<dof_id_type> dof_indices_v;
    std::vector<dof_id_type> dof_indices_w;
    std::vector<dof_id_type> dof_indices_p;

    std::vector<dof_id_type> face_dof_indices;


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
    const std::vector<Real>& JxW      = fe->get_JxW();
    const std::vector<Real>& JxW_face = fe_face->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real> >& phi      = fe->get_phi();
    const std::vector<std::vector<Real> >& phi_face = fe_face->get_phi();

    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<RealGradient> >& dphi      = fe->get_dphi();
    const std::vector<std::vector<RealGradient> >& dphi_face = fe_face->get_dphi();
    const std::vector<Point> & normals = fe_face->get_normals();

    // The XY locations of the quadrature points used for face integration
    const std::vector<Point>& qface_points = fe_face->get_xyz();

    // A reference to the \p DofMap object for this system.  The \p DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the \p DofMap
    // in future examples.
    const DofMap& dof_map      = system.get_dof_map();
    const DofMap& dof_map_flow = flow_system.get_dof_map();

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_face;

    // Here we extract the velocity & parameters that we put in the
    // EquationSystems object.
    //const RealVectorValue velocity = es.parameters.get<RealVectorValue> ("velocity");
    const Real dt = es.parameters.get<Real> ("dt");

    const Real fopc = es.parameters.get<Real> ("fopc");
    const Real theta = es.parameters.get<Real> ("theta");
    const Real k = es.parameters.get<Real> ("Diffusivity");
    const Real Us = es.parameters.get<Real> ("Us");
    const Real ex = es.parameters.get<Real> ("ex");
    const Real ey = es.parameters.get<Real> ("ey");
    const Real ez = es.parameters.get<Real> ("ez");
    const Real Rp = es.parameters.get<Real> ("erosion/Rp");
    const Real dt_stab = es.parameters.get<Real> ("dt_stab");
    
    RealVectorValue e(ex, ey, ez);
                  
    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  Since the mesh
    // will be refined we want to only consider the ACTIVE elements,
    // hence we use a variant of the \p active_elem_iterator.
    MeshBase::const_element_iterator el          = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for (; el != end_el; ++el) {
        // Store a pointer to the element we are currently
        // working on.  This allows for nicer syntax later.
        const Elem* elem = *el;

        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        dof_map.dof_indices(elem, dof_indices);

        int n_dofs = dof_indices.size();

        dof_map_flow.dof_indices(elem, dof_indices_u, u_var);
        dof_map_flow.dof_indices(elem, dof_indices_v, v_var);
        dof_map_flow.dof_indices(elem, dof_indices_w, w_var);
        dof_map_flow.dof_indices(elem, dof_indices_p, p_var);

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe->reinit(elem);

        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).
        Ke.resize(dof_indices.size(),
                dof_indices.size());

        Fe.resize(dof_indices.size());


        // Compute SUPG stabilization parameters:
        // The carateristic height of the element

        // loop over quadrature points
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
            // Values to hold the old solution & its gradient.
            Number s_old = 0.0, s = 0.0;
            Gradient grad_s_old, grad_s;
            Gradient grad_u, grad_v, grad_w, grad_p;
            Number u = 0.0, v = 0.0, w = 0.0;
            Number u_old = 0.0, v_old = 0.0, w_old = 0.0;

            // Compute the old solution & its gradient.
            for (unsigned int l = 0; l < phi.size(); l++) {

                s_old += phi[l][qp] * system.old_solution(dof_indices[l]);
                grad_s_old.add_scaled(dphi[l][qp], system.old_solution(dof_indices[l]));

                s += phi[l][qp] * system.current_solution(dof_indices[l]);
                grad_s.add_scaled(dphi[l][qp], system.current_solution(dof_indices[l]));

                u += phi[l][qp]*(flow_system.current_solution(dof_indices_u[l]));
                v += phi[l][qp]*(flow_system.current_solution(dof_indices_v[l]));
                w += phi[l][qp]*(flow_system.current_solution(dof_indices_w[l]));

                u_old += phi[l][qp]*(flow_system.old_solution(dof_indices_u[l]));
                v_old += phi[l][qp]*(flow_system.old_solution(dof_indices_v[l]));
                w_old += phi[l][qp]*(flow_system.old_solution(dof_indices_w[l]));

                grad_u.add_scaled(dphi[l][qp], flow_system.current_solution(dof_indices_u[l]));
                grad_v.add_scaled(dphi[l][qp], flow_system.current_solution(dof_indices_v[l]));
                grad_w.add_scaled(dphi[l][qp], flow_system.current_solution(dof_indices_w[l]));
                grad_p.add_scaled(dphi[l][qp], flow_system.current_solution(dof_indices_p[l]));


            }


            
            RealGradient g = compute_g(fe.get(), dim, qp);
            RealTensor G = compute_G(fe.get(), dim, qp);

            RealVectorValue U(u, v, w);
            Point f = e*s_old;

            Real tau_vms = compute_tau_M(g, G, U, k, dt, dt_stab);

            // Now compute the element matrix and RHS contributions.

            

            Real rint_x = (u - u_old) / dt + (U * grad_u) + grad_p(0) - f(0);
            Real rint_y = (v - v_old) / dt + (U * grad_v) + grad_p(1) - f(1);
            Real rint_z = (w - w_old) / dt + (U * grad_w) + grad_p(1) - f(2);

            RealVectorValue velocity(u - tau_vms*rint_x + e(0)*Us, v - tau_vms * rint_y + e(1)*Us, w - tau_vms * rint_z + e(2)*Us);

            Real tau = compute_tau_M(g, G, velocity, k, dt, dt_stab);

/*
            const Real res = (s - s_old) / dt + velocity*grad_s;
            const Real unorm = velocity.norm();
            const Real aux2 = pow(2.0 * unorm / h_caract, 2.0) + aux1;
            const Real tau = pow(aux2, -0.5);


            // CAU stabilization parameter
            Real res_mass = (s - s_old) / dt;
            Real res_adv = (velocity * grad_s);
            Real residuo = res_mass + res_adv;
            Real gcnorm = grad_s.norm();
            gcnorm = std::max(1.0E-10, gcnorm);
            Real ogcnorm = 1.0 / gcnorm;
            Real aux3 = res_adv / (ogcnorm * ogcnorm);
            RealVectorValue b(grad_s(0) * aux3, grad_s(1) * aux3, grad_s(2) * aux3);
            Real bnorm = b.norm();
            bnorm = std::max(bnorm, 1.0E-10);
            Real bdb = k * bnorm*bnorm;
            bdb = std::max(bdb, 1.0E-10);
            Real Pe_p = 0.5 * h_caract * (bnorm * bnorm * bnorm) / bdb;
            Real alpha_c = std::min(0.25 * Pe_p, 0.70);
            Real delta_sco = 0.5 * h_caract * alpha_c * residuo * ogcnorm * fopc;
*/

            // Now compute the element matrix and RHS contributions.
            for (unsigned int i = 0; i < phi.size(); i++) {
                // The RHS contribution
                // Galerkin term
                Fe(i) += JxW[qp]*(s_old * phi[i][qp] //mass term
                        -(1.0 - theta) * dt * (
                        // Convection term
                        (grad_s_old * velocity) * phi[i][qp] +
                        // Diffusion term
                        k * (grad_s_old * dphi[i][qp]))
                        );
                // SUPG term
                Fe(i) += JxW[qp] * tau * (s_old * (dphi[i][qp] * velocity)
                        -(1 - theta) * dt * (grad_s_old * velocity)*(velocity * dphi[i][qp])
                        );
                for (unsigned int j = 0; j < phi.size(); j++) {
                    // The Galerkin contribution
                    Ke(i, j) += JxW[qp]*(
                            phi[i][qp] * phi[j][qp] + // Mass-matrix
                            theta * dt * ((velocity * dphi[j][qp]) * phi[i][qp] + // Convection
                            k * (dphi[i][qp] * dphi[j][qp]) // Diffusion
                            )
                            );

                    // The SUPG contribution
                    Ke(i, j) += JxW[qp] * tau * (
                            phi[j][qp]*(velocity * dphi[i][qp]) +
                            theta * dt * (velocity * dphi[j][qp])*(velocity * dphi[i][qp])
                            );
                    // CAU
                    //Ke(i, j) += JxW[qp] * delta_sco * dt * (dphi[i][qp] * dphi[j][qp]);
                }
            }
        }

        // At this point the interior element integration has
        // been completed.  However, we have not yet addressed
        // boundary conditions.  For this example we will only
        // consider simple Dirichlet boundary conditions imposed
        // via the penalty method.
        {


            for (unsigned int s = 0; s < elem->n_sides(); s++)
                if (elem->neighbor(s) == NULL) {

                    fe_face->reinit(elem, s);

                    // Applying flux advective boundary condition
                    if (mesh.boundary_info->boundary_id(elem, s) == this->erosion_bc_id) {

                        for (unsigned int qp = 0; qp < qface.n_points(); qp++) {
                            Number s_old = 0.0;
                            Gradient grad_u, grad_v, grad_w;

                            Gradient grad_s;

                            for (unsigned int i = 0; i < phi_face.size(); i++) {
                                s_old += phi_face[i][qp] * system.old_solution(dof_indices[i]);
                                grad_u.add_scaled(dphi_face[i][qp], flow_system.current_solution(dof_indices_u[i]));
                                grad_v.add_scaled(dphi_face[i][qp], flow_system.current_solution(dof_indices_v[i]));
                                grad_w.add_scaled(dphi_face[i][qp], flow_system.current_solution(dof_indices_w[i]));
                            }

                            /*
                            DenseMatrix<Number> G;

                            G(0, 0) = k * 2.0 * grad_u(0);
                            G(0, 1) = G(1, 0) = k * (grad_u(1) + grad_v(0));
                            G(0, 2) = G(2, 0) = k * (grad_u(2) + grad_w(0));
                            G(1, 1) = k * 2.0 * grad_v(1);
                            G(1, 2) = G(2, 1) = k * (grad_v(2) + grad_w(1));
                            G(2, 2) = k * 2.0 * grad_w(2);

                            RealVectorValue taub;
                            taub(0) = normals[qp](0) * G(0, 0) + normals[qp](1) * G(1, 0) + normals[qp](2) * G(2, 0);
                            taub(1) = normals[qp](0) * G(0, 1) + normals[qp](1) * G(1, 1) + normals[qp](2) * G(2, 1);
                            taub(2) = normals[qp](0) * G(0, 2) + normals[qp](1) * G(1, 2) + normals[qp](2) * G(2, 2);

                            Number taub_norm = taub.norm();
                            const Real A = 1.3E-07;
                            const Real Z = std::sqrt(taub_norm) * std::pow(Rp, 0.6) / Us;
                            const Real Z5 = std::pow(Z, 5.0);

                            Number E = (A * Z5) / (1.0 + A * Z5 / 0.3);
                            */
                            
                            Number tmp = (k/Us);
                            // RHS contribution
                            for (unsigned int i = 0; i < phi_face.size(); i++) {
                                Fe(i) += JxW_face[qp]*tmp*s_old*phi_face[i][qp];
                                for (unsigned int j = 0; j < phi_face.size(); j++)
                                    Ke(i, j) += JxW_face[qp]*tmp*phi_face[i][qp] * phi_face[j][qp];

                            }

                        }
                    }

                    if (mesh.boundary_info->boundary_id(elem, s) == this->noflux_bc_id) {
                        for (unsigned int qp = 0; qp < qface.n_points(); qp++) {
                            Number s = 0.0;

                            for (unsigned int i = 0; i < phi_face.size(); i++) {
                                s += phi_face[i][qp] * system.old_solution(dof_indices[i]);

                            }

                            Number tmp = (1.0 - theta) * Us * dt * JxW_face[qp];
                            // RHS contribution
                            for (unsigned int i = 0; i < phi_face.size(); i++)
                                Fe(i) += tmp * s * phi_face[i][qp];

                            //Matrix contribution
                            for (unsigned int i = 0; i < phi_face.size(); i++)
                                for (unsigned int j = 0; j < phi_face.size(); j++)
                                    Ke(i, j) += tmp * phi_face[i][qp] * phi_face[j][qp];

                        }

                    }
                }

        }

        
        // If this assembly program were to be used on an adaptive mesh,
        // we would have to apply any hanging node constraint equations
        dof_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The \p SparseMatrix::add_matrix()
        // and \p NumericVector::add_vector() members do this for us.
        system.matrix->add_matrix(Ke, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);
    }
    // That concludes the system matrix assembly routine.
    

    perf_log->stop_event("Assembly", "Transport");
    perf_log->restart_event("Solver", "Transport");

}
