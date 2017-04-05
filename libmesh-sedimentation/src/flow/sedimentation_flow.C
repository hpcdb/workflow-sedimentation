#include "sedimentation_flow.h"

#include "libmesh/parsed_function.h"
#include "libmesh/zero_function.h"
#include "libmesh/perf_log.h"

#include "define.h"

#include "stab_helper.h"

#include "mesh_moviment.h"



// Boundary conditions for the 2D test case
class HydrostaticPressureFunction: public FunctionBase<Number>
{
public:
    
    HydrostaticPressureFunction (int dim, Real rho, Real gravity):_dim(dim), _rho(rho), _gravity(gravity) 
    { 
        this->_initialized = true; 
    }
    

    virtual Number operator() (const Point &, const Real = 0)
    { libmesh_not_implemented(); }

    virtual void operator() (const Point & p,
                             const Real,
                             DenseVector<Number> & output)
    {
        output(0) = _rho * _gravity * p(_dim-1);
    }

  virtual UniquePtr<FunctionBase<Number> > clone() const
  { 
      return UniquePtr<FunctionBase<Number> > (new HydrostaticPressureFunction(_dim, _rho, _gravity)); 
  }

private:
    int _dim;
    Real _rho;
    Real _gravity;
};


void SedimentationFlow::init() {

    MeshBase& mesh = this->es.get_mesh();
    this->dim = mesh.mesh_dimension();
    TransientLinearImplicitSystem & flow_system = this->es.add_system<TransientLinearImplicitSystem> ("flow");

    unsigned int u_var = flow_system.add_variable("u", FIRST);
    unsigned int v_var = flow_system.add_variable("v", FIRST);
    unsigned int w_var;
    if (dim == 3)
        w_var = flow_system.add_variable("w", FIRST);

    unsigned int p_var = flow_system.add_variable("p", FIRST);

}

void SedimentationFlow::setup(GetPot &infile) 
{
    
    MeshBase& mesh = this->es.get_mesh();

    this->dim = mesh.mesh_dimension();

    TransientLinearImplicitSystem & flow_system = this->es.get_system<TransientLinearImplicitSystem> ("flow");

    flow_system.attach_assemble_object(*this);

    unsigned int u_var = flow_system.variable_number("u");
    unsigned int v_var = flow_system.variable_number("v");
    unsigned int w_var = 0;
    if (dim == 3)
        w_var = flow_system.variable_number("w");

    unsigned int p_var = flow_system.variable_number("p");


    std::set<boundary_id_type> slipx;   // 1
    std::set<boundary_id_type> slipy;   // 2
    std::set<boundary_id_type> slipz;   // 3
    std::set<boundary_id_type> pnull;   // 4
    std::set<boundary_id_type> noslip;  // 5
    std::set<boundary_id_type> inlet;   // 6
    

    std::vector<unsigned int> velocity_var(2);
    velocity_var[0] = u_var;
    velocity_var[1] = v_var;
    if (this->dim == 3)
        velocity_var.push_back(w_var);

    std::vector<unsigned int> velocity_x(1);
    velocity_x[0] = u_var;

    std::vector<unsigned int> velocity_y(1);
    velocity_y[0] = v_var;

    std::vector<unsigned int> velocity_z(1);
    velocity_z[0] = w_var;

    std::vector<unsigned int> pressure(1);
    pressure[0] = p_var;

    ZeroFunction<Number> zero;
    
    
    //SET Boundary Condition
    // Velocity Field

    int size = infile.vector_variable_size("flow/dirichlet/slipx");
    for (int i = 0; i < size; i++) 
    {
        int slipx_id = infile("flow/dirichlet/slipx", -1, i);
        if (slipx_id != -1) {
            slipx.insert(slipx_id);
            
        }
    }
    if (slipx.size() != 0) 
        flow_system.get_dof_map().add_dirichlet_boundary(DirichletBoundary(slipx, velocity_x, &zero));


    size = infile.vector_variable_size("flow/dirichlet/slipy");
    for (int i = 0; i < size; i++) {
        int slipy_id = infile("flow/dirichlet/flow/slipy", -1, i);
        if (slipy_id != -1) {
            slipy.insert(slipy_id);

        }
    }
    if (slipy.size() != 0)
        flow_system.get_dof_map().add_dirichlet_boundary(DirichletBoundary(slipy, velocity_y, &zero));


    size = infile.vector_variable_size("flow/dirichlet/slipz");
    for (int i = 0; i < size; i++) {
        int slipz_id = infile("flow/dirichlet/slipz", -1, i);
        if (slipz_id != -1) {
            slipz.insert(slipz_id);
        }
    }
    
    if (slipz.size() != 0)
        flow_system.get_dof_map().add_dirichlet_boundary(DirichletBoundary(slipz, velocity_z, &zero));
    
        
    size = infile.vector_variable_size("flow/dirichlet/noslip");
    for (int i = 0; i < size; i++) {
        int noslip_id = infile("flow/dirichlet/noslip", -1, i);
        if (noslip_id != -1) {
            noslip.insert(noslip_id);
    
        }
    }
    if (noslip.size() != 0)
        flow_system.get_dof_map().add_dirichlet_boundary(DirichletBoundary(noslip, velocity_var, &zero));
        
    size = infile.vector_variable_size("flow/dirichlet/inlet");
    for (int i = 0; i < size; i++) {
        int inlet_id = infile("flow/dirichlet/inlet", -1, i);
        if (inlet_id != -1) {
            inlet.insert(inlet_id);
        }
    }
    if (inlet.size() != 0) {

        if (infile.vector_variable_size("flow/dirichlet/inlet/constant/u") == 1) {
            Real u_value = infile("flow/dirichlet/inlet/constant/u", 0.0, 0);
            ConstFunction<Number> ubc(u_value);
            flow_system.get_dof_map().add_dirichlet_boundary(DirichletBoundary(inlet, velocity_x, &ubc));
        }
        if (infile.vector_variable_size("flow/dirichlet/inlet/constant/v") == 1) {
            Real v_value = infile("flow/dirichlet/inlet/constant/v", 0.0, 0);
            ConstFunction<Number> vbc(v_value);
            flow_system.get_dof_map().add_dirichlet_boundary(DirichletBoundary(inlet, velocity_y, &vbc));
        }
        if (infile.vector_variable_size("flow/dirichlet/inlet/constant/w") == 1) {
            Real w_value = infile("flow/dirichlet/inlet/constant/w", 0.0, 0);
            ConstFunction<Number> wbc(w_value);
            flow_system.get_dof_map().add_dirichlet_boundary(DirichletBoundary(inlet, velocity_z, &wbc));
        }
        
        if (infile.vector_variable_size("flow/dirichlet/inlet/functiont/u") == 1) {
            std::string f = infile("flow/dirichlet/inlet/function/u", "0.0", 0);
            flow_system.get_dof_map().add_dirichlet_boundary(DirichletBoundary(inlet, velocity_x, ParsedFunction<Number>(f)));
        }
        if (infile.vector_variable_size("flow/dirichlet/inlet/function/v") == 1) {
            std::string f = infile("flow/dirichlet/inlet/function/v", "0.0", 0);
            flow_system.get_dof_map().add_dirichlet_boundary(DirichletBoundary(inlet, velocity_y, ParsedFunction<Number>(f)));
        }
        if (infile.vector_variable_size("flow/dirichlet/inlet/function/w") == 1) {
            std::string f =infile("flow/dirichlet/inlet/function/w", "0.0", 0);
            flow_system.get_dof_map().add_dirichlet_boundary(DirichletBoundary(inlet, velocity_z, ParsedFunction<Number>(f)));
        }
    }
    
    
    size = infile.vector_variable_size("flow/dirichlet/pressure/zero");
    for (int i = 0; i < size; i++) {
        int pnull_id = infile("flow/dirichlet/pressure/zero", -1, i);
        if (pnull_id != -1) {
            pnull.insert(pnull_id);
        }
    }
    if (pnull.size() != 0)
        flow_system.get_dof_map().add_dirichlet_boundary(DirichletBoundary(pnull, pressure, &zero));
 
    pnull.clear();
    
    size = infile.vector_variable_size("flow/dirichlet/pressure/hydrostatic");
    for (int i = 0; i < size; i++) {
        int pnull_id = infile("flow/dirichlet/pressure/hydrostatic", -1, i);
        if (pnull_id != -1) {
            pnull.insert(pnull_id);
        }
    }
    if (pnull.size() != 0) {
        flow_system.get_dof_map().add_dirichlet_boundary(DirichletBoundary(pnull, pressure,HydrostaticPressureFunction(dim, rho, gravity)));
    }

    this->outbflow_id = infile("flow/nenmann/outflow", -1);
    
}

void SedimentationFlow::assemble() {
    switch (this->dim) {
        case 2:
            this->assemble2D();
            break;
        default:
            this->assemble3D();
            break;
    }
}


// The matrix assembly function to be called at each time step to
// prepare for the linear solve.

void SedimentationFlow::assemble3D() {
    
    // It is a good idea to make sure we are assembling
    // the proper system.

    PerfLog* perf_log = es.parameters.get<PerfLog*>("PerfLog");
    perf_log->pause_event("Solver", "Flow");
    perf_log->start_event("Assembly", "Flow");

    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the Stokes system object.
    TransientLinearImplicitSystem & navier_stokes_system =
            es.get_system<TransientLinearImplicitSystem> ("flow");

    TransientLinearImplicitSystem & transport_system =
            es.get_system<TransientLinearImplicitSystem> ("transport");

    // Numeric ids corresponding to each variable in the system
    const unsigned int u_var = navier_stokes_system.variable_number("u");
    const unsigned int v_var = navier_stokes_system.variable_number("v");
    const unsigned int w_var = navier_stokes_system.variable_number("w");
    const unsigned int p_var = navier_stokes_system.variable_number("p");
    const unsigned int s_var = transport_system.variable_number("s");


#ifdef MESH_MOVIMENT
    LinearImplicitSystem & mesh_system =
            es.get_system<LinearImplicitSystem> ("MeshMoving");
    const unsigned int disp_var = mesh_system.variable_number("disp-z");
#endif


    FEType fe_vel_type = navier_stokes_system.variable_type(u_var);

    UniquePtr<FEBase> fe_vel(FEBase::build(dim, fe_vel_type));
    QGauss qrule(dim, fe_vel_type.default_quadrature_order());


    // Tell the finite element objects to use our quadrature rule.
    fe_vel->attach_quadrature_rule(&qrule);
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real>& JxW = fe_vel->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();


    // A reference to the \p DofMap object for this system.  The \p DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the \p DofMap
    // in future examples.
    const DofMap & dof_map = navier_stokes_system.get_dof_map();
    const DofMap & dof_map_sed = transport_system.get_dof_map();
    
#ifdef MESH_MOVIMENT      
    const DofMap & dof_map_mesh = mesh_system.get_dof_map();
    std::vector<dof_id_type> dof_indices_mesh;
#endif
    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;


    DenseSubMatrix<Number>
            Kuu(Ke), Kuv(Ke), Kuw(Ke), Kup(Ke),
            Kvu(Ke), Kvv(Ke), Kvw(Ke), Kvp(Ke),
            Kwu(Ke), Kwv(Ke), Kww(Ke), Kwp(Ke),
            Kpu(Ke), Kpv(Ke), Kpw(Ke), Kpp(Ke);

    DenseSubVector<Number>
            Fu(Fe),
            Fv(Fe),
            Fw(Fe),
            Fp(Fe);

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_u;
    std::vector<dof_id_type> dof_indices_v;
    std::vector<dof_id_type> dof_indices_w;
    std::vector<dof_id_type> dof_indices_p;
    std::vector<dof_id_type> dof_indices_s;


    const Real dt = es.parameters.get<Real>("dt");
    const Real dt_stab = es.parameters.get<Real>("dt_stab");
    const Real one_dt = 1.0 / dt;
    Real Reynolds = es.parameters.get<Real>("Reynolds");
    Real ex = es.parameters.get<Real>("ex");
    Real ey = es.parameters.get<Real>("ey");
    Real ez = es.parameters.get<Real>("ez");
    const Real uRe = 1.0 / Reynolds;

    const NumberVectorValue norm(ex, ey, ez);
    //
    const Real one3 = 1.0 / 3.0;
    const Real invPi = 1.0 / libMesh::pi;

    // Now we will loop over all the elements in the mesh that
    // live on the local processor.
    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for (; el != end_el; ++el) {
        // Store a pointer to the element we are currently
        // working on.  This allows for nicer syntax later.
        const Elem* elem = *el;

        // Compute SUPG stabilization parameters:
        // The carateristic height of the element
        //const Real volume = elem->volume();
        //const Real h_caract = pow(6.0 * volume*invPi, one3);
        //const Real aux1 = 9.0 * pow(4.0 * uRe / (h_caract * h_caract), 2.0)+ (4.0 / (dt * dt)) * dt_stab;


        dof_map.dof_indices(elem, dof_indices);
        dof_map.dof_indices(elem, dof_indices_u, u_var);
        dof_map.dof_indices(elem, dof_indices_v, v_var);
        dof_map.dof_indices(elem, dof_indices_w, w_var);
        dof_map.dof_indices(elem, dof_indices_p, p_var);
        dof_map_sed.dof_indices(elem, dof_indices_s, s_var);

#ifdef MESH_MOVIMENT      
        dof_map_mesh.dof_indices(elem, dof_indices_mesh, disp_var);
#endif

        const unsigned int n_dofs = dof_indices.size();
        const unsigned int n_u_dofs = dof_indices_u.size();
        const unsigned int n_v_dofs = dof_indices_v.size();
        const unsigned int n_w_dofs = dof_indices_w.size();
        const unsigned int n_p_dofs = dof_indices_p.size();

        // Compute the element-specific data for the current
        // element.
        fe_vel->reinit(elem);


        // Zero the element matrix and right-hand side before
        // summing them.
        Ke.resize(n_dofs, n_dofs);
        Fe.resize(n_dofs);

        //
        // Similarly, the \p DenseSubVector.reposition () member
        // takes the (row_offset, row_size)
        Kuu.reposition(u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuv.reposition(u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
        Kuw.reposition(u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);
        Kup.reposition(u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);

        Kvu.reposition(v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
        Kvv.reposition(v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
        Kvw.reposition(v_var*n_v_dofs, w_var*n_v_dofs, n_v_dofs, n_w_dofs);
        Kvp.reposition(v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);

        Kwu.reposition(w_var*n_w_dofs, u_var*n_w_dofs, n_w_dofs, n_u_dofs);
        Kwv.reposition(w_var*n_w_dofs, v_var*n_w_dofs, n_w_dofs, n_v_dofs);
        Kww.reposition(w_var*n_w_dofs, w_var*n_w_dofs, n_w_dofs, n_w_dofs);
        Kwp.reposition(w_var*n_w_dofs, p_var*n_w_dofs, n_w_dofs, n_p_dofs);

        Kpu.reposition(p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
        Kpv.reposition(p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
        Kpw.reposition(p_var*n_u_dofs, w_var*n_u_dofs, n_p_dofs, n_w_dofs);
        Kpp.reposition(p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);

        Fu.reposition(u_var*n_u_dofs, n_u_dofs);
        Fv.reposition(v_var*n_u_dofs, n_v_dofs);
        Fw.reposition(w_var*n_u_dofs, n_w_dofs);
        Fp.reposition(p_var*n_u_dofs, n_p_dofs);

        // Now we will build the element matrix and right-hand-side.
        // Constructing the RHS requires the solution and its
        // gradient from the previous timestep.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
            // Values to hold the solution & its gradient at the previous timestep.
            Number u = 0., u_old = 0.;
            Number v = 0., v_old = 0.;
            Number w = 0., w_old = 0.;
            Number p_old;
            Number c = 0., fg = 0.;
            Number u_mesh = 0.0, v_mesh = 0.0, w_mesh = 0.0;
            Gradient grad_u, grad_v, grad_w, grad_p;
            Gradient grad_u_old, grad_v_old, grad_w_old, grad_p_old;
            Point f;
            f.zero();


            for (unsigned int l = 0; l < n_u_dofs; l++) {
                // From the old timestep:
                u_old += phi[l][qp] * navier_stokes_system.old_solution(dof_indices_u[l]);
                v_old += phi[l][qp] * navier_stokes_system.old_solution(dof_indices_v[l]);
                w_old += phi[l][qp] * navier_stokes_system.old_solution(dof_indices_w[l]);
                p_old += phi[l][qp] * navier_stokes_system.old_solution(dof_indices_p[l]);

                grad_u_old.add_scaled(dphi[l][qp], navier_stokes_system.old_solution(dof_indices_u[l]));
                grad_v_old.add_scaled(dphi[l][qp], navier_stokes_system.old_solution(dof_indices_v[l]));
                grad_w_old.add_scaled(dphi[l][qp], navier_stokes_system.old_solution(dof_indices_w[l]));
                grad_p_old.add_scaled(dphi[l][qp], navier_stokes_system.old_solution(dof_indices_p[l]));

                
#ifdef MESH_MOVIMENT
                double w = one_dt * mesh_system.current_solution(dof_indices_mesh[l]);
                w_mesh += phi[l][qp] * w;
#endif

                // From the previous Newton iterate:
                u += phi[l][qp] * navier_stokes_system.current_solution(dof_indices_u[l]);
                v += phi[l][qp] * navier_stokes_system.current_solution(dof_indices_v[l]);
                w += phi[l][qp] * navier_stokes_system.current_solution(dof_indices_w[l]);

                c += phi[l][qp] * transport_system.current_solution(dof_indices_s[l]);

                grad_u.add_scaled(dphi[l][qp], navier_stokes_system.current_solution(dof_indices_u[l]));
                grad_v.add_scaled(dphi[l][qp], navier_stokes_system.current_solution(dof_indices_v[l]));
                grad_w.add_scaled(dphi[l][qp], navier_stokes_system.current_solution(dof_indices_w[l]));
                grad_p.add_scaled(dphi[l][qp], navier_stokes_system.current_solution(dof_indices_p[l]));

            }

            f = norm*c;


            RealGradient g = compute_g(fe_vel.get(), 3, qp);
            RealTensor   G = compute_G(fe_vel.get(), 3, qp);


            // Definitions for convenience.  It is sometimes simpler to do a
            // dot product if you have the full vector at your disposal.
            const NumberVectorValue U_old(u_old, v_old, w_old);
            const NumberVectorValue U(u - u_mesh, v - v_mesh, w - w_mesh);

            Real tau_m = compute_tau_M(g, G, U, uRe, dt, dt_stab);
            

            // Residual VMS term
            Real r_m_x = (u - u_old) / dt + (U * grad_u) + grad_p(0) - f(0);
            Real r_m_y = (v - v_old) / dt + (U * grad_v) + grad_p(1) - f(1);
            Real r_m_z = (w - w_old) / dt + (U * grad_w) + grad_p(2) - f(2);
            Real r_c = grad_u(0) + grad_v(1) + grad_w(2);


            // SUPG/PSPG/LSIC parameters
            //Number unorm  = U.size();
            //Number aux2   = pow(2.0*unorm/h_caract,2.0) + aux1;
            //Number tau    = pow(aux2,-0.5);

            
            const NumberVectorValue Uvms(U(0) - tau_m*r_m_x - u_mesh, U(1) - tau_m*r_m_y - v_mesh, U(2) - tau_m * r_m_z - w_mesh);
  
            tau_m            = compute_tau_M(g, G, Uvms, uRe, dt, dt_stab);
            const Real tau_c = compute_tau_C(g, tau_m);
             

            // First, an i-loop over the velocity degrees of freedom.
            // We know that n_u_dofs == n_v_dofs so we can compute contributions
            // for both at the same time.
            for (unsigned int i = 0; i < n_u_dofs; i++) {

                const Number Udphi_i = Uvms * dphi[i][qp];

                
                Fu(i) += JxW[qp]*((phi[i][qp] + tau_m * Udphi_i) * u_old + // mass-matrix term
                        dt * (phi[i][qp] + tau_m * Udphi_i) * f(0)); // SUPG term

                Fv(i) += JxW[qp]*((phi[i][qp] + tau_m * Udphi_i) * v_old + // mass-matrix term
                        dt * (phi[i][qp] + tau_m * Udphi_i) * f(1)); // SUPG term

                Fw(i) += JxW[qp]*((phi[i][qp] + tau_m * Udphi_i) * w_old + // mass-matrix term
                        dt * (phi[i][qp] + tau_m * Udphi_i) * f(2));

                // SUPG term
                Fp(i) += JxW[qp] * tau_m * (U_old * dphi[i][qp] // PSPG term
                            + dt * (dphi[i][qp] * norm) * fg  //PSPG buoyancy vector
                        );

                // Matrix contributions for the uu and vv couplings.
                for (unsigned int j = 0; j < n_u_dofs; j++) {
                    Number DphiUiUjDphj = 0.0;
                    //DphiUiUjDphj = (Udphi_i)*(Uvms * dphi[j][qp]);
                    
                    for (int d =0; d < dim; d++)
                    {
                        DphiUiUjDphj += (dphi[i][qp](0) * Uvms(0) * Uvms(d) + dphi[i][qp](1) * Uvms(1) * Uvms(d) +
                                dphi[i][qp](2) * Uvms(2) * Uvms(d)) * dphi[j][qp](d);
                    }
                    

                    const Number Udphi_j = Uvms * dphi[j][qp];

                    // Galerkin contribution
                    // G.1: u row
                    Kuu(i, j) += JxW[qp]*(phi[i][qp] * phi[j][qp] + // mass matrix term
                            dt * uRe * (2.0 * dphi[i][qp](0) * dphi[j][qp](0) +
                            dphi[i][qp](1) * dphi[j][qp](1) +
                            dphi[i][qp](2) * dphi[j][qp](2)) + // xx diffusion term
                            dt * Udphi_j * phi[i][qp]); // convection term

                    Kuv(i, j) += JxW[qp] * dt * uRe * (dphi[i][qp](1) * dphi[j][qp](0)); // xy diffusion term
                    Kuw(i, j) += JxW[qp] * dt * uRe * (dphi[i][qp](2) * dphi[j][qp](0)); // xz diffusion term
                    Kup(i, j) += JxW[qp]*(-dt * dphi[i][qp](0) * phi[j][qp]);

                    // G.2: v row
                    Kvu(i, j) += JxW[qp] * dt * uRe * (dphi[i][qp](0) * dphi[j][qp](1)); // xy diffusion term
                    Kvv(i, j) += JxW[qp]*(phi[i][qp] * phi[j][qp] + // mass matrix term
                            dt * uRe * (2.0 * dphi[i][qp](1) * dphi[j][qp](1) +
                            dphi[i][qp](0) * dphi[j][qp](0) +
                            dphi[i][qp](2) * dphi[j][qp](2)) + // diffusion term
                            dt * (Udphi_j) * phi[i][qp]); // convection term
                    Kvw(i, j) += JxW[qp] * dt * uRe * (dphi[i][qp](2) * dphi[j][qp](1)); // xy diffusion term
                    Kvp(i, j) += JxW[qp]*(-dt * dphi[i][qp](1) * phi[j][qp]);

                    // G.3: w row
                    Kwu(i, j) += JxW[qp] * dt * uRe * (dphi[i][qp](0) * dphi[j][qp](2)); // xy diffusion term
                    Kwv(i, j) += JxW[qp] * dt * uRe * (dphi[i][qp](1) * dphi[j][qp](2)); // xy diffusion term
                    Kww(i, j) += JxW[qp]*(phi[i][qp] * phi[j][qp] + // mass matrix term
                            dt * uRe * (2.0 * dphi[i][qp](2) * dphi[j][qp](2) +
                            dphi[i][qp](1) * dphi[j][qp](1) +
                            dphi[i][qp](0) * dphi[j][qp](0)) + // diffusion term
                            dt * (Udphi_j) * phi[i][qp]); // convection term

                    Kwp(i, j) += JxW[qp]*(-dt * dphi[i][qp](2) * phi[j][qp]);

                    // G.4: p row
                    Kpu(i, j) += JxW[qp] * (dt * phi[i][qp] * dphi[j][qp](0)); // X divergent operator
                    Kpv(i, j) += JxW[qp] * (dt * phi[i][qp] * dphi[j][qp](1)); // Y divergent operator
                    Kpw(i, j) += JxW[qp] * (dt * phi[i][qp] * dphi[j][qp](2)); // Z divergent operator


                    // Stabilization Contribution
                    // S1. u row
                    Kuu(i, j) += JxW[qp] * tau_m * (Udphi_i * phi[j][qp] + // mass-matrix term
                            dt * DphiUiUjDphj); // advection term
                    Kup(i, j) += JxW[qp] * (dt * tau_m * Udphi_i * dphi[j][qp](0)); // X pressure gradient term

                    //S2. v row
                    Kvv(i, j) += JxW[qp] * tau_m * (Udphi_i * phi[j][qp] + // mass-matrix term
                            dt * DphiUiUjDphj); // advection term
                    Kvp(i, j) += JxW[qp] * (dt * tau_m * Udphi_i * dphi[j][qp](1)); // X pressure gradient term

                    //S3. w row
                    Kww(i, j) += JxW[qp] * tau_m * ((Udphi_i) * phi[j][qp] + // mass-matrix term
                            dt * DphiUiUjDphj); // advection term
                    Kwp(i, j) += JxW[qp] * (dt * tau_m * Udphi_i * dphi[j][qp](2)); // Z pressure gradient term

                    // PSPG
                    // ----
                    //P.1) P row
                    Kpu(i, j) += JxW[qp] * (tau_m * (dphi[i][qp](0) * phi[j][qp] + // mass-matrix term
                            dt * dphi[i][qp](0) * Udphi_j)); // advection term
                    Kpv(i, j) += JxW[qp] * (tau_m * (dphi[i][qp](1) * phi[j][qp] + // mass-matrix term
                            dt * dphi[i][qp](1) * Udphi_j)); // advection term
                    Kpw(i, j) += JxW[qp] * (tau_m * (dphi[i][qp](2) * phi[j][qp] + // mass-matrix term
                            dt * dphi[i][qp](2) * Udphi_j)); // advection term
                    Kpp(i, j) += JxW[qp] * dt * tau_m * (dphi[i][qp] * dphi[j][qp]); // gradient operator

                    // LSIC
                    //L.1) U row
                    Kuu(i, j) += JxW[qp] * tau_c * dt * dphi[i][qp](0) * dphi[j][qp](0); // XX divergent operator
                    Kuv(i, j) += JxW[qp] * tau_c * dt * dphi[i][qp](0) * dphi[j][qp](1); // XY divergent operator
                    Kuw(i, j) += JxW[qp] * tau_c * dt * dphi[i][qp](0) * dphi[j][qp](2); // XZ divergent operator

                    //L.2) V row
                    Kvu(i, j) += JxW[qp] * tau_c * dt * dphi[i][qp](1) * dphi[j][qp](0); // YX divergent operator
                    Kvv(i, j) += JxW[qp] * tau_c * dt * dphi[i][qp](1) * dphi[j][qp](1); // YY divergent operator
                    Kvw(i, j) += JxW[qp] * tau_c * dt * dphi[i][qp](1) * dphi[j][qp](2); // YZ divergent operator

                    //L.3) W row
                    Kwu(i, j) += JxW[qp] * tau_c * dt * dphi[i][qp](2) * dphi[j][qp](0); // ZX divergent operator
                    Kwv(i, j) += JxW[qp] * tau_c * dt * dphi[i][qp](2) * dphi[j][qp](1); // ZY divergent operator
                    Kww(i, j) += JxW[qp] * tau_c * dt * dphi[i][qp](2) * dphi[j][qp](2); // ZZ divergent operator

                }

            }

        } // end of the quadrature point qp-loop

        dof_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

        navier_stokes_system.matrix->add_matrix(Ke, dof_indices);
        navier_stokes_system.rhs->add_vector(Fe, dof_indices);

    } // end of element loop

    perf_log->stop_event("Assembly", "Flow");
    perf_log->restart_event("Solver", "Flow");
    // That's it.
    return;
}


// The matrix assembly function to be called at each time step to
// prepare for the linear solve.

void SedimentationFlow::assemble2D() 
{
    
    PerfLog* perf_log = es.parameters.get<PerfLog*>("PerfLog");
    perf_log->pause_event("Solver", "Flow");
    perf_log->start_event("Assembly", "Flow");

    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the Stokes system object.
    TransientLinearImplicitSystem & navier_stokes_system =
            es.get_system<TransientLinearImplicitSystem> ("flow");

    TransientLinearImplicitSystem & transport_system =
            es.get_system<TransientLinearImplicitSystem> ("transport");



    // Numeric ids corresponding to each variable in the system
    const unsigned int u_var = navier_stokes_system.variable_number("u");
    const unsigned int v_var = navier_stokes_system.variable_number("v");
    const unsigned int p_var = navier_stokes_system.variable_number("p");
    const unsigned int s_var = transport_system.variable_number("s");


#ifdef MESH_MOVIMENT
    LinearImplicitSystem & mesh_system =
            es.get_system<LinearImplicitSystem> ("MeshMoving");
    const unsigned int disp_var = mesh_system.variable_number("disp-z");
#endif


    FEType fe_vel_type = navier_stokes_system.variable_type(u_var);

    UniquePtr<FEBase> fe_vel(FEBase::build(dim, fe_vel_type));
    QGauss qrule(dim, fe_vel_type.default_quadrature_order());


    // Tell the finite element objects to use our quadrature rule.
    fe_vel->attach_quadrature_rule(&qrule);
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real>& JxW = fe_vel->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();


    // A reference to the \p DofMap object for this system.  The \p DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the \p DofMap
    // in future examples.
    const DofMap & dof_map = navier_stokes_system.get_dof_map();
    const DofMap & dof_map_sed = transport_system.get_dof_map();

#ifdef MESH_MOVIMENT
    const DofMap & dof_map_mesh = mesh_system.get_dof_map();
    std::vector<dof_id_type> dof_indices_mesh;
#endif
    
    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    DenseSubMatrix<Number>
            Kuu(Ke), Kuv(Ke), Kup(Ke),
            Kvu(Ke), Kvv(Ke), Kvp(Ke),
            Kpu(Ke), Kpv(Ke), Kpp(Ke);

    DenseSubVector<Number>
            Fu(Fe),
            Fv(Fe),
            Fp(Fe);


    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_u;
    std::vector<dof_id_type> dof_indices_v;
    std::vector<dof_id_type> dof_indices_p;
    std::vector<dof_id_type> dof_indices_s;
    


    const Real dt = es.parameters.get<Real>("dt");
    const Real dt_stab = es.parameters.get<Real>("dt_stab");
    Real Reynolds = es.parameters.get<Real>("Reynolds");
    Real ex = es.parameters.get<Real>("ex");
    Real ey = es.parameters.get<Real>("ey");
    const Real uRe = 1.0 / Reynolds;

    const NumberVectorValue norm(ex, ey);
    //

    //const Real one_dt = 1.0 / dt;

    // Now we will loop over all the elements in the mesh that
    // live on the local processor.
    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for (; el != end_el; ++el) {
        // Store a pointer to the element we are currently
        // working on.  This allows for nicer syntax later.
        const Elem* elem = *el;

        // Compute SUPG stabilization parameters:
        // The carateristic height of the element
        //const Real volume      = elem->volume();
        //const Real h_caract    = 2.0*pow(volume*invPi,0.5);
        //const Real aux1        = 9.0*pow(4.0*uRe/(h_caract*h_caract),2.0)+ 4.0/(dt*dt);


        dof_map.dof_indices(elem, dof_indices);
        dof_map.dof_indices(elem, dof_indices_u, u_var);
        dof_map.dof_indices(elem, dof_indices_v, v_var);
        dof_map.dof_indices(elem, dof_indices_p, p_var);
        dof_map_sed.dof_indices(elem, dof_indices_s, s_var);

#ifdef MESH_MOVIMENT      
        dof_map_mesh.dof_indices(elem, dof_indices_mesh, disp_var);
#endif

        const unsigned int n_dofs = dof_indices.size();
        const unsigned int n_u_dofs = dof_indices_u.size();
        const unsigned int n_v_dofs = dof_indices_v.size();
        const unsigned int n_p_dofs = dof_indices_p.size();

        // Compute the element-specific data for the current
        // element.
        fe_vel->reinit(elem);


        // Zero the element matrix and right-hand side before
        // summing them.
        Ke.resize(n_dofs, n_dofs);
        Fe.resize(n_dofs);

        //
        // Similarly, the \p DenseSubVector.reposition () member
        // takes the (row_offset, row_size)
        Kuu.reposition(u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuv.reposition(u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
        Kup.reposition(u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);

        Kvu.reposition(v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
        Kvv.reposition(v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
        Kvp.reposition(v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);

        Kpu.reposition(p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
        Kpv.reposition(p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
        Kpp.reposition(p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);

        Fu.reposition(u_var*n_u_dofs, n_u_dofs);
        Fv.reposition(v_var*n_u_dofs, n_v_dofs);
        Fp.reposition(p_var*n_u_dofs, n_p_dofs);

        // Now we will build the element matrix and right-hand-side.
        // Constructing the RHS requires the solution and its
        // gradient from the previous timestep.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
            // Values to hold the solution & its gradient at the previous timestep.
            Number u = 0., u_old = 0.;
            Number v = 0., v_old = 0.;
            Number c = 0., fg = 0.;
            Number u_mesh = 0.0, v_mesh = 0.0;
            Gradient grad_u, grad_v, grad_p;
            Point f;
            f.zero();

            RealGradient g = compute_g(fe_vel.get(), 2, qp);
            RealTensor G   = compute_G(fe_vel.get(), 2, qp);


            for (unsigned int l = 0; l < n_u_dofs; l++) {
                // From the old timestep:
                u_old += phi[l][qp] * navier_stokes_system.old_solution(dof_indices_u[l]);
                v_old += phi[l][qp] * navier_stokes_system.old_solution(dof_indices_v[l]);

#ifdef MESH_MOVIMENT
                double v = one_dt * mesh_system.current_solution(dof_indices_mesh[l]);
                v_mesh += phi[l][qp] * v;
#endif              

                // From the previous Newton iterate:
                u += phi[l][qp] * navier_stokes_system.current_solution(dof_indices_u[l]);
                v += phi[l][qp] * navier_stokes_system.current_solution(dof_indices_v[l]);

                c += phi[l][qp] * transport_system.current_solution(dof_indices_s[l]);

                grad_u.add_scaled(dphi[l][qp], navier_stokes_system.current_solution(dof_indices_u[l]));
                grad_v.add_scaled(dphi[l][qp], navier_stokes_system.current_solution(dof_indices_v[l]));
                grad_p.add_scaled(dphi[l][qp], navier_stokes_system.current_solution(dof_indices_p[l]));

            }

            f = norm*c;

            // Definitions for convenience.  It is sometimes simpler to do a
            // dot product if you have the full vector at your disposal.
            const NumberVectorValue U_old(u_old, v_old);
            const NumberVectorValue U(u - u_mesh, v - v_mesh);


            // Residual VMS term
            Real r_m_x = (u - u_old) / dt + (U * grad_u) + grad_p(0) - f(0);
            Real r_m_y = (v - v_old) / dt + (U * grad_v) + grad_p(1) - f(1);

            // SUPG/PSPG/LSIC parameters
            //Number unorm = U.size();

            Real tau_m = compute_tau_M(g, G, U, uRe, dt, dt_stab);
           
            const NumberVectorValue Uvms(u - tau_m*r_m_x - u_mesh, v - tau_m * r_m_y - v_mesh);
            tau_m     = compute_tau_M(g, G, Uvms, uRe, dt, dt_stab); 
            Real tau_c = compute_tau_C(g, tau_m);
            

            for (unsigned int i = 0; i < n_u_dofs; i++) {

                const Number Udphi_i = Uvms * dphi[i][qp];

                Fu(i) += JxW[qp]*(u_old * (phi[i][qp] + tau_m * Udphi_i) + // mass-matrix term
                        dt * f(0)*(phi[i][qp] + tau_m * Udphi_i)); // SUPG term

                Fv(i) += JxW[qp]*(v_old * (phi[i][qp] + tau_m * Udphi_i) + // mass-matrix term
                        dt * f(1)*(phi[i][qp] + tau_m * Udphi_i)); // SUPG term


                Fp(i) += JxW[qp] * tau_m * (U * dphi[i][qp]
                        + dt * (dphi[i][qp] * norm) * fg
                        );

                // Matrix contributions for the uu and vv couplings.
                for (unsigned int j = 0; j < n_u_dofs; j++) 
                {

                    const Number Udphi_j = Uvms * dphi[j][qp];

                    Number DphiUiUjDphj = 0.0;
                    //DphiUiUjDphj = (U*dphi[i][qp])*Udphi_j;


                    for (int d = 0; d < dim; d++) {
                        DphiUiUjDphj += (dphi[i][qp](0) * Uvms(0) * Uvms(d) +
                                dphi[i][qp](1) * Uvms(1) * Uvms(d)) * dphi[j][qp](d);
                    }

                    // Galerkin contribution
                    Kuu(i, j) += JxW[qp]*(phi[i][qp] * phi[j][qp] + // mass matrix term
                            dt * uRe * (2.0 * dphi[i][qp](0) * dphi[j][qp](0) +
                            dphi[i][qp](1) * dphi[j][qp](1)) + // xx diffusion term
                            dt * Udphi_j * phi[i][qp]); // convection term

                    Kuv(i, j) += JxW[qp] * dt * uRe * (dphi[i][qp](1) * dphi[j][qp](0)); // XY diffusion term
                    Kup(i, j) += JxW[qp]*(-dt * phi[j][qp] * dphi[i][qp](0));


                    Kvv(i, j) += JxW[qp]*(phi[i][qp] * phi[j][qp] + // mass matrix term
                            dt * uRe * (2.0 * dphi[i][qp](1) * dphi[j][qp](1) +
                            dphi[i][qp](0) * dphi[j][qp](0)) + // diffusion term
                            dt * (Udphi_j) * phi[i][qp]); // convection term
                    Kvu(i, j) += JxW[qp] * dt * uRe * (dphi[i][qp](0) * dphi[j][qp](1)); // XY diffusion term

                    Kvp(i, j) += JxW[qp]*(-dt * phi[j][qp] * dphi[i][qp](1));



                    Kpu(i, j) += JxW[qp] * (dt * phi[i][qp] * dphi[j][qp](0)); // X divergent operator
                    Kpv(i, j) += JxW[qp] * (dt * phi[i][qp] * dphi[j][qp](1)); // Y divergent operator


                    // SUPG Contribution
                    Kuu(i, j) += JxW[qp] * tau_m * ((Udphi_i) * phi[j][qp] + // mass-matrix term
                            dt * DphiUiUjDphj); // advection term
                    Kup(i, j) += JxW[qp] * (dt * tau_m * Udphi_i * dphi[j][qp](0)); // X pressure gradient term


                    Kvv(i, j) += JxW[qp] * tau_m * ((Udphi_i) * phi[j][qp] + // mass-matrix term
                            dt * DphiUiUjDphj); // advection term
                    Kvp(i, j) += JxW[qp] * (dt * tau_m * Udphi_i * dphi[j][qp](1)); // X pressure gradient term


                    //  Stabilization terms
                    // ----
                    //P.1) P row
                    Kpu(i, j) += JxW[qp] * (tau_m * (dphi[i][qp](0) * phi[j][qp] + // mass-matrix term
                            dt * dphi[i][qp](0) * (U * dphi[j][qp]))); // advection term
                    Kpv(i, j) += JxW[qp] * (tau_m * (dphi[i][qp](1) * phi[j][qp] + // mass-matrix term
                            dt * dphi[i][qp](1) * (U * dphi[j][qp]))); // advection term

                    Kpp(i, j) += JxW[qp] * dt * tau_m * (dphi[i][qp] * dphi[j][qp]); // gradient operator

                    // LSIC
                    // L.1) U row
                    Kuu(i, j) += JxW[qp] * tau_c * dt * dphi[i][qp](0) * dphi[j][qp](0); // XX divergent operator
                    Kuv(i, j) += JxW[qp] * tau_c * dt * dphi[i][qp](0) * dphi[j][qp](1); // XY divergent operator

                    //L.2) V row
                    Kvu(i, j) += JxW[qp] * tau_c * dt * dphi[i][qp](1) * dphi[j][qp](0); // YX divergent operator
                    Kvv(i, j) += JxW[qp] * tau_c * dt * dphi[i][qp](1) * dphi[j][qp](1); // YY divergent operator

                }

            }

        } // end of the quadrature point qp-loop

        dof_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

        navier_stokes_system.matrix->add_matrix(Ke, dof_indices);
        navier_stokes_system.rhs->add_vector(Fe, dof_indices);

    } // end of element loop

    perf_log->stop_event("Assembly", "Flow");
    perf_log->restart_event("Solver", "Flow");

    // That's it.
    return;
}
