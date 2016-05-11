#include "sedimentation_flow.h"

#include "libmesh/zero_function.h"

#include "define.h"


// Define a wrapper for exact_solution that will be needed below
void inlet_wrapper (DenseVector<Number>& output,
                             const Point& p,
                             const Real)
{
  output(0) = 0.5;
  output(1) = 0.0;
  output(2) = 0.0;
}

 void SedimentationFlow::set_gravitational_current_on()
 {
      this->boussinesq = 1.0;
 }

void SedimentationFlow::setup()
{
   MeshBase& mesh = this->es.get_mesh();

   this->dim = mesh.mesh_dimension();

    // Add a transient system to the EquationSystems
  // object named "Convection-Diffusion".
  TransientLinearImplicitSystem & flow_system = this->es.add_system<TransientLinearImplicitSystem> ("flow");

  flow_system.attach_assemble_object(*this);

  unsigned int u_var = flow_system.add_variable("u", FIRST);
  unsigned int v_var = flow_system.add_variable("v", FIRST);
  unsigned int w_var;
  if(dim == 3)
      w_var = flow_system.add_variable("w", FIRST);

  unsigned int p_var = flow_system.add_variable("p", FIRST);

  std::set<boundary_id_type> noslip;
  std::set<boundary_id_type> slipy;
  std::set<boundary_id_type> slipx;
  std::set<boundary_id_type> slipz;
  std::set<boundary_id_type> pnull;
  this->boundary_id[DEP]   = BOUNDARY_DEPOSITION;
  this->boundary_id[PNULL] = BOUNDARY_PNULL;

  if(this->dim == 3) {

      this->boundary_id[MIN_X] = 4;
      this->boundary_id[MIN_Y] = 1;
      this->boundary_id[MIN_Z] = 0;
      this->boundary_id[MAX_X] = 2;
      this->boundary_id[MAX_Y] = 3;
      this->boundary_id[MAX_Z] = 5;

      noslip.insert(this->boundary_id[MIN_Z]);
      //noslip.insert(this->boundary_id[MAX_X]);
      noslip.insert(this->boundary_id[MAX_Z]);

      slipy.insert(this->boundary_id[MIN_Y]);
      slipy.insert(this->boundary_id[MAX_Y]);
      slipx.insert(this->boundary_id[MIN_X]);
      slipx.insert(this->boundary_id[MAX_X]);
      pnull.insert(this->boundary_id[PNULL]);

  }
  else if (dim == 2)
  {
      this->boundary_id[MIN_X] = 3;
      this->boundary_id[MIN_Y] = 0;
      this->boundary_id[MAX_X] = 1;
      this->boundary_id[MAX_Y] = 2;

      noslip.insert(this->boundary_id[MIN_Y]);
      noslip.insert(this->boundary_id[MAX_Y]);
      //noslip.insert(this->boundary_id[MAX_Y]);
      //slipy.insert(this->boundary_id[MAX_Y]);
      slipx.insert(this->boundary_id[MIN_X]);
      slipx.insert(this->boundary_id[MAX_X]);

      pnull.insert(this->boundary_id[PNULL]);

  }


  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();

  int boundary_dep =  (this->dim > 2)?this->boundary_id[MIN_Z]:this->boundary_id[MIN_Y];
  int boundary_top = (this->dim > 2)?this->boundary_id[MAX_Z]:this->boundary_id[MAX_Y];

  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
   for ( ; el != end_el; ++el)
   {
              const Elem* elem = *el;

              unsigned int side_dep     = 0;
              unsigned int side_max_x   = 0;
              unsigned int side_top     = 0;
              bool found_side_dep       = false;
              bool found_side_max_x     = false;
              bool found_side_top       = false;

              for(unsigned int side=0; side<elem->n_sides(); side++)
                {
                  if( mesh.get_boundary_info().has_boundary_id(elem, side, boundary_dep))
                    {
                      side_dep = side;
                      found_side_dep = true;
                    }

                    if( mesh.get_boundary_info().has_boundary_id(elem, side, boundary_top))
                    {
                      side_top = side;
                      found_side_top = true;
                    }

                    if( mesh.get_boundary_info().has_boundary_id(elem, side, this->boundary_id[MAX_X]))
                    {
                      side_max_x = side;
                      found_side_max_x = true;
                    }


                }
                if(found_side_dep)
                {
                  for(unsigned int n=0; n<elem->n_nodes(); n++)
                    {
                      if (elem->is_node_on_side(n, side_dep))
                        {
                          mesh.get_boundary_info().add_node(elem->get_node(n), this->boundary_id[DEP]);
                        }
                    }
                }
                if(found_side_max_x && found_side_top)
                {
                    for(unsigned int n=0; n<elem->n_nodes(); n++)
                    {
                      if (elem->is_node_on_side(n, side_max_x) &&
                          elem->is_node_on_side(n, side_top  )
                          )
                        {
                          mesh.get_boundary_info().add_node(elem->get_node(n), this->boundary_id[PNULL]);
                        }
                    }
                }

    }

  std::vector<unsigned int> velocity_var(2);
  velocity_var[0] = u_var;
  velocity_var[1] = v_var;
  if(this->dim == 3)
     velocity_var.push_back(w_var);


  std::vector<unsigned int> velocity_y(1);
  velocity_y[0]    = v_var;

  std::vector<unsigned int> velocity_x(1);
  velocity_x[0]    = u_var;

  std::vector<unsigned int> pressure(1);
  pressure[0] = p_var;

  ZeroFunction<Number> zero;

  DirichletBoundary dirichlet_noslip(noslip,velocity_var, &zero);
  DirichletBoundary dirichlet_slipy (slipy ,velocity_y  , &zero);
  DirichletBoundary dirichlet_slipx (slipx ,velocity_x  , &zero);
  DirichletBoundary dirichlet_pnull (pnull ,pressure    , &zero);

  flow_system.get_dof_map().add_dirichlet_boundary(dirichlet_noslip);
  flow_system.get_dof_map().add_dirichlet_boundary(dirichlet_slipx);
  flow_system.get_dof_map().add_dirichlet_boundary(dirichlet_slipy);
  flow_system.get_dof_map().add_dirichlet_boundary(dirichlet_pnull);

}


void SedimentationFlow::assemble()
{
    switch(this->dim){
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
void SedimentationFlow::assemble3D()
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "flow");

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Stokes system object.
  TransientLinearImplicitSystem & navier_stokes_system =
    es.get_system<TransientLinearImplicitSystem> ("flow");

  TransientLinearImplicitSystem & transport_system =
    es.get_system<TransientLinearImplicitSystem> ("sediment");

  // Numeric ids corresponding to each variable in the system
  const unsigned int u_var = navier_stokes_system.variable_number ("u");
  const unsigned int v_var = navier_stokes_system.variable_number ("v");
  const unsigned int w_var = navier_stokes_system.variable_number ("w");
  const unsigned int p_var = navier_stokes_system.variable_number ("p");

  const unsigned int s_var = transport_system.variable_number("s");

  FEType fe_vel_type = navier_stokes_system.variable_type(u_var);

  UniquePtr<FEBase> fe_vel   (FEBase::build(dim, fe_vel_type));
  QGauss qrule (dim, fe_vel_type.default_quadrature_order());


  // Tell the finite element objects to use our quadrature rule.
  fe_vel->attach_quadrature_rule (&qrule);
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
  const DofMap & dof_map     = navier_stokes_system.get_dof_map();
  const DofMap & dof_map_sed = transport_system.get_dof_map();

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



  const Real dt       = es.parameters.get<Real>("dt");

  Real Reynolds      = es.parameters.get<Real>("Reynolds");
  Real ex            = es.parameters.get<Real>("ex");
  Real ey            = es.parameters.get<Real>("ey");
  Real ez            = es.parameters.get<Real>("ez");
  const Real uRe      = 1.0/Reynolds;

  const NumberVectorValue norm(ex,ey,ez);
  //
  const Real one3  = 1.0/3.0;
  const Real invPi = 1.0/libMesh::pi;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      // Compute SUPG stabilization parameters:
      // The carateristic height of the element
      const Real volume      = elem->volume();
      const Real h_caract    = pow(6.0*volume*invPi, one3);
      const Real aux1        = 9.0*pow(4.0*uRe/(h_caract*h_caract),2.0)+ 4.0/(dt*dt);


      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (elem, dof_indices_w, w_var);
      dof_map.dof_indices (elem, dof_indices_p, p_var);
      dof_map_sed.dof_indices (elem, dof_indices_s, s_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size();
      const unsigned int n_v_dofs = dof_indices_v.size();
      const unsigned int n_w_dofs = dof_indices_w.size();
      const unsigned int n_p_dofs = dof_indices_p.size();

      // Compute the element-specific data for the current
      // element.
      fe_vel->reinit  (elem);


      // Zero the element matrix and right-hand side before
      // summing them.
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      //
      // Similarly, the \p DenseSubVector.reposition () member
      // takes the (row_offset, row_size)
      Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
      Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);
      Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);

      Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
      Kvw.reposition (v_var*n_v_dofs, w_var*n_v_dofs, n_v_dofs, n_w_dofs);
      Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);

      Kwu.reposition (w_var*n_w_dofs, u_var*n_w_dofs, n_w_dofs, n_u_dofs);
      Kwv.reposition (w_var*n_w_dofs, v_var*n_w_dofs, n_w_dofs, n_v_dofs);
      Kww.reposition (w_var*n_w_dofs, w_var*n_w_dofs, n_w_dofs, n_w_dofs);
      Kwp.reposition (w_var*n_w_dofs, p_var*n_w_dofs, n_w_dofs, n_p_dofs);

      Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
      Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
      Kpw.reposition (p_var*n_u_dofs, w_var*n_u_dofs, n_p_dofs, n_w_dofs);
      Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);

      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);
      Fw.reposition (w_var*n_u_dofs, n_w_dofs);
      Fp.reposition (p_var*n_u_dofs, n_p_dofs);

      // Now we will build the element matrix and right-hand-side.
      // Constructing the RHS requires the solution and its
      // gradient from the previous timestep.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Values to hold the solution & its gradient at the previous timestep.
          Number   u = 0., u_old = 0.;
          Number   v = 0., v_old = 0.;
          Number   w = 0., w_old = 0.;
          Number   c = 0., fg = 0.;
          Gradient grad_u, grad_v, grad_w, grad_p;
          Point f;
          f.zero();

          for (unsigned int l=0; l<n_u_dofs; l++)
            {
              // From the old timestep:
              u_old += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_u[l]);
              v_old += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_v[l]);
              w_old += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_w[l]);


              // From the previous Newton iterate:
              u += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_u[l]);
              v += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_v[l]);
              w += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_w[l]);

              c += phi[l][qp]*transport_system.current_solution (dof_indices_s[l]);

              grad_u.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_u[l]));
              grad_v.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_v[l]));
              grad_w.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_w[l]));
              grad_p.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_p[l]));

            }

           fg = c;
           f(0) = norm(0)*fg;
           f(1) = norm(1)*fg;
           f(2) = norm(2)*fg;

          // Definitions for convenience.  It is sometimes simpler to do a
          // dot product if you have the full vector at your disposal.
          const NumberVectorValue U_old (u_old, v_old, w_old);
          const NumberVectorValue U     (u,     v,  w);

          // Residual VMS term
          Real rint_x = (u-u_old)/dt + (U*grad_u) + grad_p(0) - f(0);
          Real rint_y = (v-v_old)/dt + (U*grad_v) + grad_p(1) - f(1);
          Real rint_z = (w-w_old)/dt + (U*grad_w) + grad_p(2) - f(2);


	        // SUPG/PSPG/LSIC parameters
          Number unorm  = U.size();
          Number aux2   = pow(2.0*unorm/h_caract,2.0) + aux1;
          Number tau    = pow(aux2,-0.5);

          const Number tau_s = tau;
          const NumberVectorValue Uvms (u-tau_s*rint_x, v-tau_s*rint_y, w-tau_s*rint_z);

          unorm = Uvms.size();
          aux2  = pow(2.0*unorm/h_caract,2.0) + aux1;
          tau    = pow(aux2,-0.5);
          Number teta   = unorm * h_caract * 0.5;

          // First, an i-loop over the velocity degrees of freedom.
          // We know that n_u_dofs == n_v_dofs so we can compute contributions
          // for both at the same time.
          for (unsigned int i=0; i<n_u_dofs; i++)
            {

             const Number Udphi_i = Uvms*dphi[i][qp];

              Fu(i) += JxW[qp]*(u_old*(phi[i][qp] + tau*Udphi_i) +                  // mass-matrix term
			        dt*f(0)*(phi[i][qp]+ tau*Udphi_i));                 // SUPG term

              Fv(i) += JxW[qp]*(v_old*(phi[i][qp] + tau*Udphi_i) +                  // mass-matrix term
			        dt*f(1)*(phi[i][qp]+ tau*Udphi_i));                 // SUPG term

              Fw(i) += JxW[qp]*(w_old*(phi[i][qp] + tau*Udphi_i) +                  // mass-matrix term
			        dt*f(2)*(phi[i][qp]+ tau*Udphi_i));                 // SUPG term

              Fp(i) += JxW[qp] * tau*(U_old*dphi[i][qp]                             // PSPG term
                                      + dt*(dphi[i][qp]*norm)*fg                    //PSPG buoyancy vector
	                             );

              // Matrix contributions for the uu and vv couplings.
              for (unsigned int j=0; j<n_u_dofs; j++)
              {
		              const Number Udphi_j = Uvms*dphi[j][qp];

	                // Galerkin contribution
                  Kuu(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp] +                     // mass matrix term
                                       dt*uRe*2.0*(dphi[i][qp]*dphi[j][qp]) +      // xx diffusion term
                                       dt*Udphi_j*phi[i][qp]) ;                    // convection term

		              Kuv(i,j) += JxW[qp]*dt*uRe*(dphi[i][qp](1)*dphi[j][qp](0));      // xy diffusion term
                  Kuw(i,j) += JxW[qp]*dt*uRe*(dphi[i][qp](2)*dphi[j][qp](0));      // xz diffusion term
		              Kup(i,j) += JxW[qp]*(-dt*phi[j][qp]*dphi[i][qp](0));


                  Kvv(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp] +                    // mass matrix term
                                       dt*uRe*2.0*(dphi[i][qp]*dphi[j][qp]) +     // diffusion term
                                       dt*(Udphi_j)*phi[i][qp] )   ;              // convection term

                  Kvu(i,j) += JxW[qp]*dt*uRe*(dphi[i][qp](0)*dphi[j][qp](1));     // xy diffusion term
                  Kvw(i,j) += JxW[qp]*dt*uRe*(dphi[i][qp](2)*dphi[j][qp](1));     // xy diffusion term
                  Kvp(i,j) += JxW[qp]*(-dt*phi[j][qp]*dphi[i][qp](1));

                  Kww(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp] +                    // mass matrix term
                                       dt*uRe*2.0*(dphi[i][qp]*dphi[j][qp]) +     // diffusion term
                                       dt*(Udphi_j)*phi[i][qp] );                 // convection term


                  Kwu(i,j) += JxW[qp]*dt*uRe*(dphi[i][qp](0)*dphi[j][qp](2));     // xy diffusion term
                  Kwv(i,j) += JxW[qp]*dt*uRe*(dphi[i][qp](1)*dphi[j][qp](2));     // xy diffusion term
		              Kwp(i,j) += JxW[qp]*(-dt*phi[j][qp]*dphi[i][qp](2));


		              Kpu(i,j) += JxW[qp] * (dt*phi[i][qp]*dphi[j][qp](0));           // X divergent operator
                  Kpv(i,j) += JxW[qp] * (dt*phi[i][qp]*dphi[j][qp](1));           // Y divergent operator
                  Kpw(i,j) += JxW[qp] * (dt*phi[i][qp]*dphi[j][qp](2));           // Z divergent operator

		              // SUPG Contribution
		              Kuu(i,j)+= JxW[qp]*tau*((Udphi_i)*phi[j][qp]  +                 // mass-matrix term
                                          dt*Udphi_i*Udphi_j);                    // advection term
                  Kup(i,j)+= JxW[qp] * (dt*tau*Udphi_i*dphi[j][qp](0));           // X pressure gradient term


		              Kvv(i,j)+= JxW[qp]*tau*((Udphi_i)*phi[j][qp]  +                 // mass-matrix term
                                          dt*Udphi_i*Udphi_j);                    // advection term
                  Kvp(i,j)+= JxW[qp] * (dt*tau*Udphi_i*dphi[j][qp](1));           // X pressure gradient term

                  Kww(i,j)+= JxW[qp]*tau*((Udphi_i)*phi[j][qp]  +                 // mass-matrix term
                                          dt*Udphi_i*Udphi_j);                    // advection term
                  Kwp(i,j)+= JxW[qp] * (dt*tau*Udphi_i*dphi[j][qp](2));           // Z pressure gradient term

		              // PSPG
                  // ----
                  //P.1) P row
                  Kpu(i,j)+= JxW[qp] * (tau * (dphi[i][qp](0)*phi[j][qp] +          // mass-matrix term
                                        dt*dphi[i][qp](0) * Udphi_j ));             // advection term
                  Kpv(i,j)+= JxW[qp] * (tau * (dphi[i][qp](1)*phi[j][qp] +          // mass-matrix term
                                         dt*dphi[i][qp](1) * Udphi_j ));            // advection term
                  Kpw(i,j)+= JxW[qp] * (tau * (dphi[i][qp](2)*phi[j][qp] +          // mass-matrix term
                                         dt*dphi[i][qp](2) * Udphi_j ));            // advection term
                  Kpp(i,j)+= JxW[qp] *  dt*tau*(dphi[i][qp]*dphi[j][qp]);           // gradient operator

                  // LSIC
                  /*
                  //L.1) U row
                  Kuu(i,j)+= JxW[qp] * teta * dt * dphi[i][qp](0) * dphi[j][qp](0); // XX divergent operator
                  Kuv(i,j)+= JxW[qp] * teta * dt * dphi[i][qp](0) * dphi[j][qp](1); // XY divergent operator
                  Kuw(i,j)+= JxW[qp] * teta * dt * dphi[i][qp](0) * dphi[j][qp](2); // XZ divergent operator
                  //L.2) V row
                  Kvu(i,j)+= JxW[qp] * teta * dt * dphi[i][qp](1) * dphi[j][qp](0); // YX divergent operator
                  Kvv(i,j)+= JxW[qp] * teta * dt * dphi[i][qp](1) * dphi[j][qp](1); // YY divergent operator
                  Kvw(i,j)+= JxW[qp] * teta * dt * dphi[i][qp](1) * dphi[j][qp](2); // YZ divergent operator
                  //L.3) W row
                  Kwu(i,j)+= JxW[qp] * teta * dt * dphi[i][qp](2) * dphi[j][qp](0); // ZX divergent operator
                  Kwv(i,j)+= JxW[qp] * teta * dt * dphi[i][qp](2) * dphi[j][qp](1); // ZY divergent operator
                  Kww(i,j)+= JxW[qp] * teta * dt * dphi[i][qp](2) * dphi[j][qp](2); // ZZ divergent operator
		  */
                }

            }

        } // end of the quadrature point qp-loop

      dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      navier_stokes_system.matrix->add_matrix (Ke, dof_indices);
      navier_stokes_system.rhs->add_vector    (Fe, dof_indices);

    } // end of element loop

  // That's it.
  return;
}


// The matrix assembly function to be called at each time step to
// prepare for the linear solve.
void SedimentationFlow::assemble2D()
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "flow");

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Stokes system object.
  TransientLinearImplicitSystem & navier_stokes_system =
    es.get_system<TransientLinearImplicitSystem> ("flow");

  TransientLinearImplicitSystem & transport_system =
    es.get_system<TransientLinearImplicitSystem> ("sediment");

  // Numeric ids corresponding to each variable in the system
  const unsigned int u_var = navier_stokes_system.variable_number ("u");
  const unsigned int v_var = navier_stokes_system.variable_number ("v");
  const unsigned int p_var = navier_stokes_system.variable_number ("p");

  const unsigned int s_var = transport_system.variable_number("s");

  FEType fe_vel_type = navier_stokes_system.variable_type(u_var);

  UniquePtr<FEBase> fe_vel   (FEBase::build(dim, fe_vel_type));
  QGauss qrule (dim, fe_vel_type.default_quadrature_order());


  // Tell the finite element objects to use our quadrature rule.
  fe_vel->attach_quadrature_rule (&qrule);
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
  const DofMap & dof_map     = navier_stokes_system.get_dof_map();
  const DofMap & dof_map_sed = transport_system.get_dof_map();

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


  const Real dt       = es.parameters.get<Real>("dt");

  Real Reynolds      = es.parameters.get<Real>("Reynolds");
  Real ex            = es.parameters.get<Real>("ex");
  Real ey            = es.parameters.get<Real>("ey");
  Real ez            = es.parameters.get<Real>("ez");
  const Real uRe      = 1.0/Reynolds;

  const NumberVectorValue norm(ex,ez);
  //
  const Real one3  = 1.0/3.0;
  const Real invPi = 1.0/libMesh::pi;
  const Real theta = 0.5;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      // Compute SUPG stabilization parameters:
      // The carateristic height of the element
      const Real volume      = elem->volume();
      const Real h_caract    = 2.0*pow(volume*invPi,0.5);
      const Real aux1        = 9.0*pow(4.0*uRe/(h_caract*h_caract),2.0)+ 4.0/(dt*dt);


      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (elem, dof_indices_p, p_var);
      dof_map_sed.dof_indices (elem, dof_indices_s, s_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size();
      const unsigned int n_v_dofs = dof_indices_v.size();
      const unsigned int n_p_dofs = dof_indices_p.size();

      // Compute the element-specific data for the current
      // element.
      fe_vel->reinit  (elem);


      // Zero the element matrix and right-hand side before
      // summing them.
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      //
      // Similarly, the \p DenseSubVector.reposition () member
      // takes the (row_offset, row_size)
      Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
      Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);

      Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
      Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);

      Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
      Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
      Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);

      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);
      Fp.reposition (p_var*n_u_dofs, n_p_dofs);

      // Now we will build the element matrix and right-hand-side.
      // Constructing the RHS requires the solution and its
      // gradient from the previous timestep.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {
          // Values to hold the solution & its gradient at the previous timestep.
          Number   u = 0., u_old = 0.;
          Number   v = 0., v_old = 0.;
          Number   c = 0., fg = 0.;
          Gradient grad_u, grad_v, grad_p;
          Point f;
          f.zero();

          for (unsigned int l=0; l<n_u_dofs; l++)
          {
              // From the old timestep:
              u_old += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_u[l]);
              v_old += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_v[l]);

              // From the previous Newton iterate:
              u += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_u[l]);
              v += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_v[l]);

              c += phi[l][qp]*transport_system.current_solution (dof_indices_s[l]);

              grad_u.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_u[l]));
              grad_v.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_v[l]));
              grad_p.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_p[l]));

          }

          fg = c;
          f(0) = norm(0)*fg;
          f(1) = norm(1)*fg;

          // Definitions for convenience.  It is sometimes simpler to do a
          // dot product if you have the full vector at your disposal.
          const NumberVectorValue U_old (u_old, v_old);
          const NumberVectorValue U     (u,     v);

          // Residual VMS term
          Real rint_x = (u-u_old)/dt + (U*grad_u) + grad_p(0) - f(0);
          Real rint_y = (v-v_old)/dt + (U*grad_v) + grad_p(1) - f(1);


	        // SUPG/PSPG/LSIC parameters

          Number unorm  = U.size();
          Number aux2   = pow(2.0*unorm/h_caract,2.0) + aux1;
          Number tau    = pow(aux2,-0.5);

          const Number tau_s = tau;
          const NumberVectorValue Uvms (u-tau_s*rint_x, v-tau_s*rint_y);

          unorm  = Uvms.size();
          aux2   = pow(2.0*unorm/h_caract,2.0) + aux1;
          tau    = pow(aux2,-0.5);
          Number teta  = unorm * h_caract * 0.5;
          teta = 0.0;


          // First, an i-loop over the velocity degrees of freedom.
          // We know that n_u_dofs == n_v_dofs so we can compute contributions
          // for both at the same time.
          for (unsigned int i=0; i<n_u_dofs; i++)
          {

              const Number Udphi_i = Uvms*dphi[i][qp];

              Fu(i) += JxW[qp]*(u_old*(phi[i][qp] + tau*Udphi_i) +              // mass-matrix term
			                          dt*f(0)*(phi[i][qp]+ tau*Udphi_i));             // SUPG term

              Fv(i) += JxW[qp]*(v_old*(phi[i][qp] + tau*Udphi_i) +              // mass-matrix term
			                          dt*f(1)*(phi[i][qp]+ tau*Udphi_i));             // SUPG term


              Fp(i) += JxW[qp] * tau*(U*dphi[i][qp]
                             + dt*(dphi[i][qp]*norm)*fg
                           );                          // PSPG term);

              // Matrix contributions for the uu and vv couplings.
              for (unsigned int j=0; j<n_u_dofs; j++)
              {
                  const Number Udphi_j = Uvms*dphi[j][qp];

	                // Galerkin contribution
                  Kuu(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]                +   // mass matrix term
                                       dt*uRe*2.0*(dphi[i][qp]*dphi[j][qp]) +   // xx diffusion term
                                       dt*Udphi_j*phi[i][qp]) ;                 // convection term

		              Kuv(i,j) += JxW[qp]*dt*uRe*(dphi[i][qp](1)*dphi[j][qp](0));   // xy diffusion term
		              Kup(i,j) += JxW[qp]*(-dt*phi[j][qp]*dphi[i][qp](0));


                  Kvv(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp] +                  // mass matrix term
                                       dt*uRe*2.0*(dphi[i][qp]*dphi[j][qp]) +   // diffusion term
                                       dt*(Udphi_j)*phi[i][qp] )   ;            // convection term
		              Kvu(i,j) += JxW[qp]*dt*uRe*(dphi[i][qp](0)*dphi[j][qp](1));   // xy diffusion term

                  Kvp(i,j) += JxW[qp]*(-dt*phi[j][qp]*dphi[i][qp](1));



		              Kpu(i,j) += JxW[qp] * (dt*phi[i][qp]*dphi[j][qp](0));         // X divergent operator
                  Kpv(i,j) += JxW[qp] * (dt*phi[i][qp]*dphi[j][qp](1));         // Y divergent operator


		              // SUPG Contribution
		              Kuu(i,j)+= JxW[qp]*tau*((Udphi_i)*phi[j][qp]  +               // mass-matrix term
                                          dt*Udphi_i*Udphi_j);                  // advection term
                  Kup(i,j)+= JxW[qp] * (dt*tau*Udphi_i*dphi[j][qp](0));         // X pressure gradient term


		              Kvv(i,j)+= JxW[qp]*tau*((Udphi_i)*phi[j][qp]  +               // mass-matrix term
                                          dt*Udphi_i*Udphi_j);                  // advection term
                  Kvp(i,j)+= JxW[qp] * (dt*tau*Udphi_i*dphi[j][qp](1));         // X pressure gradient term


		              // PSPG
                  // ----
                  //P.1) P row
                  Kpu(i,j)+= JxW[qp] * (tau * (dphi[i][qp](0)*phi[j][qp] +      // mass-matrix term
                                        dt*dphi[i][qp](0) * Udphi_j ));         // advection term
                  Kpv(i,j)+= JxW[qp] * (tau * (dphi[i][qp](1)*phi[j][qp] +      // mass-matrix term
                                         dt*dphi[i][qp](1) * Udphi_j ));        // advection term

                  Kpp(i,j)+= JxW[qp] *  dt*tau*(dphi[i][qp]*dphi[j][qp]);       // gradient operator

                  // LSIC

                  //L.1) U row
                  Kuu(i,j)+= JxW[qp] * teta * dt * dphi[i][qp](0) * dphi[j][qp](0); // XX divergent operator
                  Kuv(i,j)+= JxW[qp] * teta * dt * dphi[i][qp](0) * dphi[j][qp](1); // XY divergent operator

                  //L.2) V row
                  Kvu(i,j)+= JxW[qp] * teta * dt * dphi[i][qp](1) * dphi[j][qp](0); // YX divergent operator
                  Kvv(i,j)+= JxW[qp] * teta * dt * dphi[i][qp](1) * dphi[j][qp](1); // YY divergent operator

                }

            }

      } // end of the quadrature point qp-loop

      dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      navier_stokes_system.matrix->add_matrix (Ke, dof_indices);
      navier_stokes_system.rhs->add_vector    (Fe, dof_indices);

    } // end of element loop

  // That's it.
  return;
}
