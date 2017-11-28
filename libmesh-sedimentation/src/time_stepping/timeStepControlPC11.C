/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   timeStepControlPID.cpp
 * Author: arossa
 * 
 * Created on 29 de Junho de 2016, 16:31
 */

#include <iostream>
#include <sstream>
#include <ostream>
#include <fstream>

#include "timeStepControlPC11.h"

using namespace libMesh;
using namespace std;

timeStepControlPC11::timeStepControlPC11( double dt_init, double dt_min, double dt_max, unsigned int nsa_max, double tol_U, double tol_S, double theta, double alpha, double k_exp, double s_min, double s_max, bool complete):
                                        timeStepControlBase(dt_init, dt_min, dt_max, dt_init, 0, 2, nsa_max, "PC11"),
                                        tol_U(tol_U),
                                        tol_S(tol_S),
                                        rn_1(1.0),
                                        theta(theta),
                                        alpha(alpha),
                                        k_exp(k_exp),
                                        s_min(s_min),
                                        s_max(s_max),
                                        complet_flow_norn(complete)
{
    cout<<"PC11 time-step control will be adopted"<<((complet_flow_norn)?" with complete solution vector norm.\n":" with only velocity solution vector norm.\n"); 
}

timeStepControlPC11::timeStepControlPC11(const timeStepControlPC11& orig):
                                       timeStepControlBase() 
{
    cout<<"contrutor de copia\n";
}

timeStepControlPC11::~timeStepControlPC11() {
}

void  timeStepControlPC11::computeSolutionChangeInTime(EquationSystems & es) {
    
    // solution vector and its variation in respect to last NLS L2 norms for the current time step
    Real r_flow = 1.0e-10, flow_norm, flow_diff_norm;
    Real r_transport = 1.0e-10;
    
    // Get a reference to the Flow system object.
    // PS: There is always a "flow" system 
    TransientLinearImplicitSystem & flow_system_reference =
    es.get_system<TransientLinearImplicitSystem> ("flow");   
    
    // Here I get a copy of the solution vector at the end of the nonlinear solution process.
    // This vector will be used to measure the change of the quantities of interest in time, i.e.,
    // the variation on velocity vector and pressure over two successive solutions over time.
    UniquePtr<NumericVector<Number> > current_flow_soln (flow_system_reference.current_local_solution->clone()); // "current_local_solution" was already updated after solve

    // Compute the difference between the current and last (in time) solutions.
    current_flow_soln->add (-1., *flow_system_reference.old_local_solution);
         
    // Close the vector before computing its norm
    current_flow_soln->close();
    
    if (this->complet_flow_norn) { // to compute the complete solution vector L2 norm (includes pressure)
        // Compute the l2 norm of the flow solution (includes pressure)
        flow_norm = flow_system_reference.solution->l2_norm(); 

        // Compute the l2 norm of the difference between two successive solutions
        flow_diff_norm = current_flow_soln->l2_norm();
    
        //cout<<"Complete Solution L2 norm = "<<flow_norm<<" Complete Diff Solution L2 norm = "<<flow_diff_norm<<endl;
    
        // compute the measure of the change in time for flow
        if (flow_diff_norm>1.0e-10) { // avoid 'e_flow' becomes zero!
            r_flow = flow_diff_norm / flow_norm;
            // divide by tolerance here to compare this value with that floe transport problem
            r_flow /= this->tol_U;
        }
        //cout<<" rflow_complete = "<< r_flow<<endl;
    } else { // to compute the only velocity solution vector L2 norm (excludes pressure)

        // The dimension that we are running
        const unsigned int dim = es.get_mesh().mesh_dimension();

        MPI_Comm comm = MPI_COMM_WORLD;

        double global_norm, global_dif_norm, local_sqr_norm = 0.0;
        double local_sqr_dif_norm = 0.0;
        int ndata = 2;
        double rms[ndata], rms_sum[ndata];

        for (int i=0; i<ndata; i++) {
           rms[i]= 0.0;
           rms_sum[i]= 0.0;
        }

        // counting on only velocity components to compute L2 norm  
        for (numeric_index_type i = flow_system_reference.solution->first_local_index(); i<flow_system_reference.solution->last_local_index()-dim; i+=dim+1) {
            for ( numeric_index_type k = i; k < i+dim; k++ ) {
                local_sqr_norm += ( flow_system_reference.solution->el(k)*flow_system_reference.solution->el(k) );
                local_sqr_dif_norm += (current_flow_soln->el(k)*current_flow_soln->el(k) );
            }
        }

        rms[0] = local_sqr_norm;
        rms[1] = local_sqr_dif_norm;

        MPI_Allreduce(rms, rms_sum, ndata, MPI_DOUBLE, MPI_SUM, comm);

        global_norm = std::pow(rms_sum[0],0.5);
        global_dif_norm = std::pow(rms_sum[1],0.5);

        // compute the measure of the change in time for flow
        if (global_dif_norm>1.0e-10) { // avoid 'r_flow' to become zero!
            r_flow = global_dif_norm / global_norm;
            r_flow /= this->tol_U;
        }
    }    
    
    // Get a reference to the Transport system object.
    // PS: Must be checked cause "transport" system may not be defined
    if ( es.has_system("transport") ) {
        
        //std::cout << " \nChange Flow solution in time = " << r_flow;
        
        TransientLinearImplicitSystem & transport_system_reference =
        es.get_system<TransientLinearImplicitSystem> ("transport");
        
        Real transport_norm, transport_diff_norm;
        
        // Compute the l2 norm of the transport solution
        transport_norm = transport_system_reference.solution->l2_norm();    

        // Here I get a copy of the solution vector at the end of the nonlinear solution process.
        // This vector will be used to measure the change of the quantities of interest in time, i.e.,
        // the variation on scalar variable over two successive solutions over time.
        UniquePtr<NumericVector<Number> > current_transport_soln (transport_system_reference.current_local_solution->clone()); // "current_local_solution" was already updated after solve

        // Compute the difference between the current and last (in time) solutions.
        current_transport_soln->add (-1., *transport_system_reference.old_local_solution);

        // Close the vector before computing its norm
        current_transport_soln->close();

        // Compute the l2 norm of the difference between two successive solutions
        transport_diff_norm = current_transport_soln->l2_norm();
        
        // compute the measure of the change in time for transport
        if (transport_diff_norm>1.0e-10) { // avoid 'e_transport' to become zero!
            r_transport = transport_diff_norm / transport_norm;
            r_transport /= this->tol_S;
        }
        
        //std::cout << " \nChange Transport solution in time = " << r_transport << std::endl;
    }
            
    // compute the maximum measure of the change in time for the systems
    this->rn = std::max(r_flow, r_transport);

    std::cout << "\n Measure of the solution change in current TS = " << this->rn;
    //std::cout << "\n Measure of the solution change in previous TS = " <<this->rn_1;
}

void timeStepControlPC11::checkTimeStepAcceptance(Real dt, int flow_nonlinear_iteractions, int transport_nonlinear_iteractions, bool& accepted)
{
    accepted = true;
    if ( (this->rn > 1.0  || flow_nonlinear_iteractions > this->nsa_max || transport_nonlinear_iteractions > this->nsa_max ) && dt > this->dt_min) {
        cout << "\n Current solution was REJECTED";        
        accepted = false;
        // updating number of rejected time-steps
        this->n_rejected_time_step++;
    }
    
    // check whether the current time-step was accepted because it was already the minimum 
    if (accepted) {
        std::cout << "\n Current solution was ACCEPTED";        
        this->ckeckKeepMinTimeStep(dt, flow_nonlinear_iteractions, transport_nonlinear_iteractions);
        // updating number of accepted time-steps
        this->n_accepted_time_step++;        
    }      
}

void timeStepControlPC11::computeTimeStep(bool accepted, Real time, Real tmax, Real& dt)
{
    
    if(!accepted)
    {
        this->dt_avg -= dt;
       
        // For while next time-step will be reduced by an alpha factor [0,1]
        dt = this->alpha * dt;  
        
        //std::cout << "\n New calculated time step = " << dt << std::endl;
        
        // Assure time-step size not smaller than dt_min!
        dt = max(dt, this->dt_min);        
        
    }else
    {      
        if(this->keep_dt_min) {
            std::cout <<" because time-step is minimum";
            // since current time-step value was accepted only because it is equal the minimum, it's value will be repeated

            //storing last accepted time-step
            this->dt_last = dt;

        } else {
            //std::cout <<"!";
            double exp = 1.0/(this->k_exp+1.0);
            double dt_star = this->theta * pow( this->rn_1/(this->rn*this->rn ),exp) * dt*dt/this->dt_last;
            
            //storing last accepted time-step
            this->dt_last = dt;

            dt = std::min(this->s_max*dt, std::max(this->s_min*dt, dt_star) );

            //std::cout << "\n New calculated time step = " << dt;

            // check time step range
            dt = max(dt,this->dt_min);
            dt = min(dt,this->dt_max);
        }

        // check whether simulation has reached the end to avoid exceed maximum simulation time
        if (fabs(time - tmax) > 1.0e-08)
            if (time+dt>tmax || tmax-(time+dt) < this->dt_min)
                dt = tmax - time;
            
        // storing error measure
        this->storeSolutionChangeinTime();     
            
    }

    std::cout << "\n Next adopted time step = " << dt <<". Var = "<<(dt-this->dt_last)/this->dt_last*100 <<"%"<< std::endl;
}

void timeStepControlPC11::storeSolutionChangeinTime() {
    this->rn_1 = this->rn;
}

void timeStepControlPC11::printSelf (ostream& os, const char* indent) const {
    // Write performance parameters to output file
    os 	<< "\n"
        << indent <<"PC11 Time-Step Control Performance:"
        << "\nTotal number of time-steps: "
        << this->n_accepted_time_step + this->n_rejected_time_step
        << "\nNumber of accepted time-steps: "
        << this->n_accepted_time_step               
        << "\nNumber of rejected time-steps: "
        << this->n_rejected_time_step            
        << "\nRate Rejected/Total time-steps: "
        << 100.0 * this->n_rejected_time_step/(this->n_accepted_time_step + this->n_rejected_time_step)<<"%"
        << "\nEffective time-step size average: "
        << this->dt_avg/this->n_accepted_time_step << std::endl;    
}

void timeStepControlPC11::ckeckKeepMinTimeStep( Real dt, int flow_nonlinear_iteractions, int transport_nonlinear_iteractions) 
{
    if ( (this->rn > 1.0 || flow_nonlinear_iteractions > (this->nsa_max) || transport_nonlinear_iteractions > (this->nsa_max) ) && abs(dt - this->dt_min) < 1.0e-10 )
        this->keep_dt_min = true;
    else
        this->keep_dt_min = false;       
}
    