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

#include "timeStepControlPID.h"

using namespace libMesh;
using namespace std;

timeStepControlPID::timeStepControlPID( double dt_init, double dt_min, double dt_max, unsigned int nsa_max, double tol_U, double tol_S, double kp, double ki, double kd):
                                        timeStepControlBase(dt_init, dt_min, dt_max, dt_init, 0, 3, nsa_max, "PID"),
                                        tol_U(tol_U),
                                        tol_S(tol_S),
                                        kp(kp),
                                        ki(ki),
                                        kd(kd),
                                        en_1(1.0),
                                        en_2(1.0),
                                        dt_prev(dt_init)
{
    cout<<"PID time-step control will be adopted\n"; 
}

timeStepControlPID::timeStepControlPID(const timeStepControlPID& orig):
                                       timeStepControlBase() 
{
    cout<<"contrutor de copia\n";
}

timeStepControlPID::~timeStepControlPID() {
}

void  timeStepControlPID::computeSolutionChangeInTime(EquationSystems & es) {
    
    // solution vector and its variation in respect to last NLS L2 norms for the current time step
    Real e_flow = 1.0e-10, flow_norm, flow_diff_norm;
    Real e_transport = 1.0e-10;
    
    // Get a reference to the Flow system object.
    // PS: There is always a "flow" system 
    TransientLinearImplicitSystem & flow_system_reference =
    es.get_system<TransientLinearImplicitSystem> ("flow");
    
    
    const unsigned int dim = es.get_mesh().mesh_dimension();
    
    // Here I get a copy of the solution vector at the end of the nonlinear solution process.
    // This vector will be used to measure the change of the quantities of interest in time, i.e.,
    // the variation on velocity vector and pressure over two successive solutions over time.
    UniquePtr<NumericVector<Number> > current_flow_soln (flow_system_reference.current_local_solution->clone()); // "current_local_solution" was already updated after solve

    // Compute the difference between the current and last (in time) solutions.
    current_flow_soln->add (-1., *flow_system_reference.old_local_solution);
         
    // Close the vector before computing its norm
    current_flow_soln->close();

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
     
    global_norm     = std::pow(rms_sum[0],0.5);
    global_dif_norm = std::pow(rms_sum[1],0.5);

    // compute the measure of the change in time for flow
    if (global_dif_norm>1.0e-10) { // avoid 'e_flow' to become zero!
        e_flow = global_dif_norm / global_norm;
        e_flow /= this->tol_U;
    }
     
    
    // Get a reference to the Transport system object.
    // PS: Must be checked cause "transport" system may not be defined
    if ( es.has_system("transport") ) {
        
        //std::cout << " \nChange Flow solution in time = " << e_flow;
        
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
            e_transport = transport_diff_norm / transport_norm;
            e_transport /= this->tol_S;
        }
        
        //std::cout << " \nChange Transport solution in time = " << e_transport << std::endl;        
    }
            
    // compute the maximum measure of the change in time for the systems
    this->en = std::max(e_flow, e_transport);

    std::cout << "\n Measure of the solution change in current TS = " << this->en;
    //std::cout << "\n Measure of the solution change in previous TS = " << this->en_1;
    //std::cout << "\n Measure of the solution change in older    TS = " << this->en_2;
}

void timeStepControlPID::checkTimeStepAcceptance(Real dt, int flow_nonlinear_iteractions, int transport_nonlinear_iteractions, bool& accepted)
{
    accepted = true;
    if ( (this->en > 1.0 || flow_nonlinear_iteractions > this->nsa_max || transport_nonlinear_iteractions > this->nsa_max ) && dt > this->dt_min) {
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

void timeStepControlPID::computeTimeStep(bool accepted, Real time, Real tmax, Real &dt)
{
    if(!accepted) 
    {
        this->dt_avg -= dt;
        Real factor = 1.0/this->en;

        if (factor > 0.8)
            factor = 0.8;
        
        dt = factor*dt;
        //std::cout << " New calculated time step = "<< dt<<endl;
                
        // check time-step range
        dt = max(dt,this->dt_min);
        
        // to compute next time-step
        this->dt_prev = dt*dt/this->dt_prev;
            
    } else 
    {   //current time-step is accepted
        //storing last accepted time-step
        this->dt_last = dt;
                
        if(this->keep_dt_min) {
            std::cout <<" because time-step is minimum";
            // since current time-step value was accepted only because it is equal the minimum, it's value will be repeated
        } else {
            //std::cout <<"!";
            dt = pow(this->en_1/this->en,this->kp)*pow(1.0/this->en,this->ki)*pow(this->en_1*this->en_1/(this->en*this->en_2),this->kd) * this->dt_prev;

            //std::cout << "\n New calculated time step = " << dt << std::endl;

            // check time step range
            dt = max(dt,this->dt_min);
            dt = min(dt,this->dt_max);
        }
        

        // check whether simulation has reached the end to avoid exceed maximum simulation time
        if (fabs(time - tmax) > 1.0e-08)
            if ( time+dt>tmax || tmax-(time+dt) < this->dt_min )
                dt = tmax - time;        
           
        // to compute next time-step
        this->dt_prev = dt;
        // storing error measure
        this->storeSolutionChangeinTime();
    }

    std::cout << "\n Next adopted time step = " << dt <<". Var = "<<(dt-this->dt_last)/this->dt_last*100<<"%"<< std::endl;
    
}

void timeStepControlPID::storeSolutionChangeinTime() {
    this->en_2 = this->en_1;
    this->en_1 = this->en;
}

void timeStepControlPID::printSelf (ostream& os, const char* indent) const {
    // Write performance parameters to output file
    os 	<< "\n"
        << indent <<"PID Time-Step Control Performance:"
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

void timeStepControlPID::ckeckKeepMinTimeStep( Real dt, int flow_nonlinear_iteractions, int transport_nonlinear_iteractions) {
    if ( (this->en > 1.0 || flow_nonlinear_iteractions > (this->nsa_max) || transport_nonlinear_iteractions> (this->nsa_max) ) && abs(dt - this->dt_min) < 1.0e-10 )
        this->keep_dt_min = true;
    else
        this->keep_dt_min = false;       
}   