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
                                        timeStepControlBase(dt_init, dt_min, dt_max, dt_init, 3*dt_init, 3, nsa_max, "PID"),
                                        tol_U(tol_U),
                                        tol_S(tol_S),
                                        kp(kp),
                                        ki(ki),
                                        kd(kd),
                                        en_1(1.0),
                                        en_2(1.0),
                                        dt_prev(dt_init)
{
    cout<<"PID time-step control will be adopted\n\n"; 
}

timeStepControlPID::timeStepControlPID(const timeStepControlPID& orig):
                                       timeStepControlBase() 
{
    cout<<"contrutor de copia\n";
}

timeStepControlPID::~timeStepControlPID() {
}
//
//void timeStepControl::setInitDt(double idt) {
//    this->init_dt = idt;
//}

void  timeStepControlPID::computeSolutionChangeInTime(EquationSystems & es) {
    
    // solution vector and its variation in respect to last NLS L2 norms for the current time step
    Real e_flow = 1.0e-10, flow_norm, flow_diff_norm;
    Real e_transport = 1.0e-10;
    
    // Get a reference to the Flow system object.
    // PS: There is always a "flow" system 
    TransientLinearImplicitSystem & flow_system_reference =
    es.get_system<TransientLinearImplicitSystem> ("flow");
    
    // Compute the l2 norm of the flow solution (includes pressure)
    //flow_norm = flow_system_reference.solution->l2_norm();
    
    // Here I get a copy of the solution vector at the end of the nonlinear solution process.
    // This vector will be used to measure the change of the quantities of interest in time, i.e.,
    // the variation on velocity vector and pressure over two successive solutions over time.
    UniquePtr<NumericVector<Number> > current_flow_soln (flow_system_reference.current_local_solution->clone()); // "current_local_solution" was already updated after solve

    // Compute the difference between the current and last (in time) solutions.
    current_flow_soln->add (-1., *flow_system_reference.old_local_solution);
         
    // Close the vector before computing its norm
    current_flow_soln->close();

    // Compute the l2 norm of the difference between two successive solutions
    //flow_diff_norm = current_flow_soln->l2_norm();
    
    //cout<<"norma LM = "<<flow_norm<<" norm_delta = "<<flow_diff_norm<<endl;
    
    // compute the measure of the change in time for flow
    //if (flow_diff_norm>1.0e-10) { // avoid 'e_flow' become zero!
    //    e_flow = flow_diff_norm / flow_norm;
    //    e_flow /= this->tol_U;
    //}
    
   // cout<<" eflow_before = "<< e_flow<<endl;
    
    // Get a conSTant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

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
        local_sqr_norm += (flow_system_reference.solution->el(i)*flow_system_reference.solution->el(i) + flow_system_reference.solution->el(i+1)*flow_system_reference.solution->el(i+1));
        local_sqr_dif_norm += (current_flow_soln->el(i)*current_flow_soln->el(i) + current_flow_soln->el(i+1)*current_flow_soln->el(i+1));
    }

    rms[0] = local_sqr_norm;
    rms[1] = local_sqr_dif_norm;

    MPI_Allreduce(rms, rms_sum, ndata, MPI_DOUBLE, MPI_SUM, comm);
     
    global_norm = std::pow(rms_sum[0],0.5);
    global_dif_norm = std::pow(rms_sum[1],0.5);
     
    //cout<<" minha norma = "<<global_norm<<" norm_delta = "<<global_dif_norm<<endl;
     

    // compute the measure of the change in time for flow
    if (global_dif_norm>1.0e-10) { // avoid 'e_flow' become zero!
        e_flow = global_dif_norm / global_norm;
        e_flow /= this->tol_U;
        //cout<<" eflow_after = "<< e_flow<<endl;
    }
     
    
    // Get a reference to the Transport system object.
    // PS: Must be checked cause "transport" system may not be defined
    if ( es.has_system("transport") ) {
        
        std::cout << " \nChange Flow solution in time = " << e_flow;
        
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
        if (transport_diff_norm>1.0e-10) { // avoid 'e_transport' become zero!
            e_transport = transport_diff_norm / transport_norm;
            e_transport /= this->tol_S;
        }
        
        std::cout << " \nChange Transport solution in time = " << e_transport << std::endl;        
    }
            
    // compute the maximum measure of the change in time for the systems
    this->en = std::max(e_flow, e_transport);

    std::cout << " \nMeasure of the solution change in time = " << this->en << std::endl;
}

void timeStepControlPID::checkTimeStepAcceptance(EquationSystems & es, double& dt, unsigned int& t_step, bool& accepted) {
    
    unsigned int flow_nli_counter = es.parameters.get<unsigned int>("n_non_linear_iter_flow");
    unsigned int transport_nli_counter = es.parameters.get<unsigned int>("n_non_linear_iter_transport");
    
    if ( (this->en>1.0 || flow_nli_counter>(this->nsa_max) || transport_nli_counter>(this->nsa_max) ) && dt > this->dt_min)     { // current time-step is refused  
        
        std::cout << " Current solution was refused!\n Reverting time data" << std::endl;
        
        // Get a reference to the Flow and Transport systems object and restore the solution vectors from the previous time step
        TransientLinearImplicitSystem & flow_system_reference = es.get_system<TransientLinearImplicitSystem> ("flow");
        *flow_system_reference.current_local_solution = *flow_system_reference.old_local_solution;    
        if ( es.has_system("Transport") ) {        
            TransientLinearImplicitSystem & transport_system_reference =  es.get_system<TransientLinearImplicitSystem> ("transport");
            *transport_system_reference.current_local_solution = *transport_system_reference.old_local_solution;
        }

        es.parameters.set<Real> ("time") -= dt;
        this->dt_avg -= dt;
        t_step--;

        Real factor = 1.0/this->en;

        if (factor>0.8)
            factor = 0.8;
            std::cout << " dt*factor = "<< factor*dt<<endl;
                
            dt = max(factor*dt,this->dt_min);
            std::cout << " New calculated time step = " << dt << std::endl;
            this->dt_avg += dt;

            this->dt_prev = dt*dt / this->dt_prev;
                
//            // updating number of linear iterations rejected
//            this->n_rejected_linear_iterations_total+=this->n_rejected_linear_iterations_per_ts;
//            // updating number of non-linear iterations rejected
//            this->n_rejected_nonlinear_flow_iterations_total+=flow_nli_counter;
//            this->n_rejected_nonlinear_transport_iterations_total+=transport_nli_counter;                     
            // updating number of rejected time-steps
            this->n_rejected_time_step++;
            accepted = false;

        } else { //current time-step is accepted

            //storing last accepted time-step
            this->dt_last = dt;

            std::cout << " Current solution was accepted";
                
            if (fabs((dt-this->dt_min)<1.0e-10) && (this->en>1.0 || flow_nli_counter>(this->nsa_max) || transport_nli_counter>(this->nsa_max) ) )
                std::cout << " because dt=dt_min!"<<std::endl;
            else {
                std::cout << "!"<<endl;

                dt = pow(this->en_1/this->en,this->kp)*pow(1.0/this->en,this->ki)*pow(this->en_1*this->en_1/(this->en*this->en_2),this->kd) * this->dt_prev;

                std::cout << " New calculated time step = " << dt << std::endl;

                // check time step range
                dt = max(dt,this->dt_min);
                dt = min(dt,this->dt_max);
                
                // check whether simulation has reached the end to avoid exceed maximum simulation time
                double time = es.parameters.get<Real>("time");
                double tmax = es.parameters.get<Real>("tmax");
                if (time+dt>tmax || tmax-(time+dt) < this->dt_min)
                    dt = tmax - time; 
                    
                std::cout << " New adopted time step = " << dt <<". Variation = "<<(dt-this->dt_last)/this->dt_last*100<<"%"<< std::endl;
            }

            // updating parameters
            this->dt_avg += dt;
            this->dt_prev = dt;
            this->en_2 = this->en_1;
            this->en_1 = this->en;
            accepted = true;
            //cout<<" # accepted TS = "<<this->n_accepted_time_step<<endl;
            // updating number of accepted time-steps
            this->n_accepted_time_step++;   
        }  // end "else" time-step accepted
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
//        << "\n\n Sum effective time-steps: "
//        << this->dt_m            
        << "\n\n Effective time-step size average: "
        << this->dt_avg/this->n_accepted_time_step            
		<< std::endl;    
}
    