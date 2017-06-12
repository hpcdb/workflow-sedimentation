/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   timeStepControlResidual.cpp
 * Author: arossa
 * 
 * Created on 29 de Junho de 2016, 16:31
 */

#if 0 

#include <iostream>
#include <sstream>
#include <ostream>
#include <fstream>

#include "timeStepControlResidual.h"

using namespace libMesh;
using namespace std;

timeStepControlResidual::timeStepControlResidual( double dt_init, double dt_min, double dt_max, unsigned int n_target_flow, unsigned int n_target_transp, unsigned int max_nsa, unsigned int n_limit_flow, unsigned int n_limit_transp,
                                                double factor_max, double factor_min, double factor_red):
                                                timeStepControlBase(dt_init, dt_min, dt_max, dt_init, dt_init, 1, max_nsa, "Residual"),
                                                nsa_target_flow(n_target_flow),
                                                nsa_target_transp(n_target_transp),
                                                nsa_limit_flow(n_limit_flow),
                                                nsa_limit_transp(n_limit_transp),
                                                dt_factor_max(factor_max),
                                                dt_factor_min(factor_min),
                                                dt_factor_red(factor_red)
{
    cout<<"Residual time-step control will be adopted\n\n"; 
}

timeStepControlResidual::timeStepControlResidual(const timeStepControlResidual& orig):
                                                 timeStepControlBase() 
{
    cout<<"contrutor de copia\n";
}

timeStepControlResidual::~timeStepControlResidual() {
}


void  timeStepControlResidual::computeSolutionChangeInTime(EquationSystems & es) {
    
    unsigned int current_n_flow_iter = es.parameters.get<unsigned int> ("n_non_linear_iter_flow");
    unsigned int old_n_flow_iter = es.parameters.get<unsigned int> ("old_n_non_linear_iter_flow");
    unsigned int current_n_trans_iter = es.parameters.get<unsigned int> ("n_non_linear_iter_transport");
    unsigned int old_n_trans_iter = es.parameters.get<unsigned int> ("old_n_non_linear_iter_transport");
    
    double slope, b;
    char str0[256], str1[256]="\n";
       
    double m_factor_flow, m_factor_transport = this->dt_factor_max;
    
    // First range:
    // -----------
    // Flow problem    
    if (current_n_flow_iter <= this->nsa_target_flow) { // current number of successive approximations is smaller than "target" for flow.
        slope = (1.0 - this->dt_factor_max)/(this->nsa_target_flow-1);
        b = 1.0 - slope * this->nsa_target_flow;
        m_factor_flow = slope * current_n_flow_iter + b;
        sprintf(str0,"\nTarget reached for Flow Problem! Current # of NL iterations (%d) is smaller/equal than # NLI target (%d). \n", current_n_flow_iter, this->nsa_target_flow);
        sprintf(str1,"TS is gonna be multiplied by %f !\n", m_factor_flow);
    // Second range:
    } else if (current_n_flow_iter > this->nsa_target_flow && current_n_flow_iter <= this->nsa_limit_flow) { // current number of NLI iterations is greater than target but smaller than limit
        slope = (this->dt_factor_min - 1.0)/(this->nsa_limit_flow-this->nsa_target_flow);
        b = 1.0 - slope * this->nsa_target_flow;
        m_factor_flow = slope * current_n_flow_iter + b;
        sprintf(str0,"\nFlow: Current number of NL iterations (%d) approached the limit (%d).", current_n_flow_iter, this->nsa_limit_flow);
        //sprintf(str1,"\n iternew > iternew_target (%d) && <= iternew_lim (%d)\n", iternew_target, iternew_lim);
        sprintf(str1,"\nTS is gonna be reduced by %f %%!\n", (1.0-m_factor_flow)*100);
    // Limit exceeded
    } else { // current number of NLI iterations exceeded the limit: time step will be reduced by a user input factor
        m_factor_flow = this->dt_factor_red; // to reduct TS size. May be linearized too (according the # nsa)!!!
        sprintf(str0,"\nFlow: Current number of NL iterations (%d) exceeded the limit (%d).", current_n_flow_iter, this->nsa_limit_flow);
        sprintf(str1,"\nTS is gonna be reduced by %f %%!\n", (1.0-m_factor_flow)*100);
    }
    cout<<str0<<"\n"<<str1<<endl;
    
    // Get a reference to the Transport system object.
    // PS: Must be checked cause "transport" system may not be defined
    if ( es.has_system("transport") ) {
        char str0[256], str1[256]="\n";
        // Transport problem 
        if (current_n_trans_iter <= this->nsa_target_transp) { // current number of successive approximations is smaller than "target" for flow.
            slope = (1.0 - this->dt_factor_max)/(this->nsa_target_flow-1);
            b = 1.0 - slope * this->nsa_target_transp;
            m_factor_transport = slope * current_n_trans_iter + b;
            sprintf(str0,"\nTarget reached for Transport Problem! Current # of NL iterations (%d) is smaller/equal than # NLI target (%d).\n", current_n_trans_iter, this->nsa_target_transp);
            sprintf(str1,"TS is gonna be multiplied by %f !\n", m_factor_transport);
        // Second range:
        } else if (current_n_trans_iter > this->nsa_target_transp && current_n_flow_iter <= this->nsa_limit_transp) { // current number of NLI iterations is greater than target but smaller than limit
            slope = (this->dt_factor_min - 1.0)/(this->nsa_limit_transp-this->nsa_target_transp);
            b = 1.0 - slope * this->nsa_target_transp;
            m_factor_transport = slope * current_n_trans_iter + b;
            sprintf(str0,"\nTransport: Current number of NL iterations (%d) approached the limit (%d).", current_n_trans_iter, this->nsa_limit_transp);
            //sprintf(str1,"\n iternew > iternew_target (%d) && <= iternew_lim (%d)\n", iternew_target, iternew_lim);
            sprintf(str1,"\nTS is gonna be reduced by %f %%!\n", (1.0-m_factor_transport)*100);
        // Limit exceeded
        } else { // current number of NLI iterations exceeded the limit: time step will be reduced by a user input factor
            m_factor_transport= this->dt_factor_red; // to reduct TS size. May be linearized too (according the # nsa)!!!
            sprintf(str0,"\nTransport: Current number of NL iterations (%d) exceeded the limit (%d)", current_n_trans_iter, this->nsa_limit_transp);
            sprintf(str1,"\nTS is gonna be reduced by %f %%!\n", (1.0-m_factor_transport)*100);
        }
        cout<<str0<<"\n"<<str1<<endl;
    }
    // Actual multiplication factor will be the smaller between flow and transport problems
    this->mult_factor = min(m_factor_flow,m_factor_transport);
    
    // But whether current NNLI has grown and it became greater than the target, the next time step will be reduced according current and last NNLI ratio
    double n_flow_iter_ratio =  2.0;
    if (current_n_flow_iter>0)
        n_flow_iter_ratio = double(old_n_flow_iter)/current_n_flow_iter;
    double n_trans_iter_ratio =  2.0;
    if ( es.has_system("transport") )
        if (current_n_trans_iter>0)
            n_trans_iter_ratio = double(old_n_trans_iter)/current_n_trans_iter;
    
    if ( (n_flow_iter_ratio < 1.0 || n_trans_iter_ratio < 1.0) && (current_n_flow_iter > this->nsa_target_flow || current_n_trans_iter > this->nsa_target_transp) ) { // if (current # NLI) > (last # NLI), ensure to use a smaller TS than the last one
        sprintf(str0,"\n# Non-Linear iterations has grown: current number (%d) > last number (%d)", current_n_flow_iter, old_n_flow_iter );
        this->mult_factor = min(n_flow_iter_ratio,n_trans_iter_ratio);
        sprintf(str1,"\nTS is gonna be reduced by %f !\n", this->mult_factor);
        cout<<str0<<"\n"<<str1<<endl;        
    }
}

void timeStepControlResidual::checkTimeStepAcceptance(EquationSystems & es, double& dt, unsigned int& t_step, bool& accepted) {
    
    unsigned int flow_nli_counter = es.parameters.get<unsigned int>("n_non_linear_iter_flow");
    unsigned int transport_nli_counter = es.parameters.get<unsigned int>("n_non_linear_iter_transport");
        
    if ( flow_nli_counter>(this->nsa_max) || transport_nli_counter>(this->nsa_max)) { // current time-step is refused  
        
        std::cout << " Current solution was refused!\n Reverting time data" << std::endl;
        
        // Get a reference to the Flow and Transport systems object and restore the solution vectors from the previous time step
        TransientLinearImplicitSystem & flow_system_reference = es.get_system<TransientLinearImplicitSystem> ("flow");
        *flow_system_reference.current_local_solution = *flow_system_reference.old_local_solution;    
        if ( es.has_system("transport") ) {        
            TransientLinearImplicitSystem & transport_system_reference =  es.get_system<TransientLinearImplicitSystem> ("transport");
            *transport_system_reference.current_local_solution = *transport_system_reference.old_local_solution;
        }

        es.parameters.set<Real> ("time") -= dt;
        this->dt_avg -= dt;
        t_step--;
                
        // updating number of linear iterations rejected
        this->n_rejected_linear_iterations_total+=this->n_rejected_linear_iterations_per_ts;
        // updating number of non-linear iterations rejected
        this->n_rejected_nonlinear_flow_iterations_total+=flow_nli_counter;
        this->n_rejected_nonlinear_transport_iterations_total+=transport_nli_counter;
        // updating number of rejected time-steps
         this->n_rejected_time_step++;
        accepted = false;
                
    } else { //current time-step is accepted

        std::cout << " Current solution was accepted!\n";
        //cout<<" # accepted TS = "<<this->n_accepted_time_step<<endl;
        // updating number of accepted time-steps
        this->n_accepted_time_step++;
        accepted = true;
    }  // end "else" time-step accepted
    
    // compute new time-step
    dt = this->mult_factor*dt;
    std::cout << " New calculated time step = " << dt << std::endl;

    // check time step range
    dt = max(dt,this->dt_min);
    dt = min(dt,this->dt_max);
    
    // check whether simulation has reached the end to avoid exceed maximum simulation time
    double time = es.parameters.get<Real>("time");
    double tmax = es.parameters.get<Real>("tmax");
    if (time+dt>tmax || tmax-(time+dt) < 1.0e-03)
        dt = tmax - time;
    //cout<<" dt_avg = "<<this->dt_avg<<endl;
    // to compute average time-step
    this->dt_avg+=dt;    
    std::cout << " New adopted time step = " << dt << std::endl;
}

void timeStepControlResidual::printSelf (ostream& os, const char* indent) const {
    // Write performance parameters to output file
    os 	<< "\n"
        << indent <<"Residual Time-Step Control Performance:"
        << "\nTotal number of time-steps: "
        << this->n_accepted_time_step + this->n_rejected_time_step
        << "\nNumber of accepted time-steps: "
        << this->n_accepted_time_step           
        << "\nNumber of rejected time-steps: "
        << this->n_rejected_time_step            
        << "\nRate Rejected/Total time-steps: "
        << 100.0 * this->n_rejected_time_step/(this->n_accepted_time_step + this->n_rejected_time_step)<<"%"
        << "\n\n Effective time-step size average: "
        << this->dt_avg/this->n_accepted_time_step            
		<< std::endl;    
}
  
#endif