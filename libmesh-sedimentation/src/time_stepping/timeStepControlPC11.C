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
#include <iomanip>      // std::setprecision

#include "timeStepControlPC11.h"

using namespace libMesh;
using namespace std;

timeStepControlPC11::timeStepControlPC11( double dt_init, double dt_min, double dt_max, unsigned int nsa_max, double tol_U, double tol_S, double theta, double alpha, double k_exp, double s_min, double s_max, bool complete):
                                        timeStepControlBase(dt_init, dt_min, dt_max, dt_init, dt_init*2, 2, nsa_max, "PC11"),
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
//
//void timeStepControl::setInitDt(double idt) {
//    this->init_dt = idt;
//}

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
    
        cout<<"Complete Solution L2 norm = "<<flow_norm<<" Complete Diff Solution L2 norm = "<<flow_diff_norm<<endl;
    
        // compute the measure of the change in time for flow
        if (flow_diff_norm>1.0e-10) { // avoid 'e_flow' becomes zero!
            r_flow = flow_diff_norm / flow_norm;
            // divide by tolerance here to compare this value with that floe transport problem
            r_flow /= this->tol_U;
        }
        cout<<" rflow_complete = "<< r_flow<<endl;
    } else { // to compute the only velocity solution vector L2 norm (excludes pressure)
    
        // Get a conSTant reference to the mesh object.
        const MeshBase& mesh = es.get_mesh();

        // The dimension that we are running
        const unsigned int dim = mesh.mesh_dimension();

        MPI_Comm comm = MPI_COMM_WORLD;

        double global_norm, global_dif_norm, local_sqr_norm = 0.0;
        double local_sqr_dif_norm = 0.0;
        int ndata = 2;
        double rms[2], rms_sum[2];

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

        //cout<<"Velocity L2 norm = "<<global_norm<<" Diff velocity L2 norm = "<<global_dif_norm<<endl;

        // compute the measure of the change in time for flow
        if (global_dif_norm>1.0e-10) { // avoid 'r_flow' become zero!
            r_flow = global_dif_norm / global_norm;
            r_flow /= this->tol_U;
            //cout<<" r_flow velocity = "<< r_flow<<endl;
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
        if (transport_diff_norm>1.0e-10) { // avoid 'e_transport' becomes zero!
            r_transport = transport_diff_norm / transport_norm;
            // divide by tolerance here to compare this value and that from flow problem
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
    if ( (this->rn>1.0  || flow_nonlinear_iteractions > this->nsa_max || transport_nonlinear_iteractions > this->nsa_max ) && dt > this->dt_min)
        accepted = false;
    
    // check whether the current time-step was accepted because it was already the minimum 
    if (accepted) { 
        this->ckeckKeepMinTimeStep(dt, flow_nonlinear_iteractions, transport_nonlinear_iteractions);
    }      
    
}

void timeStepControlPC11::computeTimeStep(bool accepted, Real time, Real tmax, Real& dt)
{
    
    if(!accepted)
    {
        cout << "\n Current solution was REJECTED";
        
        this->dt_avg -= dt;
       
        // For while next time-step will be reduced by an alpha factor [0,1]
        dt = this->alpha * dt;  
        
        //std::cout << " New calculated time step = " << dt << std::endl;
        
        // For while next time-step will be reduced by an alpha factor [0,1] but not smaller than dt_min!
        dt = max(dt, this->dt_min);        
            
        // updating number of rejected time-steps
        this->dt_avg += dt;
        this->n_rejected_time_step++;
        
    }else
    {   
        std::cout << "\n Current solution was ACCEPTED";
                
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

        // updating number of accepted time-steps
        // to compute average time-step
        this->dt_avg+=dt;
        this->n_accepted_time_step++;
        // storing error measure
        this->storeSolutionChangeinTime();     
            
    }

    // check whether simulation has reached the end to avoid exceed maximum simulation time
    if (time+dt>tmax || tmax-(time+dt) < this->dt_min)
        dt = tmax - time;

    std::cout << "\n Next adopted time step = " << dt <<". Var = "<<(dt-this->dt_last)/this->dt_last*100 <<"%"<< std::endl;
}



/*
void timeStepControlPC11::checkTimeStepAcceptance(EquationSystems & es, double& dt, bool& accepted) {
    
    
    unsigned int flow_nli_counter      = es.parameters.get<unsigned int>("n_non_linear_iter_flow");
    unsigned int transport_nli_counter = es.parameters.get<unsigned int>("n_non_linear_iter_transport"); 
    
   // bool min_dt = false;
    double dt_star;// = this->dt_min;
    
    if ( (this->rn>1.0 || flow_nli_counter>(this->nsa_max-1) || transport_nli_counter>(this->nsa_max-1) ) && dt > this->dt_min) { // current time-step is refused  
        
        std::cout << " Current solution was refused!\n Reverting time data" << std::endl;
        
        // Get a reference to the Flow and Transport systems object and restore the solution vectors from the previous time step
        TransientLinearImplicitSystem & flow_system_reference = es.get_system<TransientLinearImplicitSystem> ("flow");
        *flow_system_reference.current_local_solution = *flow_system_reference.old_local_solution;    
        if ( es.has_system("Transport") ) {        
            TransientLinearImplicitSystem & transport_system_reference =  es.get_system<TransientLinearImplicitSystem> ("Transport");
            *transport_system_reference.current_local_solution = *transport_system_reference.old_local_solution;
        }

        es.parameters.set<Real> ("time") -= dt;
        this->dt_avg -= dt;
        t_step--;
        
        // For while next time-step will be reduced by an alpha factor [0,1] but not smaller than dt_min!
        dt = max(this->alpha * dt, this->dt_min);
        std::cout << " New calculated time step = " << dt << std::endl;
            
        // updating number of rejected time-steps
        this->n_rejected_time_step++;
        this->dt_avg += dt;
        accepted = false;

    } else { //current time-step is accepted

        std::cout << " Current solution was accepted";
             
        if (fabs((dt-this->dt_min)<1.0e-10) && (this->rn>1.0 || flow_nli_counter>(this->nsa_max-1) || transport_nli_counter>(this->nsa_max-1) ) ) {
            std::cout << " because dt=dt_min!"<<std::endl;
            //storing last accepted time-step
            this->dt_last = dt;
        } else {
            std::cout << "!"<<endl;
            // computing the new new time step whatever the current time step was accepted or not
            double exp = 1.0/(this->k_exp+1.0);
//            cout<<setprecision(9)<<" k_exp = "<<this->k_exp<<endl;
//            cout<<" teta = "<<this->theta<<endl;
//            cout<<" rn_1 = "<<this->rn_1<<endl;
//            cout<<" rn = "<<this->rn<<endl;
//            cout<<" dt = "<<dt<<endl;
//            cout<<" dt_last = "<<this->dt_last<<endl;
//            cout<<" s_max = "<<this->s_max<<endl;
//            cout<<" s_min = "<<this->s_min<<endl;

            dt_star = this->theta * pow( this->rn_1/(this->rn*this->rn ),exp) * dt*dt/this->dt_last;
            cout<<" dt_star = "<<dt_star<<endl;

            //storing last accepted time-step
            this->dt_last = dt;

            dt = std::min(this->s_max*dt, std::max(this->s_min*dt, dt_star) );

            std::cout << " New calculated time step = " << dt <<std::endl;

            // check time step range
            dt = max(dt,this->dt_min);
            dt = min(dt,this->dt_max);

            // check whether simulation has reached the end to avoid exceed maximum simulation time
            double time = es.parameters.get<Real>("time");
            double tmax = es.parameters.get<Real>("tmax");
            if (time+dt>tmax || tmax-(time+dt) < this->dt_min)
                dt = tmax - time;
            
            // to compute average time-step
            this->dt_avg+=dt;
            
        }  // end "else" time-step accepted
        
        std::cout << " New adopted time step = " << dt <<". Variation = "<<(dt-this->dt_last)/this->dt_last*100 <<"%"<< std::endl;
        // updating number of accepted time-steps
        this->n_accepted_time_step++;
        this->rn_1 = this->rn;        
        accepted = true;     
    }
    
    // computing the new new time step whatever the current time step was accepted or not
//    double exp = 1.0/(this->k_exp+1.0);
//    cout<<setprecision(9)<<" k_exp = "<<this->k_exp<<endl;
//    cout<<" teta = "<<this->teta<<endl;
//    cout<<" rn_1 = "<<this->rn_1<<endl;
//    cout<<" rn = "<<this->rn<<endl;
//    cout<<" dt_last = "<<this->dt_last<<endl;
//    cout<<" s_max = "<<this->s_max<<endl;
//    cout<<" s_min = "<<this->s_min<<endl;
//    if (!min_dt) {
//        dt_star = this->teta * pow( this->rn_1/(this->rn*this->rn ),exp) * dt*dt/this->dt_last;
//        cout<<" dt_star = "<<dt_star<<endl; 
//    }
//    // updating parameters
//    if (accepted) {
//        this->dt_last = dt;
//        this->rn_1 = this->rn;
//    }  
//    
//    cout<<" dt_last = "<<this->dt_last<<endl;
//    
//    dt = std::min(this->s_max*dt, std::max(this->s_min*dt, dt_star) );
//    
//    std::cout << " New calculated time step = " << dt << std::endl;
//
//    // check time step range
//    dt = max(dt,this->dt_min);
//    dt = min(dt,this->dt_max);
//                
//    // check whether simulation has reached the end to avoid exceed maximum simulation time
//    double time = es.parameters.get<Real>("time");
//    double tmax = es.parameters.get<Real>("tmax");
//    if (time+dt>tmax || tmax-(time+dt) < 1.0e-03)
//        dt = tmax - time; 
                    
//    std::cout << " New adopted time step = " << dt << std::endl;       
} 
 */

/*
void timeStepControlPC11::checkTimeStepAcceptance(EquationSystems & es, double& dt, unsigned int& t_step, bool& accepted) {
    
    
    unsigned int flow_nli_counter = es.parameters.get<unsigned int>("n_non_linear_iter_flow");
    unsigned int transport_nli_counter = es.parameters.get<unsigned int>("n_non_linear_iter_transport"); 
    
   // bool min_dt = false;
    double dt_star, dt_temp = this->dt_min;
    
   // if (dt>this->dt_min ) {
    
        // compute new time time-step whatever the solution measure
        double exp = 1.0/(this->k_exp+1.0);
        cout<<setprecision(9)<<" k_exp = "<<this->k_exp<<endl;
        cout<<" teta = "<<this->theta<<endl;
        cout<<" rn_1 = "<<this->rn_1<<endl;
        cout<<" rn = "<<this->rn<<endl;
        cout<<" dt = "<<dt<<endl;
        cout<<" dt_last = "<<this->dt_last<<endl;
        cout<<" s_max = "<<this->s_max<<endl;
        cout<<" s_min = "<<this->s_min<<endl;
        dt_star = this->theta * pow( this->rn_1/(this->rn*this->rn ),exp) * dt*dt/this->dt_last;
        cout<<" dt_star = "<<dt_star<<endl; 

        dt_temp = std::min(this->s_max*dt, std::max(this->s_min*dt, dt_star) );

        std::cout << " New calculated time step = " << dt_temp << std::endl;

        // check time step range
        dt_temp = max(dt_temp,this->dt_min);
        dt_temp = min(dt_temp,this->dt_max);

        // check whether simulation has reached the end to avoid exceed maximum simulation time
        double time = es.parameters.get<Real>("time");
        double tmax = es.parameters.get<Real>("tmax");
        if (time+dt>tmax || tmax-(time+dt_temp) < this->dt_min)
            dt_temp = tmax - time;
    //}
    
    std::cout << " New time step after range check= " << dt_temp << std::endl;
    
    if ( (this->rn>1.0 || flow_nli_counter>(this->nsa_max-1) || transport_nli_counter>(this->nsa_max-1) ) && dt > this->dt_min) { // current time-step is refused  
        
        std::cout << " Current solution was refused!\n Reverting time data" << std::endl;
        
        // Get a reference to the Flow and Transport systems object and restore the solution vectors from the previous time step
        TransientLinearImplicitSystem & flow_system_reference = es.get_system<TransientLinearImplicitSystem> ("Flow");
        *flow_system_reference.current_local_solution = *flow_system_reference.old_local_solution;    
        if ( es.has_system("Transport") ) {        
            TransientLinearImplicitSystem & transport_system_reference =  es.get_system<TransientLinearImplicitSystem> ("Transport");
            *transport_system_reference.current_local_solution = *transport_system_reference.old_local_solution;
        }

        es.parameters.set<Real> ("time") -= dt;
        this->dt_avg -= dt;
        t_step--;
        
        // For while next time-step will be reduced by theta factor but not smaller than dt_min!
        if (dt_temp>dt)
            dt_temp = max(this->theta * dt, this->dt_min);
        std::cout << " New adopted time step = " << dt_temp << std::endl;
            
        // updating number of rejected time-steps
        this->n_rejected_time_step++;
        this->dt_avg += dt;
        accepted = false;

    } else { //current time-step is accepted
        std::cout << " Current solution was accepted";
        if (fabs((dt-this->dt_min)<1.0e-10) && (this->rn>1.0 || flow_nli_counter>(this->nsa_max-1) || transport_nli_counter>(this->nsa_max-1) ) ) {
            std::cout << " because dt=dt_min!"<<std::endl;
        } else {
            std::cout << "!"<<endl;
        }  // end "else" time-step accepted
        
        std::cout << " New adopted time step = " << dt_temp << std::endl;
        // updating number of accepted time-steps
        this->n_accepted_time_step++;
        this->rn_1 = this->rn;
        this->dt_last = dt;
        // to compute average time-step
        this->dt_avg+=dt;
        accepted = true;     
    }
    dt = dt_temp;   
}
 */

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
        << "\n\n Effective time-step size average: "
        << this->dt_avg/this->n_accepted_time_step            
		<< std::endl;    
}

void timeStepControlPC11::ckeckKeepMinTimeStep( Real dt, int flow_nonlinear_iteractions, int transport_nonlinear_iteractions) 
{
    if ( (this->rn > 1.0 || flow_nonlinear_iteractions > (this->nsa_max) || transport_nonlinear_iteractions > (this->nsa_max) ) && abs(dt - this->dt_min) < 1.0e-10 )
        this->keep_dt_min = true;
    else
        this->keep_dt_min = false;       
}
    