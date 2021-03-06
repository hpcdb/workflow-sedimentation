/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   timeStepControlBase.h
 * Author: arossa
 *
 * Created on 29 de Junho de 2016, 16:31
 */

#ifndef TIMESTEPCONTROLBASE_H
#define TIMESTEPCONTROLBASE_H

#include <iostream>
#include <sstream>
#include <ostream>
#include <fstream>

#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"

using namespace libMesh;
using namespace std;

class timeStepControlBase {
public:
    timeStepControlBase(double init_dt, double dt_min, double dt_max, double dt_prev, double dt_avg, unsigned int n_start_control, unsigned int nsa_max, std::string name );
    timeStepControlBase();
    timeStepControlBase(const timeStepControlBase& orig);
    virtual ~timeStepControlBase();
    
//    void setInitDt(double idt);
//    void setDtMin(double dt_min);
//    void setDtMax(double dt_max);
//    void setTolU(double tol_u);
//    void setNumberSucessiveAproxMax(double nsa_max);
    
//    double getInitDt(){
//        return this->init_dt;
//    }
//    double getDtMin(){
//        return this->dt_min;
//    }
//    double getDtMax(){
//        return this->dt_max;
//    }
//    double getTolU(){
//        return this->tol_U;
//    }
//    double getTolS(){
//        return this->tol_S;
//    }    
//    double getNumberSucessiveAproxMax(){
//        return this->nsa_max;
//    }
    
    virtual void computeSolutionChangeInTime(EquationSystems & es) = 0;
    
    virtual double getCurrentError() = 0;
    
    virtual double getOldError() = 0;
    
    virtual double getOlderError() = 0;
    
    virtual void setCurrentError( double currentError) = 0;
    
    virtual void setOldError(double oldError) = 0;
    
    virtual void setOlderError(double olderError) = 0;
    
    virtual void checkTimeStepAcceptance(Real dt, int flow_nonlinear_iteractions, int transport_nonlinear_iteractions, bool& accepted) = 0;
    
    virtual void computeTimeStep(bool accepted, Real time, Real tmax, Real &dt) = 0;
    
    virtual void storeSolutionChangeinTime() = 0;
    
    virtual void ckeckKeepMinTimeStep( Real dt, int flow_nonlinear_iteractions, int transport_nonlinear_iteractions) = 0;
    
    virtual void printSelf (ostream& os, const char* indent) const = 0;
    
    virtual double getPreviousDt() = 0;
    
    virtual void setPreviousDt(double previousDt) = 0;
    
    unsigned int getStartTimeStepControl() const {
        return this->start_control;
    }
    
    double getLastAcceptedTS() const {
        return this->dt_last;
    }
    
    void setLastAcceptedTS( double lastDt) {
        this->dt_last = lastDt;
    }
    
    double getAverageTS() const {
        return this->dt_avg;
    }
    
    void setAverageTS( double avgDt) {
        this->dt_avg = avgDt;
    }
    
    void updateAverageTS( double avgDt) {
        this->dt_avg += avgDt;
    }
    
    unsigned int getNumberAcceptedTS() {
        return this->n_accepted_time_step;
    }
    
    void setNumberAcceptedTS (unsigned int n_accepted_ts) {
        this->n_accepted_time_step = n_accepted_ts;
    }
    
    unsigned int getNumberRejectedTS() {
        return this->n_rejected_time_step;
    }
    
    void setNumberRejectedTS (unsigned int n_rejected_ts) {
        this->n_rejected_time_step = n_rejected_ts;
    }    

    friend ostream& operator<<(ostream& os, const timeStepControlBase& o) {
        o.printSelf(os, "\t");
        return os;
    } 
    
protected:
    // model name
    std::string model_name;
    // initial time-step
    double init_dt;
    // minimum time-step
    double dt_min;
    // maximum time-step
    double dt_max;
    // old time step (last accepted time-step size)
    double dt_last;
    // average time-step
    double dt_avg;    
    // time-step control statistic parameters
    unsigned int n_rejected_nonlinear_flow_iterations_total;
    unsigned int n_rejected_nonlinear_transport_iterations_total;
    unsigned int n_rejected_linear_iterations_total;
    unsigned int n_rejected_linear_iterations_per_ts;
    // updating number of rejected time-steps
    unsigned int n_rejected_time_step;    
    // updating number of accepted time-steps
    unsigned int n_accepted_time_step;
    // number of time-steps before start to control 
    unsigned int start_control;
    // maximum number of non linear iterations (successive approximation)
    unsigned int nsa_max;
    // flag to check whether a time-step solution was accepted only because minimum time-step value is already reached
    bool keep_dt_min;
};

#endif /* TIMESTEPCONTROLBASE_H */