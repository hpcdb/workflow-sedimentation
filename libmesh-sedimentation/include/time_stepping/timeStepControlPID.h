/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   timeStepControlPID.h
 * Author: arossa
 *
 * Created on 29 de Junho de 2016, 16:31
 */

#ifndef TIMESTEPCONTROLPID_H
#define TIMESTEPCONTROLPID_H

#include <iostream>
#include <sstream>
#include <ostream>
#include <fstream>

#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/transient_system.h"
#include "libmesh/linear_implicit_system.h"
#include "timeStepControlBase.h"

using namespace libMesh;
using namespace std;

class timeStepControlPID: public timeStepControlBase {
    
public:
    timeStepControlPID(double dt_init, double dt_min, double dt_max,  unsigned int nsa_max, double tol_u, double tol_s, double kp, double ki, double kd);
    timeStepControlPID(const timeStepControlPID& orig);
    ~timeStepControlPID();
    
    void computeSolutionChangeInTime(EquationSystems & es);
     
    void checkTimeStepAcceptance(Real dt, int flow_nonlinear_iteractions, int transport_nonlinear_iteractions, bool& accepted);
    
    void computeTimeStep(bool accepted, Real time, Real tmax, Real &dt);

    void storeSolutionChangeinTime();
    
    void ckeckKeepMinTimeStep( Real dt, int flow_nonlinear_iteractions, int transport_nonlinear_iteractions);
    
    double getCurrentError(){
        return this->en;
    };
    
    double getOldError(){
        return this->en_1;
    };
    
    double getOlderError(){
        return this->en_2;
    };    
    
    void setCurrentError( double currentError) {
        this->en = currentError;
    }
    
    void setOldError( double oldError) {
        this->en_1 = oldError;
    }
    
    void setOlderError( double olderError) {
        this->en_2 = olderError;
    }
    
    double getPreviousDt() {
        return this->dt_prev;
    }
    
    void setPreviousDt(double previousDt) {
        this->dt_prev = previousDt;
    }    
    
    void printSelf (ostream& os, const char* indent) const;
    
    friend ostream& operator<<(ostream& os, const timeStepControlPID& o) {
        o.printSelf(os, "\t");
        return os;
    } 
    
private:
    // Flow tolerance between two consecutive solution norm
    double tol_U;
    // Transport tolerance between two consecutive solution norm
    double tol_S;    
    // solution vector change between two consecutive time-step
    double en;
    double en_1;
    double en_2;
    // PID parameters
    double kp, ki, kd;
    // previous TS (updated also when current TS is refused)
    double dt_prev;
};

#endif /* TIMESTEPCONTROLPID_H */