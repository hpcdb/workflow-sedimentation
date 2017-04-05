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
    
//    void setInitDt(double idt);
//    void setDtMin(double dt_min);
//    void setDtMax(double dt_max);
//    void setTolU(double tol_u);
//    void setNumberSucessiveAproxMax(double nsa_max);
//    
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
    
    void computeSolutionChangeInTime(EquationSystems & es);
    
    void checkTimeStepAcceptance(EquationSystems & es, double& dt, unsigned int& t_step, bool& accepted);

    virtual void storeSolutionChangeinTime();
    
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