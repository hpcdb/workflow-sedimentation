/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   timeStepControlPC11.h
 * Author: arossa
 *
 * Created on 29 de Junho de 2016, 16:31
 */

#ifndef TIMESTEPCONTROLPC11_H
#define TIMESTEPCONTROLPC11_H

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

class timeStepControlPC11: public timeStepControlBase {
    
public:
    timeStepControlPC11(double dt_init, double dt_min, double dt_max,  unsigned int nsa_max, double tol_u, double tol_s, double theta, double alpha, double k_exp, double s_min, double s_max, bool complete);
    timeStepControlPC11(const timeStepControlPC11& orig);
    ~timeStepControlPC11();
    
    void computeSolutionChangeInTime(EquationSystems & es);
    
    double getError(){
        return this->rn;
    };
    
    //void checkTimeStepAcceptance(EquationSystems & es, double& dt, unsigned int& t_step, bool& accepted);
    
    void storeSolutionChangeinTime();
    
    void checkTimeStepAcceptance(Real dt, int flow_nonlinear_iteractions, int transport_nonlinear_iteractions, bool& accepted);
    
    void computeTimeStep(bool accepted, Real time, Real tmax, Real &dt);
    
    void ckeckKeepMinTimeStep( Real dt, int flow_nonlinear_iteractions, int transport_nonlinear_iteractions);
    
    void printSelf (ostream& os, const char* indent) const;
    
    friend ostream& operator<<(ostream& os, const timeStepControlPC11& o) {
        o.printSelf(os, "\t");
        return os;
    } 
    
private:
    // Flow tolerance between two consecutive solution norm
    double tol_U;
    // Transport tolerance between two consecutive solution norm
    double tol_S;    
    // solution vector change between two consecutive time-step
    double rn, rn_1;
    // exponent to compute new time step
    double k_exp;
    // Safety factor to reduce probability of rejecting the new time step
    double theta;
    // Factor to reduce time-step in case of rejecting the new time step
    double alpha;    
    // Parameters to avoid a strong increase or decrease of subsequent time steps
    double s_min, s_max;
    // flag to compute complete flow solution vector norm (velocity + pressure) or only of the velocity components
    bool complet_flow_norn;
};

#endif /* TIMESTEPCONTROLPC11_H */