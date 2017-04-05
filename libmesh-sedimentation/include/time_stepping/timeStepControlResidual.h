/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   timeStepControlResidual.h
 * Author: arossa
 *
 * Created on 30 de Junho de 2016, 16:31
 */

#ifndef TIMESTEPCONTROLRESIDUAL_H
#define TIMESTEPCONTROLRESIDUAL_H

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

class timeStepControlResidual: public timeStepControlBase {
    
public:
    timeStepControlResidual(double dt_init, double dt_min, double dt_max, unsigned int n_target_flow, unsigned int n_target_transport, 
                            unsigned int n_limit_flow, unsigned int n_limit_transport, unsigned int nsa_max, double factor_max, double factor_min, double factor_red);
    timeStepControlResidual(const timeStepControlResidual& orig);
    ~timeStepControlResidual();
    
    void computeSolutionChangeInTime(EquationSystems & es);
    
    void checkTimeStepAcceptance(EquationSystems & es, double& dt, unsigned int& t_step, bool& accepted);
    
    void storeSolutionChangeinTime (){ };
    
    void printSelf (ostream& os, const char* indent) const;
    
    friend ostream& operator<<(ostream& os, const timeStepControlResidual& o) {
        o.printSelf(os, "\t");
        return os;
    } 
    
private:
    // number of non-linear iteration target for flow problem
    unsigned int nsa_target_flow;
    // number of non-linear iteration target for transport problem
    unsigned int nsa_target_transp;    
    // number of non-linear iterations limit before time-step reduction for flow problem
    unsigned int nsa_limit_flow;
    // number of non-linear iterations limit before time-step reduction for transport problem
    unsigned int nsa_limit_transp;    
    // maximum multiplication factor over time-step
    double dt_factor_max;
    // minimum multiplication factor over time-step    
    double dt_factor_min;
    // reduction factor over time-step
    double dt_factor_red;
    // time-step multiplication factor
    double mult_factor;
};

#endif /* TIMESTEPCONTROLRESIDUAL_H */