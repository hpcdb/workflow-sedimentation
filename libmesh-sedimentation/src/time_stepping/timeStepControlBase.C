/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   timeStepControlBase.cpp
 * Author: arossa
 * 
 * Created on 29 de Junho de 2016, 16:31
 */

#include "timeStepControlBase.h"
//#include "libmesh/libmesh.h"

using namespace libMesh;
using namespace std;

timeStepControlBase::timeStepControlBase( double init_dt, double dt_min, double dt_max, double dt_last, double dt_avg, unsigned int n_start_control, unsigned int max_nsa, std::string name ):
                                init_dt(init_dt),
                                dt_min(dt_min),
                                dt_max(dt_max),
                                dt_last(dt_last),
                                dt_avg(dt_avg),
                                start_control(n_start_control),
                                model_name(name),
                                nsa_max(max_nsa),
                                n_accepted_time_step(n_start_control),
                                n_rejected_time_step(0)
        
{
    cout<<"\nTime step control Base parameters defined\n"; 
}

timeStepControlBase::timeStepControlBase() {
    cout<<"construtor timeStepControlBase vazio\n";
}

timeStepControlBase::timeStepControlBase(const timeStepControlBase& orig) {
    cout<<"construtor timeStepControlBase de copia\n";
}

timeStepControlBase::~timeStepControlBase() {
}
    