/* 
 * File:   stab_helper.h
 * Author: camata
 *
 * Created on January 26, 2017, 2:33 PM
 */

#ifndef STAB_HELPER_H
#define	STAB_HELPER_H

#include "libmesh/libmesh.h"
#include "libmesh/fe_base.h"
#include "libmesh/elem.h"
#include "libmesh/vector_value.h"
#include "libmesh/numeric_vector.h"
using namespace libMesh;

RealGradient compute_g(FEBase* fe,
    const unsigned short dim,
    const unsigned int qp);

// G for stabilization tau
RealTensor compute_G(FEBase* fe,
    const unsigned short dim,
    const unsigned int qp);

Real compute_tau_M(RealGradient g, RealTensor G, Gradient U, Real k, Real dt, Real dt_stab);

Real compute_tau_C(RealGradient g, Real tauM);

#endif	/* STAB_HELPER_H */

