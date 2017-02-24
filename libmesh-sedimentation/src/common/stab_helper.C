// g for stabilization tau


#include "stab_helper.h"


RealGradient compute_g(FEBase* fe,
    const unsigned short dim,
    const unsigned int qp)
{
    RealGradient g(fe->get_dxidx()[qp]);

    if (dim > 1)
    {
        g(0) += fe->get_detadx()[qp];
        g(1) = fe->get_dxidy()[qp] + fe->get_detady()[qp];
    }

    if (dim == 3)
    {
        g(0) += fe->get_dzetadx()[qp];
        g(1) += fe->get_dzetady()[qp];
        g(2) = fe->get_dxidz()[qp] + fe->get_detadz()[qp] + fe->get_dzetadz()[qp];
    }

    return g;
}


// G for stabilization tau
RealTensor compute_G(FEBase* fe,
    const unsigned short dim,
    const unsigned int qp)
{
    Real dxidx = fe->get_dxidx()[qp];
    RealTensor G(dxidx*dxidx);

    if (dim > 1)
    {
        Real dxidy = fe->get_dxidy()[qp];
        Real detadx = fe->get_detadx()[qp];
        Real detady = fe->get_detady()[qp];

        G(0, 0) += detadx*detadx;
        G(0, 1) = G(1, 0) = dxidx*dxidy + detadx*detady;
        G(1, 1) = dxidy*dxidy + detady*detady;

        if (dim == 3)
        {
            Real dxidz = fe->get_dxidz()[qp];
            Real detadz = fe->get_detadz()[qp];
            Real dzetadx = fe->get_dzetadx()[qp];
            Real dzetady = fe->get_dzetady()[qp];
            Real dzetadz = fe->get_dzetadz()[qp];

            G(0, 0) += dzetadx*dzetadx;
            G(0, 1) += dzetadx*dzetady;
            G(0, 2) = dxidx*dxidz + detadx*detadz + dzetadx*dzetadz;
            G(1, 0) += dzetady*dzetadx;
            G(1, 1) += dzetady*dzetady;
            G(1, 2) = dxidy*dxidz + detady*detadz + dzetady*dzetadz;
            G(2, 0) = dxidz*dxidx + detadz*detadx + dzetadz*dzetadx;
            G(2, 1) = dxidz*dxidy + detadz*detady + dzetadz*dzetady;
            G(2, 2) = dxidz*dxidz + detadz*detadz + dzetadz*dzetadz;
        }
    }

    return G;
}

Real compute_tau_M(RealGradient g, RealTensor G, Gradient U, Real k, Real dt, Real dt_stab)
{
    Real tau = (U)*(G*U) + (k*k)*(G.contract(G)) + dt_stab*4.0/(dt*dt);
    return 1.0/std::sqrt(tau);
    
}


Real compute_tau_C(RealGradient g, Real tauM)
{
    Real tauC = g*(tauM*g);
    return 1.0/tauC;
    
}
