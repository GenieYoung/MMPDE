#ifndef __MMPDE_FUNCTIONAL__
#define __MMPDE_FUNCTIONAL__

#include "matrix.h"
#include <tuple>

enum class Functional
{
    HUANG,
    WINSLOW
};

namespace MMPDE
{

inline void functional_Huang(const Matrix2d& J, real detJ, const Matrix2d& M, real& G, Matrix2d& GJ, real& GdetJ)
{
    unsigned d = 2;
    real theta = 1.0/3.0;
    real p = 1.5;
    real dp = p * d;

    real detM = M.det();
    Matrix2d M_inv = M.inverse();
    Matrix2d JT = J.transpose();
    real tr_JMJ = (J * M_inv * JT).trace();
 
    G = theta*detM*std::pow(tr_JMJ,dp/2) + (1-2*theta)*std::pow(d,dp/2)*std::pow(detM,1-p)*std::pow(detJ,p);
    GJ = M_inv*JT*dp*theta*detM*std::pow(tr_JMJ,dp/2-1);
    GdetJ = p*(1-2*theta)*std::pow(d,dp/2)*std::pow(detM,1-p)*std::pow(detJ,p-1);
}

inline void functional_Winslow(const Matrix2d& J, const Matrix2d& M, real&G, Matrix2d& GJ, real& GdetJ)
{
    Matrix2d M_inv = M.inverse();
    Matrix2d JT = J.transpose();

    G = 0.5 * (J * M_inv * JT).trace();
    GJ = M_inv * JT;
    GdetJ = 0;
}

}

#endif