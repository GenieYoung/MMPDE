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

inline std::tuple<real, Matrix2d, Matrix2d> functional_Huang(const Matrix2d& J, real detJ, const Matrix2d& M)
{

}

inline std::tuple<real, Matrix2d, Matrix2d> functional_Winslow(const Matrix2d& J)
{

}

}

#endif