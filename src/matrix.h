#ifndef __MMPDE_MATRIX__
#define __MMPDE_MATRIX__

#include <iostream>
#include <cassert>
#include <eigen3/Eigen/Dense>

namespace MMPDE
{

typedef double real;

struct Matrix2d
{
    explicit Matrix2d(real _a00=0, real _a01=0, real _a10=0, real _a11=0)
        : a00(_a00), a01(_a01), a10(_a10), a11(_a11) 
    {
    }

    Matrix2d(const Eigen::Matrix<real, 2, 2>& E)
        : Matrix2d(E(0,0), E(0,1), E(1,0), E(1,1))
    {
    }

    Matrix2d operator*(const Matrix2d& other) const
    {
        return Matrix2d(
            a00 * other.a00 + a01 * other.a10,
            a00 * other.a01 + a01 * other.a11,
            a10 * other.a00 + a11 * other.a10,
            a10 * other.a01 + a11 * other.a11);
    }

    Matrix2d operator*(real scalar) const
    {
        return Matrix2d(a00 * scalar, a01 * scalar, a10 * scalar, a11 * scalar);
    }

    Matrix2d operator+(const Matrix2d& other) const
    {
        return Matrix2d(a00 + other.a00, a01 + other.a01, a10 + other.a10, a11 + other.a11);
    }

    Matrix2d operator/(real scalar) const
    {
        assert(scalar != 0);
        return Matrix2d(a00 / scalar, a01 / scalar, a10 / scalar, a11 / scalar);
    }

    Matrix2d transpose() const
    {
        return Matrix2d(a00, a10, a01, a11);
    }

    real trace() const
    {
        return a00 + a11;
    }

    real det() const
    {
        return a00 * a11 - a01 * a10;
    }

    Matrix2d inverse() const
    {
        real d = det();
        //assert(d != 0);
        if(d == 0)
            return Matrix2d(std::numeric_limits<real>::max(), std::numeric_limits<real>::max(), 
                            std::numeric_limits<real>::max(), std::numeric_limits<real>::max());
        return Matrix2d(a11 / d, -a01 / d, -a10 / d, a00 / d);
    }

    Matrix2d sqrt() const
    {
        Eigen::Matrix<real, 2, 2> A;
        A << a00, a01, a10, a11;
        Eigen::Matrix<real, 2, 2> sqrtA = A.llt().matrixL();
        return Matrix2d(sqrtA);
    }

    friend std::ostream& operator<<(std::ostream& os, const Matrix2d& m)
    {
        os << "[" << m.a00 << ", " << m.a01 << "]\n"
           << "[" << m.a10 << ", " << m.a11 << "]";
        return os;
    }
    
    real a00, a01, a10, a11;
};

}

#endif