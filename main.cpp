#include "trimesh.h"

using namespace MMPDE;

real func(const Point2d& p)
{
    return p.x() * p.x() + p.y() * p.y();
}

real func2(const Point2d& p)
{
    return tanh(-30*(p.y()-0.5-0.25*sin(2*M_PI*p.x())));
}

int main()
{
    return 0;
}