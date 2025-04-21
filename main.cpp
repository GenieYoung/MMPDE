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
    std::vector<real> rows(50), cols(50);
    std::generate(rows.begin(), rows.end(), [start = 0.0]() mutable { return start+=0.02; });
    std::generate(cols.begin(), cols.end(), [start = 0.0]() mutable { return start+=0.02; });
    Trimesh2d mesh(rows, cols);
    mesh.set_value(func2);
    mesh.export_vtk("mesh.vtk");
    return 0;
}