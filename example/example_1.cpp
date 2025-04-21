#include "move_mesh.h"
#include "functional"

using namespace MMPDE;

real U_func(const Point2d& p)
{
    return tanh(-30*(p.y()-0.5-0.25*sin(2*M_PI*p.x())));
}

int main()
{
    // set basic parameters
    unsigned n = 1;
    unsigned times = 10;
    std::pair<real, real> tspan(0, 0.1);
    real tau = 0.01;

    // build initial mesh
    std::vector<real> rows(n+1), cols(n+1);
    real step = 1.0 / n;
    std::generate(rows.begin(), rows.end(), [start = -step, step]() mutable { return start+=step; });
    std::generate(cols.begin(), cols.end(), [start = -step, step]() mutable { return start+=step; });
    Trimesh2d mesh(rows, cols);
    
    mesh.set_value(U_func);

    mesh.compute_metric_tensor();

    for(unsigned i = 0; i < times; ++i)
    {
        mesh = move_mesh(tspan, mesh, mesh, tau, Functional::HUANG);
    }

    mesh.export_vtk("mesh.vtk");
    return 0;
}