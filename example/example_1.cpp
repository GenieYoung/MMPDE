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
    unsigned n = 10;
    unsigned times = 10;
    std::pair<real, real> tspan(0, 0.1);
    real tau = 0.01;

    // build initial mesh
    std::vector<real> rows(n+1), cols(n+1);
    real step = 1.0 / n;
    std::generate(rows.begin(), rows.end(), [start = -step, step]() mutable { return start+=step; });
    std::generate(cols.begin(), cols.end(), [start = -step, step]() mutable { return start+=step; });
    Trimesh2d ref_mesh(rows, cols);

    Trimesh2d mesh = ref_mesh;
    mesh.perturb(0.03);

    ref_mesh.export_vtk("ref_mesh.vtk");
    mesh.export_vtk("mesh.vtk");
    
    //mesh.set_value(U_func);

    mesh.compute_metric_tensor();

    Trimesh2d result_mesh;
    for(unsigned i = 0; i < times; ++i)
    {
        result_mesh = move_mesh(tspan, ref_mesh, mesh, tau, Functional::HUANG);
        result_mesh.export_vtk(std::string("result_") + std::to_string(i) + ".vtk");
    }

    return 0;
}