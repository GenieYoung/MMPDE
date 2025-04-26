#include "move_mesh.h"

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
    std::pair<real, real> tspan(0, 1);
    real tau = 0.01;

    // build initial mesh
    std::vector<real> rows(n+1), cols(n+1);
    real step = 1.0 / n;
    std::generate(rows.begin(), rows.end(), [start = -step, step]() mutable { return start+=step; });
    std::generate(cols.begin(), cols.end(), [start = -step, step]() mutable { return start+=step; });
    Trimesh2d ref_mesh;
    construct_trimesh(ref_mesh, rows, cols);

    Trimesh2d mesh = ref_mesh;

    export_vtk(ref_mesh, "ref_mesh.vtk");

    for(unsigned i = 0; i < times; ++i)
    {
        std::vector<real> values = calc_vertex_value(mesh, U_func);
        std::vector<Matrix2d> M = calc_vertex_metric(mesh, values);
        export_vtk(mesh, (std::string("new_mesh_") + std::to_string(i) + ".vtk"), values);
        mesh = move_mesh(tspan, ref_mesh, mesh, M, tau, Functional::HUANG);
    }

    return 0;
}