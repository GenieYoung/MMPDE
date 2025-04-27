#include "move_mesh.h"

using namespace MMPDE;

int main()
{
    // set basic parameters
    unsigned n = 20;
    std::pair<double, double> tspan(0, 0.1);
    real tau = 0.01;

    // build initial mesh
    std::vector<real> rows(n+1), cols(n+1);
    real step = 1.0 / n;
    std::generate(rows.begin(), rows.end(), [start = -step, step]() mutable { return start+=step; });
    std::generate(cols.begin(), cols.end(), [start = -step, step]() mutable { return start+=step; });
    Trimesh2d ref_mesh;
    construct_trimesh(ref_mesh, rows, cols);

    Trimesh2d mesh = ref_mesh;
    //perturb(mesh, step / 5.0);
    perturb(mesh, 0.01);

    export_vtk(ref_mesh, "ref_mesh.vtk");
    export_vtk(mesh, "mesh.vtk");

    std::vector<Matrix2d> M = calc_vertex_metric(mesh, MetricType::IDENTITY);

    unsigned times = 10;
    for(unsigned i = 0; i < times; ++i)
    {
        export_vtk(mesh, (std::string("new_mesh_") + std::to_string(i) + ".vtk"));
        mesh = move_mesh(tspan, ref_mesh, mesh, M, tau, Functional::HUANG);
    }

    return 0;
}