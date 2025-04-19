#include "triangle_mesh_2d.h"

int main()
{
    MMPDE::TriangleMesh2d mesh;
    MMPDE::build_triangle_mesh(mesh, {0, 1, 2, 3}, {0, 1, 2, 3, 4});
    MMPDE::write_obj(mesh, "mesh.obj");
    std::cout << MMPDE::face_edge_matrix(mesh, mesh.face_handle(1));
    return 0;
}