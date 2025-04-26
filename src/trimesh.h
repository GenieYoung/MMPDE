#ifndef __MMPDE_TRIMESH__
#define __MMPDE_TRIMESH__

#include "matrix.h"
#include <CMTL/geo2d/geo2d_surface_mesh.h>

namespace MMPDE
{
    typedef CMTL::geo2d::Point<real> Point2d;
    typedef CMTL::geo2d::SurfaceMesh<real> Trimesh2d;

    void construct_trimesh(Trimesh2d& mesh, const std::vector<real>& rows, const std::vector<real>& cols);

    void update_trimesh(Trimesh2d& mesh, const std::vector<Point2d>& new_vertices);

    std::vector<Point2d> get_vertices(const Trimesh2d& mesh);
    
    std::vector<std::array<unsigned, 3>> get_faces(const Trimesh2d& mesh);

    std::vector<unsigned> get_boundary_ids(const Trimesh2d& mesh);

    std::vector<Point2d> calc_face_barycenter(const Trimesh2d& mesh);

    std::vector<Matrix2d> calc_face_edge_matrix(const Trimesh2d& mesh);

    void perturb(Trimesh2d& mesh, real degree);

    void export_vtk(const Trimesh2d& mesh, const std::string& file, const std::vector<real>& values = {});
}

#endif