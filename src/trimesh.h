#ifndef __MMPDE_TRIMESH__
#define __MMPDE_TRIMESH__

#include "matrix.h"
#include <CMTL/geo2d/geo2d_surface_mesh.h>
#include <CMTL/io/surface_mesh/write_obj.h>
#include <functional>

namespace MMPDE
{
    typedef CMTL::geo2d::Point<real> Point2d;

    class Trimesh2d
    {
        public:
            Trimesh2d() = default;
            Trimesh2d(const std::vector<real>& rows, const std::vector<real>& cols);

        public:
            unsigned n_vertices() const { return _mesh.n_vertices(); }

            unsigned n_faces() const { return _mesh.n_faces(); }

            std::vector<real> get_vertices() const;

            std::vector<std::array<unsigned, 3>> get_faces() const;

            bool is_boundary(unsigned vid) const;

            const Matrix2d& get_metric(unsigned vid) const;

            const Point2d& get_vertex(unsigned vid) const;

            void set_vertex(unsigned vid, const Point2d& p);

            void set_value(const std::function<real(const Point2d&)>& func);

            void compute_metric_tensor();

            void perturb(real eps);

            Matrix2d face_edge_matrix(unsigned fid) const;

            Matrix2d average_metric(unsigned fid) const;

            Point2d barycenter(unsigned fid) const;

            void export_obj(const std::string& file) const;

            void export_vtk(const std::string& file) const;

        public:
            std::vector<real> _values;

            std::vector<Matrix2d> _metrics;

            // struct Traits : public CMTL::geo2d::SurfaceMeshTraits 
            // {
            //     typedef real VertexAttribute;
            // };
            typedef CMTL::geo2d::SurfaceMesh<real> MeshType;
            typedef MeshType::VertexHandle VertexHandle;
            typedef MeshType::FaceHandle FaceHandle;
            MeshType _mesh;
    };
}

#endif