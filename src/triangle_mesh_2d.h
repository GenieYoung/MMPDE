#ifndef __MMPDE_TRIANGLE_MESH_2D__
#define __MMPDE_TRIANGLE_MESH_2D__

#include "matrix.h"
#include <CMTL/geo2d/geo2d_surface_mesh.h>
#include <CMTL/io/surface_mesh/write_obj.h>

namespace MMPDE
{
    class TriangleMesh2dTraits : public CMTL::geo2d::SurfaceMeshTraits
    {
    };

    typedef CMTL::geo2d::SurfaceMesh<real, TriangleMesh2dTraits> TriangleMesh2d;

    inline void build_triangle_mesh(TriangleMesh2d& mesh, const std::vector<real>& rows, const std::vector<real>& cols)
    {
        mesh.clear();

        assert(rows.size() > 1 && cols.size() > 1);

        std::vector<TriangleMesh2d::VertexHandle> vhs;
        vhs.reserve(rows.size() * cols.size());
        for(unsigned i = 0; i < rows.size(); ++i)
        {
            for(unsigned j = 0; j < cols.size(); ++j)
            {
                vhs.push_back(mesh.add_vertex(CMTL::geo2d::Point<real>(rows[i], cols[j])));
            }
        }

        for(unsigned i = 0; i < rows.size() - 1; ++i)
        {
            for(unsigned j = 0; j < cols.size() - 1; ++j)
            {
                TriangleMesh2d::VertexHandle v0 = vhs[i * cols.size() + j];
                TriangleMesh2d::VertexHandle v1 = vhs[i * cols.size() + j + 1];
                TriangleMesh2d::VertexHandle v2 = vhs[(i + 1) * cols.size() + j + 1];
                TriangleMesh2d::VertexHandle v3 = vhs[(i + 1) * cols.size() + j];
                mesh.add_face(v0, v1, v2);
                mesh.add_face(v0, v2, v3);
            }
        }
    }

    inline Matrix2d face_edge_matrix(const TriangleMesh2d& mesh, TriangleMesh2d::FaceHandle fh)
    {
        std::vector<CMTL::geo2d::Point<real>> pts;
        for(auto fv = mesh.fv_begin(fh); fv != mesh.fv_end(fh); ++fv)
        {
            pts.push_back(mesh.point(*fv));
        }
        return Matrix2d(pts[1].x() - pts[0].x(), pts[2].x() - pts[2].x(),
                        pts[1].y() - pts[0].y(), pts[2].y() - pts[0].y());
    }

    inline void write_obj(const TriangleMesh2d& mesh, const std::string& file)
    {
        CMTL::io::write_obj(mesh, file);
    }
}

#endif