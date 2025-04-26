#include "trimesh.h"
#include <fstream>

namespace MMPDE
{
    void construct_trimesh(Trimesh2d& mesh, const std::vector<real>& rows, const std::vector<real>& cols)
    {
        mesh.clear();
        assert(rows.size() > 1 && cols.size() > 1);

        std::vector<Trimesh2d::VertexHandle> vhs;
        vhs.reserve(rows.size() * cols.size());
        for(unsigned i = 0; i < cols.size(); ++i)
        {
            for(unsigned j = 0; j < rows.size(); ++j)
            {
                vhs.push_back(mesh.add_vertex(CMTL::geo2d::Point<real>(rows[j], cols[i])));
            }
        }

        for(unsigned i = 0; i < rows.size() - 1; ++i)
        {
            for(unsigned j = 0; j < cols.size() - 1; ++j)
            {
                Trimesh2d::VertexHandle v0 = vhs[i * rows.size() + j];
                Trimesh2d::VertexHandle v1 = vhs[i * rows.size() + j + 1];
                Trimesh2d::VertexHandle v2 = vhs[(i + 1) * rows.size() + j + 1];
                mesh.add_face(v0, v1, v2);
            }
        }

        for(unsigned i = 0; i < rows.size() - 1; ++i)
        {
            for(unsigned j = 0; j < cols.size() - 1; ++j)
            {
                Trimesh2d::VertexHandle v0 = vhs[i * rows.size() + j];
                Trimesh2d::VertexHandle v2 = vhs[(i + 1) * rows.size() + j + 1];
                Trimesh2d::VertexHandle v3 = vhs[(i + 1) * rows.size() + j];
                mesh.add_face(v0, v2, v3);
            }
        }

        // for(unsigned i = 0; i < rows.size() - 1; ++i)
        // {
        //     for(unsigned j = 0; j < cols.size() - 1; ++j)
        //     {
        //         Trimesh2d::VertexHandle v0 = vhs[i * rows.size() + j];
        //         Trimesh2d::VertexHandle v1 = vhs[i * rows.size() + j + 1];
        //         Trimesh2d::VertexHandle v2 = vhs[(i + 1) * rows.size() + j + 1];
        //         Trimesh2d::VertexHandle v3 = vhs[(i + 1) * rows.size() + j];
        //         mesh.add_face(v0, v1, v2);
        //         mesh.add_face(v0, v2, v3);
        //     }
        // }
    }

    void update_trimesh(Trimesh2d& mesh, const std::vector<Point2d>& new_vertices)
    {
        assert(new_vertices.size() == mesh.n_vertices());
        for(unsigned i = 0; i < new_vertices.size(); ++i)
        {
            mesh.point(mesh.vertex_handle(i)) = new_vertices[i];
        }
    }

    std::vector<Point2d> get_vertices(const Trimesh2d& mesh)
    {
        std::vector<Point2d> vertices;
        vertices.reserve(mesh.n_vertices());
        for(auto vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit)
            vertices.push_back(mesh.point(*vit));
        return vertices;
    }

    std::vector<std::array<unsigned, 3>> get_faces(const Trimesh2d& mesh)
    {
        std::vector<std::array<unsigned, 3>> faces;
        faces.reserve(mesh.n_faces());
        for(auto fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit)
        {
            std::array<unsigned, 3> face;
            unsigned i = 0;
            for(auto fvit = mesh.fv_begin(*fit); fvit != mesh.fv_end(*fit); ++fvit)
            {
                face[i++] = fvit->idx();
            }
            faces.push_back(face);
        }
        return faces;
    }

    std::vector<unsigned> get_boundary_ids(const Trimesh2d& mesh)
    {
        std::vector<unsigned> boundary_ids;
        for(auto vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit)
        {
            if(mesh.is_boundary(*vit))
                boundary_ids.push_back(vit->idx());
        }
        return boundary_ids;
    }

    std::vector<Point2d> calc_face_barycenter(const Trimesh2d& mesh)
    {
        std::vector<Point2d> barycenters(mesh.n_faces());
        for(auto fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit)
        {
            barycenters[fit->idx()] = mesh.barycenter(*fit);
        }
        return barycenters;
    }

    std::vector<Matrix2d> calc_face_edge_matrix(const Trimesh2d& mesh)
    {
        std::vector<Matrix2d> face_edge_matrix(mesh.n_faces());
        for(auto fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit)
        {
            std::vector<Point2d> pts;
            for(auto fvit = mesh.fv_begin(*fit); fvit != mesh.fv_end(*fit); ++fvit)
            {
                pts.push_back(mesh.point(*fvit));
            }
            face_edge_matrix[fit->idx()] = Matrix2d(pts[1].x() - pts[0].x(), pts[1].y() - pts[0].y(),
                                                    pts[2].x() - pts[0].x(), pts[2].y() - pts[0].y());
        }
        return face_edge_matrix;
    }

    void perturb(Trimesh2d& mesh, real degree)
    {
        for(auto vh = mesh.vertices_begin(); vh != mesh.vertices_end(); ++vh)
        {
            if(mesh.is_boundary(*vh))
                continue;
            Point2d p = mesh.point(*vh);
            p.x() += degree * (rand() % 200 / 100.0 - 1.0);
            p.y() += degree * (rand() % 200 / 100.0 - 1.0);
            mesh.point(*vh) = p;
        }
    }

    void export_vtk(const Trimesh2d& mesh, const std::string& file, const std::vector<real>& values)
    {
        std::ofstream fout;
        fout.open(file.c_str(), std::ofstream::trunc);

        fout << "# vtk DataFile Version 3.0" <<'\n';
        fout << "Date calculated by OctreeMesh"    <<'\n';
        fout << "ASCII"                      <<'\n';
        fout << "DATASET UNSTRUCTURED_GRID"  <<'\n' << '\n';

        fout << "POINTS " << mesh.n_vertices() << " double" << '\n';
        for(auto vh = mesh.vertices_begin(); vh != mesh.vertices_end(); ++vh)
        {
            fout << mesh.point(*vh).x() << " " << mesh.point(*vh).y() << " 0" << '\n';
        }
        fout << '\n';

        fout << "CELLS " << mesh.n_faces() << " " << mesh.n_faces() * 4 << '\n';
        for(auto fh = mesh.faces_begin(); fh != mesh.faces_end(); ++fh)
        {
            fout << "3 ";
            for(auto fv = mesh.fv_begin(*fh); fv != mesh.fv_end(*fh); ++fv)
            {
                fout << fv->idx() << " ";
            }
            fout << '\n';
        }
        fout << '\n';

        fout << "CELL_TYPES " << mesh.n_faces() << '\n';
        for ( unsigned i = 0; i < mesh.n_faces(); ++i )
            fout << 5 << std::endl; 
    
        fout << std::endl;

        if(values.size() == mesh.n_vertices())
        {
            fout << "POINT_DATA " << mesh.n_vertices() << '\n';
            fout << "SCALARS value double" << '\n';
            fout << "LOOKUP_TABLE default" << '\n';
            for(auto vh = mesh.vertices_begin(); vh != mesh.vertices_end(); ++vh)
            {
                //fout << mesh.attribute(*vh) << '\n';
                fout << values[vh->idx()] << '\n';
            }
            fout << '\n';
        }

        fout.close();
    }
}