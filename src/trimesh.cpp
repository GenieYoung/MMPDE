#include "trimesh.h"

namespace MMPDE
{
    Trimesh2d::Trimesh2d(const std::vector<real>& rows, const std::vector<real>& cols)
    {
        _mesh.clear();

        assert(rows.size() > 1 && cols.size() > 1);

        std::vector<VertexHandle> vhs;
        vhs.reserve(rows.size() * cols.size());
        for(unsigned i = 0; i < cols.size(); ++i)
        {
            for(unsigned j = 0; j < rows.size(); ++j)
            {
                vhs.push_back(_mesh.add_vertex(CMTL::geo2d::Point<real>(rows[j], cols[i])));
            }
        }

        for(unsigned i = 0; i < rows.size() - 1; ++i)
        {
            for(unsigned j = 0; j < cols.size() - 1; ++j)
            {
                VertexHandle v0 = vhs[i * rows.size() + j];
                VertexHandle v1 = vhs[i * rows.size() + j + 1];
                VertexHandle v2 = vhs[(i + 1) * rows.size() + j + 1];
                VertexHandle v3 = vhs[(i + 1) * rows.size() + j];
                _mesh.add_face(v0, v1, v2);
                _mesh.add_face(v0, v2, v3);
            }
        }
    }

    std::vector<real> Trimesh2d::get_vertices() const
    {
        std::vector<real> vertices;
        vertices.reserve(_mesh.n_vertices() * 2);
        for(auto vh = _mesh.vertices_begin(); vh != _mesh.vertices_end(); ++vh)
        {
            vertices.push_back(_mesh.point(*vh).x());
            vertices.push_back(_mesh.point(*vh).y());
        }
        return vertices;
    }

    std::vector<std::array<unsigned, 3>> Trimesh2d::get_faces() const
    {
        std::vector<std::array<unsigned, 3>> faces;
        faces.reserve(_mesh.n_faces());
        for(auto fh = _mesh.faces_begin(); fh != _mesh.faces_end(); ++fh)
        {
            std::array<unsigned, 3> face;
            unsigned i = 0;
            for(auto fv = _mesh.fv_begin(*fh); fv != _mesh.fv_end(*fh); ++fv)
            {
                face[i++] = fv->idx();
            }
            faces.push_back(face);
        }
        return faces;
    }

    bool Trimesh2d::is_boundary(unsigned vid) const
    {
        return _mesh.is_boundary(_mesh.vertex_handle(vid));
    }

    const Matrix2d& Trimesh2d::get_metric(unsigned vid) const
    {
        return _metrics[vid];
    }

    void Trimesh2d::set_vertex(unsigned vid, const Point2d& p)
    {
        _mesh.point(_mesh.vertex_handle(vid)) = p;
    }

    void Trimesh2d::set_value(const std::function<real(const Point2d&)>& func)
    {
        _values.resize(_mesh.n_vertices());
        for(auto vh = _mesh.vertices_begin(); vh != _mesh.vertices_end(); ++vh)
        {
            _values[vh->idx()] = func(_mesh.point(*vh));
        }
    }

    void Trimesh2d::compute_metric_tensor()
    {
        _metrics.clear();
        _metrics.resize(_mesh.n_vertices(), Matrix2d(1,0,0,1));
    }

    Matrix2d Trimesh2d::face_edge_matrix(unsigned fid) const
    {
        assert(fid < _mesh.n_faces());
        FaceHandle fh = _mesh.face_handle(fid);
        std::vector<CMTL::geo2d::Point<real>> pts;
        for(auto fv = _mesh.fv_begin(fh); fv != _mesh.fv_end(fh); ++fv)
        {
            pts.push_back(_mesh.point(*fv));
        }
        return Matrix2d(pts[1].x() - pts[0].x(), pts[2].x() - pts[0].x(),
                        pts[1].y() - pts[0].y(), pts[2].y() - pts[0].y());
    }

    void Trimesh2d::perturb(real eps)
    {
        for(auto vh = _mesh.vertices_begin(); vh != _mesh.vertices_end(); ++vh)
        {
            if(_mesh.is_boundary(*vh))
                continue;
            Point2d p = _mesh.point(*vh);
            p.x() += eps * (rand() % 100 / 100.0);
            p.y() += eps * (rand() % 100 / 100.0);
            _mesh.point(*vh) = p;
        }
    }

    Matrix2d Trimesh2d::average_metric(unsigned fid) const
    {
        assert(fid < _mesh.n_faces());
        FaceHandle fh = _mesh.face_handle(fid);
        Matrix2d m;
        for(auto fv = _mesh.fv_begin(fh); fv != _mesh.fv_end(fh); ++fv)
        {
            m = m + _metrics[fv->idx()];
        }
        m = m / 3.0;
        return m;
    }

    Point2d Trimesh2d::barycenter(unsigned fid) const
    {
        return _mesh.barycenter(_mesh.face_handle(fid));
    }

    void Trimesh2d::export_obj(const std::string& file) const
    {
        CMTL::io::write_obj(_mesh, file);
    }

    void Trimesh2d::export_vtk(const std::string& file) const
    {
        std::ofstream fout;
        fout.open(file.c_str(), std::ofstream::trunc);

        fout << "# vtk DataFile Version 3.0" <<'\n';
        fout << "Date calculated by OctreeMesh"    <<'\n';
        fout << "ASCII"                      <<'\n';
        fout << "DATASET UNSTRUCTURED_GRID"  <<'\n' << '\n';

        fout << "POINTS " << _mesh.n_vertices() << " double" << '\n';
        for(auto vh = _mesh.vertices_begin(); vh != _mesh.vertices_end(); ++vh)
        {
            fout << _mesh.point(*vh).x() << " " << _mesh.point(*vh).y() << " 0" << '\n';
        }
        fout << '\n';

        fout << "CELLS " << _mesh.n_faces() << " " << _mesh.n_faces() * 4 << '\n';
        for(auto fh = _mesh.faces_begin(); fh != _mesh.faces_end(); ++fh)
        {
            fout << "3 ";
            for(auto fv = _mesh.fv_begin(*fh); fv != _mesh.fv_end(*fh); ++fv)
            {
                fout << fv->idx() << " ";
            }
            fout << '\n';
        }
        fout << '\n';

        fout << "CELL_TYPES " << _mesh.n_faces() << '\n';
        for ( unsigned i = 0; i < _mesh.n_faces(); ++i )
            fout << 5 << std::endl; 
    
        fout << std::endl;

        if(_values.size() == _mesh.n_vertices())
        {
            fout << "POINT_DATA " << _mesh.n_vertices() << '\n';
            fout << "SCALARS value double" << '\n';
            fout << "LOOKUP_TABLE default" << '\n';
            for(auto vh = _mesh.vertices_begin(); vh != _mesh.vertices_end(); ++vh)
            {
                //fout << _mesh.attribute(*vh) << '\n';
                fout << _values[vh->idx()] << '\n';
            }
            fout << '\n';
        }

        fout.close();
    }
}