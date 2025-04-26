#include "interpolate.h"
#include "CMTL/geometry_algorithm/predicate.h"

namespace MMPDE
{
    Interpolate::Interpolate(const Trimesh2d& mesh) : Interpolate(get_vertices(mesh), get_faces(mesh))
    {
    }

    Interpolate::Interpolate(const std::vector<Point2d>& nodes, const std::vector<std::array<unsigned, 3>>& faces)
    {
        _nodes = nodes;
        _faces = faces;
        
        for(unsigned i = 0; i < _nodes.size(); ++i)
            _nodes_map.insert({_nodes[i], i});

        for(unsigned i = 0; i < _faces.size(); ++i)
        {
            std::array<unsigned, 3>& tri = _faces[i];
            real x_min = std::min({_nodes[tri[0]].x(), _nodes[tri[1]].x(), _nodes[tri[2]].x()});
            real x_max = std::max({_nodes[tri[0]].x(), _nodes[tri[1]].x(), _nodes[tri[2]].x()});
            real y_min = std::min({_nodes[tri[0]].y(), _nodes[tri[1]].y(), _nodes[tri[2]].y()});
            real y_max = std::max({_nodes[tri[0]].y(), _nodes[tri[1]].y(), _nodes[tri[2]].y()});
            real min_loc[2] = {x_min, y_min};
            real max_loc[2] = {x_max, y_max};
            _rtree.Insert(min_loc, max_loc, (void*)(&tri));
        }
    }

    template<typename T>
    std::vector<T> Interpolate::operator()(const std::vector<T>& values, const std::vector<Point2d>& query_points) const
    {
        std::vector<T> results(query_points.size());

        for(unsigned i = 0; i < query_points.size(); ++i)
        {
            const Point2d& qp = query_points[i];
            auto find_exact = _nodes_map.find(qp);
            if(find_exact != _nodes_map.end())
            {
                results[i] = values[find_exact->second];
            }
            else
            {
                real min_loc[2] = {qp.x(), qp.y()};
                real max_loc[2] = {qp.x(), qp.y()};
                std::vector<std::array<unsigned,3>*> candidate_tris;
                _rtree.Search(min_loc, max_loc, _candidate_search_call_back, &candidate_tris);
                assert(!candidate_tris.empty());
                for(unsigned j = 0; j < candidate_tris.size(); ++j)
                {
                    unsigned tid0 = (*candidate_tris[j])[0];
                    unsigned tid1 = (*candidate_tris[j])[1];
                    unsigned tid2 = (*candidate_tris[j])[2];
                    const Point2d& t0 = _nodes[tid0];
                    const Point2d& t1 = _nodes[tid1]; 
                    const Point2d& t2 = _nodes[tid2];
                    if(CMTL::geometry_algorithm::in_triangle(t0, t1, t2, qp) == CMTL::ORIENTATION::INSIDE)
                    {
                        real x10 = t1.x() - t0.x();
                        real y10 = t1.y() - t0.y();
                        real x20 = t2.x() - t0.x();
                        real y20 = t2.y() - t0.y();
                        real xq0 = qp.x() - t0.x();
                        real yq0 = qp.y() - t0.y();
                        real det = x10 * y20 - x20 * y10;
                        real w1 = (xq0 * y20 - x20 * yq0) / det;
                        real w2 = (x10 * yq0 - xq0 * y10) / det;
                        real w0 = 1 - w1 - w2;
                        results[i] = values[tid0] * w0 + values[tid1] * w1 + values[tid2] * w2;
                        break;
                    }
                }
            }
        }

        return results;
    }

    bool Interpolate::_candidate_search_call_back(void* mp, void* arg)
    {
        std::vector<std::array<unsigned,3>*>* candidates = (std::vector<std::array<unsigned,3>*>*)arg;
        candidates->push_back((std::array<unsigned,3>*)mp);
        return true;
    }

    template std::vector<real> Interpolate::operator()(const std::vector<real>& values, const std::vector<Point2d>& query_points) const;
    template std::vector<Point2d> Interpolate::operator()(const std::vector<Point2d>& values, const std::vector<Point2d>& query_points) const;
}