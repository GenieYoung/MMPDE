#include "interpolate.h"

namespace MMPDE
{
    Interpolate::Interpolate(const Trimesh2d& mesh) : _mesh(mesh)
    {
    }

    std::vector<real> Interpolate::operator()(const std::vector<real>& values, const std::vector<Point2d>& query_points) const
    {
        return {};
    }

    std::vector<Point2d> Interpolate::operator()(const std::vector<Point2d>& positions, const std::vector<Point2d>& query_points) const
    {
        
    }
}