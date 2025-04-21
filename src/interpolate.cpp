#include "interpolate.h"

namespace MMPDE
{
    Interpolate::Interpolate(const Trimesh2d& mesh)
    {
    }

    std::vector<real> Interpolate::operator()(const std::vector<real>& values, const std::vector<Point2d>& query_points) const
    {
        return {};
    }
}