#ifndef __MMPDE_METRIC__
#define __MMPDE_METRIC__

#include "trimesh.h"
#include <functional>

namespace MMPDE
{
    enum class MetricType
    {
        IDENTITY,
        CURVATURE
    };

    std::vector<real> calc_vertex_value(const Trimesh2d& mesh, const std::function<real(const Point2d&)>& func);

    std::vector<Matrix2d> calc_vertex_metric(const Trimesh2d& mesh, const std::function<Matrix2d(const Point2d&)>& func);

    std::vector<Matrix2d> calc_vertex_metric(const Trimesh2d& mesh, MetricType metric_type);

    std::vector<Matrix2d> calc_vertex_metric(const Trimesh2d& mesh, const std::vector<real>& values);
}

#endif
