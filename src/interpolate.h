#ifndef __MMPDE_INTERPOLATE__
#define __MMPDE_INTERPOLATE__

#include "trimesh.h"

namespace MMPDE
{
    class Interpolate
    {
        public:
            Interpolate(const Trimesh2d& mesh);

        public:
            std::vector<real> operator()(const std::vector<real>& values, const std::vector<Point2d>& query_points) const;

            std::vector<Point2d> operator()(const std::vector<Point2d>& positions, const std::vector<Point2d>& query_points) const;

        private:
            const Trimesh2d& _mesh;
    };
}

#endif
