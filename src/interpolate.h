#ifndef __MMPDE_INTERPOLATE__
#define __MMPDE_INTERPOLATE__

#include "trimesh.h"
#include "RTree.h"
#include <set>

namespace MMPDE
{
    class Interpolate
    {
        public:
            Interpolate(const Trimesh2d& mesh);

            Interpolate(const std::vector<Point2d>& nodes, const std::vector<std::array<unsigned, 3>>& faces);

        public:
            template<typename T>
            std::vector<T> operator()(const std::vector<T>& values, const std::vector<Point2d>& query_points) const; 

        private:
            std::vector<Point2d> _nodes;

            std::map<Point2d, unsigned> _nodes_map;

            std::vector<std::array<unsigned, 3>> _faces;

            typedef RTree<void*, real, 2, real, 8, 4> RTriangleTree;
            mutable RTriangleTree _rtree;

            static bool _candidate_search_call_back(void* mp, void* arg);
    };
}

#endif
