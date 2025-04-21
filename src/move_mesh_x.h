#ifndef __MMPDE_MOVE_MESH_X__
#define __MMPDE_MOVE_MESH_X__

#include "trimesh.h"

namespace MMPDE
{
    Trimesh2d move_mesh_x(const std::pair<real, real>& tspan,
                          const Trimesh2d& X,
                          real tau);
}

#endif