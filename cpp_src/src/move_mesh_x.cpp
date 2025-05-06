#include "move_mesh_x.h"
#include "functional.h"
#include <boost/numeric/odeint.hpp>

namespace MMPDE
{
    Trimesh2d move_mesh_x(const std::pair<real, real>& tspan,
        const Trimesh2d& X,
        real tau)
    {
        return X;
    }
}