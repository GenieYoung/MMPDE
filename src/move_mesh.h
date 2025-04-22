#ifndef __MMPDE_MOVE_MESH__
#define __MMPDE_MOVE_MESH__

#include "trimesh.h"
#include "functional.h"

// move_trimesh_2d.h

namespace MMPDE
{
    Trimesh2d move_mesh(const std::pair<real, real>& tspan,
                        const Trimesh2d& Xi_ref,
                        const Trimesh2d& X,
                        real tau,
                        Functional Func);

    class MoveMeshRHS
    {
        public:
            MoveMeshRHS(const Trimesh2d& Xi_ref, const Trimesh2d& X, real tau, Functional Func);

            void operator()(const std::vector<real>& xi,
                            std::vector<real>& dxidt,
                            const double t) const;

        private:
            const Trimesh2d& _Xi_ref;
            const Trimesh2d& _X;
            real _tau;
            Functional _Func;
    };

    struct push_back_state_and_time
    {
        std::vector<std::vector<real>>& m_states;
        std::vector<double>& m_times;

        push_back_state_and_time(std::vector<std::vector<real>>& states,
                                 std::vector<double>& times)
            : m_states(states), m_times(times)
        {
        }

        void operator()(const std::vector<real>& x, double t)
        {
            m_states.push_back(x);
            m_times.push_back(t);
        }
    };
}

#endif