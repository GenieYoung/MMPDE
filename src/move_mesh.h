#ifndef __MMPDE_MOVE_MESH__
#define __MMPDE_MOVE_MESH__

#include "trimesh.h"
#include "functional.h"
#include "metric.h"

// move_trimesh_2d.h

namespace MMPDE
{
    Trimesh2d move_mesh(const std::pair<double, double>& tspan,
                        const Trimesh2d& Xi_ref,
                        const Trimesh2d& X,
                        const std::vector<Matrix2d>& M,
                        real tau,
                        Functional Func);

    class MoveMeshRHS
    {
        public:
            MoveMeshRHS(const Trimesh2d& X, 
                        const std::vector<Matrix2d>& M,
                        real tau, 
                        Functional Func);

            void operator()(const std::vector<real>& xi,
                            std::vector<real>& dxidt,
                            const double t) const;

        private:
            std::vector<std::array<unsigned, 3>> _tris;
            std::vector<unsigned> _boundary_ids;
            std::vector<Matrix2d> _Mk;
            std::vector<Matrix2d> _E_inv;
            std::vector<real> _detE;
            std::vector<real> _volK;
            std::vector<real> _b_factor;
            real _sigma;
            real _tau;
            Functional _func;
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