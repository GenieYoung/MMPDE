#include "move_mesh.h"
#include <boost/numeric/odeint.hpp>

namespace MMPDE
{
    Trimesh2d move_mesh(const std::pair<real, real>& tspan,
        const Trimesh2d& Xi_ref,
        const Trimesh2d& X,
        real tau,
        Functional G)
    {
        MoveMeshRHS rhs(Xi_ref, X, tau, G);

        std::vector<real> xi = Xi_ref.get_vertices();
        
        double dt = (tspan.second - tspan.first) / 10.0;
        std::vector<std::vector<real>> xi_vec;
        std::vector<double> times;
        boost::numeric::odeint::integrate(rhs, xi, tspan.first, tspan.second, dt, push_back_state_and_time(xi_vec, times));
        return X;
    }

    MoveMeshRHS::MoveMeshRHS(const Trimesh2d& Xi_ref, const Trimesh2d& X, real tau, Functional G)
        : _Xi_ref(Xi_ref), _X(X), _tau(tau), _G(G)
    {
        assert(_X.n_faces() == _Xi_ref.n_faces());
    }

    void MoveMeshRHS::operator()(const std::vector<real>& xi,
                                 std::vector<real>& dxidt,
                                 const double t) const
    {
        dxidt.resize(xi.size());

        std::vector<Matrix2d> E_inv(_X.n_faces());
        std::vector<real> detE(_X.n_faces());
        std::vector<Point2d> barycenters(_X.n_faces());
        std::vector<real> volK(_X.n_faces());

        std::vector<Matrix2d> Ec(_Xi_ref.n_faces());
        std::vector<Matrix2d> Ec_inv(_Xi_ref.n_faces());
        std::vector<real> detEc(_Xi_ref.n_faces());

        std::vector<Matrix2d> J(_X.n_faces());
        std::vector<Matrix2d> detJ(_X.n_faces());

        std::vector<Matrix2d> Mk(_X.n_faces());

        real sigma = 0;

        for(unsigned i = 0; i < _X.n_faces(); ++i)
        { 
            Matrix2d E = _X.face_edge_matrix(i);
            E_inv[i] = E.inverse();
            detE[i] = std::abs(E.det());
            barycenters[i] = _X.barycenter(i);
            volK[i] = 0.5 * detE[i];
            Mk[i] = _X.average_metric(i);
            sigma += volK[i] * Mk[i].det();
        }

        for(unsigned i = 0; i < _Xi_ref.n_faces(); ++i)
        {
            Matrix2d E = _X.face_edge_matrix(i);
            Ec[i] = E;
            Ec_inv[i] = E.inverse();
            detEc[i] = std::abs(E.det());
            J[i] = Ec[i] * E_inv[i];
            detJ[i] = detEc[i] / detE[i];
        }
    }
}