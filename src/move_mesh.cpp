#include "move_mesh.h"
#include "interpolate.h"
#include <boost/numeric/odeint.hpp>

namespace MMPDE
{
    std::vector<real> reshape(const std::vector<Point2d>& points);
    std::vector<Point2d> reshape(const std::vector<real>& point_list);

    Trimesh2d move_mesh(const std::pair<real, real>& tspan,
        const Trimesh2d& Xi_ref,
        const Trimesh2d& X,
        real tau,
        Functional Func)
    {
        MoveMeshRHS rhs(Xi_ref, X, tau, Func);

        std::vector<real> xi = Xi_ref.get_vertices();
        
        double dt = (tspan.second - tspan.first) / 10.0;
        std::vector<std::vector<real>> xi_vec;
        std::vector<double> times;
        boost::numeric::odeint::integrate(rhs, xi, tspan.first, tspan.second, dt, push_back_state_and_time(xi_vec, times));

        Interpolate intp(reshape(xi), Xi_ref.get_faces());
        std::vector<Point2d> new_X_pts = intp(reshape(X.get_vertices()), reshape(Xi_ref.get_vertices()));

        Trimesh2d new_mesh = X;
        for(unsigned i = 0; i < X.n_vertices(); ++i)
        {
            new_mesh.set_vertex(i, Point2d(xi[2*i], xi[2*i+1]));
        }

        return new_mesh;
    }

    MoveMeshRHS::MoveMeshRHS(const Trimesh2d& Xi_ref, const Trimesh2d& X, real tau, Functional Func)
        : _Xi_ref(Xi_ref), _X(X), _tau(tau), _Func(Func)
    {
        assert(_X.n_faces() == _Xi_ref.n_faces());
    }

    void MoveMeshRHS::operator()(const std::vector<real>& xi,
                                 std::vector<real>& dxidt,
                                 const double t) const
    {
        dxidt.resize(xi.size(), 0);

        std::vector<std::array<unsigned, 3>> tris = _X.get_faces();

        std::vector<Matrix2d> E_inv(_X.n_faces());
        std::vector<real> detE(_X.n_faces());
        std::vector<Point2d> barycenters(_X.n_faces());
        std::vector<real> volK(_X.n_faces());

        std::vector<Matrix2d> Ec(_Xi_ref.n_faces());
        std::vector<Matrix2d> Ec_inv(_Xi_ref.n_faces());
        std::vector<real> detEc(_Xi_ref.n_faces());

        std::vector<Matrix2d> J(_X.n_faces());
        std::vector<real> detJ(_X.n_faces());

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

        std::vector<Matrix2d> GJ(_X.n_faces());
        std::vector<real> GdetJ(_X.n_faces());

        if(_Func == Functional::HUANG)
        {
            for(unsigned i = 0; i < _X.n_faces(); ++i)
            {
                real G_, GdetJ_;
                Matrix2d GJ_;
                functional_Huang(J[i], detJ[i], Mk[i], G_, GJ_, GdetJ_);
                GJ[i] = GJ_;
                GdetJ[i] = GdetJ_ * detJ[i];
            }
        }
        else if(_Func == Functional::WINSLOW)
        {
            for(unsigned i = 0; i < _X.n_faces(); ++i)
            {
                real G_, GdetJ_;
                Matrix2d GJ_;
                functional_Winslow(J[i], Mk[i], G_, GJ_, GdetJ_);
                GJ[i] = GJ_;
                GdetJ[i] = GdetJ_ * detJ[i];
            }
        }

        for(unsigned i = 0; i < _X.n_faces(); ++i)
        {
            Matrix2d I_xi_K = (E_inv[i] * GJ[i] + Ec_inv[i] * GdetJ[i]) * volK[i];
            const std::array<unsigned,3>& tri = tris[i];
            dxidt[2*tri[1]] += I_xi_K.a00;
            dxidt[2*tri[2]] += I_xi_K.a10;
            dxidt[2*tri[0]] -= I_xi_K.a00;
            dxidt[2*tri[0]] -= I_xi_K.a10;
            dxidt[2*tri[1]+1] += I_xi_K.a01;
            dxidt[2*tri[2]+1] += I_xi_K.a11;
            dxidt[2*tri[0]+1] -= I_xi_K.a01;
            dxidt[2*tri[0]+1] -= I_xi_K.a11;
        }

        std::vector<real> b_factor(_X.n_vertices());

        if(_Func == Functional::HUANG)
        {
            real p = 1.5;
            for(unsigned i = 0; i < _X.n_vertices(); ++i)
            {
                b_factor[i] = std::pow(_X.get_metric(i).det(), 0.5*(p-1));
            }
        }
        else if(_Func == Functional::WINSLOW)
        {
            for(unsigned i = 0; i < _X.n_vertices(); ++i)
            {
                b_factor[i] = std::pow(_X.get_metric(i).det(), 0.5);
            }
        }

        for(unsigned i = 0; i < b_factor.size(); ++i)
        {
            dxidt[2*i] *= (-1 / _tau * b_factor[i]);
            dxidt[2*i+1] *= (-1 / _tau * b_factor[i]);
        }

        for(unsigned i = 0; i < _X.n_vertices(); ++i)
        {
            if(_X.is_boundary(i))
            {
                dxidt[2*i] = 0;
                dxidt[2*i+1] = 0;
            }
        }
    }

    std::vector<real> reshape(const std::vector<Point2d>& points)
    {
        std::vector<real> result;
        result.reserve(points.size() * 2);
        for(const auto& p : points)
        {
            result.push_back(p.x());
            result.push_back(p.y());
        }
        return result;
    }

    std::vector<Point2d> reshape(const std::vector<real>& point_list)
    {
        std::vector<Point2d> result;
        result.reserve(point_list.size() / 2);
        for(size_t i = 0; i < point_list.size(); i += 2)
        {
            result.emplace_back(point_list[i], point_list[i + 1]);
        }
        return result;
    }
}