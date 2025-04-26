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
        const std::vector<Matrix2d>& M,
        real tau,
        Functional Func)
    {
        assert(X.n_faces() == Xi_ref.n_faces());

        MoveMeshRHS rhs(Xi_ref, X, M, tau, Func);

        std::vector<Point2d> old_xi_vertices = get_vertices(Xi_ref);
        std::vector<real> xi = reshape(old_xi_vertices);
        
        double dt = (tspan.second - tspan.first) / 20.0;
        std::vector<std::vector<real>> xi_vec;
        std::vector<double> times;
        boost::numeric::odeint::integrate(rhs, xi, tspan.first, tspan.second, dt, push_back_state_and_time(xi_vec, times));

        Interpolate intp(reshape(xi), get_faces(Xi_ref));
        std::vector<Point2d> new_xi_vertices = intp(get_vertices(X), old_xi_vertices);

        Trimesh2d new_mesh = X;
        update_trimesh(new_mesh, new_xi_vertices);

        return new_mesh;
    }

    MoveMeshRHS::MoveMeshRHS(const Trimesh2d& Xi_ref, 
                             const Trimesh2d& X,
                             const std::vector<Matrix2d>& M,
                             real tau, 
                             Functional Func)
        : _tau(tau)
    {
        _tris = get_faces(X);
        _boundary_ids = get_boundary_ids(X);
        _barycenters = calc_face_barycenter(X);

        _Mk.resize(X.n_faces());
        for(unsigned i = 0; i < X.n_faces(); ++i)
        {
            for(unsigned j = 0; j < 3; ++j)
            {
                _Mk[i] = _Mk[i] + M[_tris[i][j]];
            }
            _Mk[i] = _Mk[i] / 3.0;
        }

        _E_inv.resize(X.n_faces());
        _detE.resize(X.n_faces());
        _volK.resize(X.n_faces());
        _sigma = 0;
        std::vector<Matrix2d> face_edge_matrix = calc_face_edge_matrix(X);
        for(unsigned i = 0; i < X.n_faces(); ++i)
        { 
            const Matrix2d& E = face_edge_matrix[i];
            _E_inv[i] = E.inverse();
            _detE[i] = std::abs(E.det());
            _volK[i] = 0.5 * _detE[i];
            _sigma += _volK[i] * _Mk[i].det();
        }

        _Ec_inv.resize(Xi_ref.n_faces());
        _detEc.resize(Xi_ref.n_faces());
        _J.resize(X.n_faces());
        _detJ.resize(X.n_faces());
        _Ec = calc_face_edge_matrix(Xi_ref);
        for(unsigned i = 0; i < Xi_ref.n_faces(); ++i)
        {
            _Ec_inv[i] = _Ec[i].inverse();
            _detEc[i] = std::abs(_Ec[i].det());
            _J[i] = _Ec[i] * _E_inv[i];
            _detJ[i] = _detEc[i] / _detE[i];
        }

        _GJ.resize(X.n_faces());
        _GdetJ.resize(X.n_faces());
        if(Func == Functional::HUANG)
        {
            for(unsigned i = 0; i < X.n_faces(); ++i)
            {
                real G_, GdetJ_;
                Matrix2d GJ_;
                functional_Huang(_J[i], _detJ[i], _Mk[i], G_, GJ_, GdetJ_);
                _GJ[i] = GJ_;
                _GdetJ[i] = GdetJ_ * _detJ[i];
            }
        }
        else if(Func == Functional::WINSLOW)
        {
            for(unsigned i = 0; i < X.n_faces(); ++i)
            {
                real G_, GdetJ_;
                Matrix2d GJ_;
                functional_Winslow(_J[i], _Mk[i], G_, GJ_, GdetJ_);
                _GJ[i] = GJ_;
                _GdetJ[i] = GdetJ_ * _detJ[i];
            }
        }

        _b_factor.resize(X.n_vertices());
        if(Func == Functional::HUANG)
        {
            real p = 1.5;
            for(unsigned i = 0; i < X.n_vertices(); ++i)
            {
                _b_factor[i] = std::pow(M[i].det(), 0.5*(p-1));
            }
        }
        else if(Func == Functional::WINSLOW)
        {
            for(unsigned i = 0; i < X.n_vertices(); ++i)
            {
                _b_factor[i] = std::pow(M[i].det(), 0.5);
            }
        }
    }

    void MoveMeshRHS::operator()(const std::vector<real>& xi,
                                 std::vector<real>& dxidt,
                                 const double t) const
    {
        dxidt.resize(xi.size(), 0);

        for(unsigned i = 0; i < _tris.size(); ++i)
        {
            Matrix2d I_xi_K = ((_E_inv[i] * _GJ[i] + _Ec_inv[i] * _GdetJ[i]) * _volK[i]).transpose();
            const std::array<unsigned,3>& tri = _tris[i];
            dxidt[2*tri[1]] += I_xi_K.a00;
            dxidt[2*tri[2]] += I_xi_K.a10;
            dxidt[2*tri[0]] -= I_xi_K.a00;
            dxidt[2*tri[0]] -= I_xi_K.a10;
            dxidt[2*tri[1]+1] += I_xi_K.a01;
            dxidt[2*tri[2]+1] += I_xi_K.a11;
            dxidt[2*tri[0]+1] -= I_xi_K.a01;
            dxidt[2*tri[0]+1] -= I_xi_K.a11;
        }

        for(unsigned i = 0; i < _b_factor.size(); ++i)
        {
            dxidt[2*i] *= (-1 / _tau * _b_factor[i]);
            dxidt[2*i+1] *= (-1 / _tau * _b_factor[i]);
        }

        for(unsigned i = 0; i < _boundary_ids.size(); ++i)
        {
            dxidt[2*_boundary_ids[i]] = 0;
            dxidt[2*_boundary_ids[i]+1] = 0;
        }

        int s;
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