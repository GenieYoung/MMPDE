#include "metric.h"

namespace MMPDE
{
    std::vector<real> calc_vertex_value(const Trimesh2d& mesh, const std::function<real(const Point2d&)>& func)
    {
        std::vector<real> values(mesh.n_vertices());
        for(auto vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit)
        {
            values[vit->idx()] = func(mesh.point(*vit));
        }
        return values;
    }

    std::vector<Matrix2d> calc_vertex_metric(const Trimesh2d& mesh, const std::function<Matrix2d(const Point2d&)>& func)
    {
        std::vector<Matrix2d> metrics(mesh.n_vertices());
        for(auto vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit)
        {
            metrics[vit->idx()] = func(mesh.point(*vit));
        }
        return metrics;
    }

    std::vector<Matrix2d> calc_vertex_metric(const Trimesh2d& mesh, MetricType metric_type)
    {
        if(metric_type == MetricType::IDENTITY)
        {
            return calc_vertex_metric(mesh, [](const Point2d& p)->Matrix2d{return Matrix2d(1, 0, 0, 1);});
        }
        else 
        {
            return {};
        }
    }

    std::vector<std::pair<real, real>> grad_recovery(const Trimesh2d& mesh, const std::vector<real>& values)
    {
        std::vector<Matrix2d> face_edge_matrix = calc_face_edge_matrix(mesh);
        std::vector<Matrix2d> E_inv(mesh.n_faces());
        std::vector<real> det_E(mesh.n_faces());
        for(unsigned i = 0; i < mesh.n_faces(); ++i)
        { 
            E_inv[i] = face_edge_matrix[i].inverse();
            det_E[i] = std::abs(face_edge_matrix[i].det()) / 2.0;
        }
        
        std::vector<std::array<real, 6>> phi(mesh.n_faces());
        for(unsigned i = 0; i < mesh.n_faces(); ++i)
        {
            phi[i][0] -= E_inv[i].a00;
            phi[i][0] -= E_inv[i].a01;
            phi[i][1] -= E_inv[i].a10;
            phi[i][1] -= E_inv[i].a11;
            phi[i][2] = E_inv[i].a00;
            phi[i][3] = E_inv[i].a10;
            phi[i][4] = E_inv[i].a01;
            phi[i][5] = E_inv[i].a11;
        }

        std::vector<std::array<unsigned, 3>> tris = get_faces(mesh);

        std::vector<std::pair<real, real>> tri_grads(mesh.n_faces(), {0, 0});
        for(unsigned i = 0; i < mesh.n_faces(); ++i)
        {
            tri_grads[i].first += (phi[i][0] * values[tris[i][0]] + 
                                   phi[i][2] * values[tris[i][1]] +
                                   phi[i][4] * values[tris[i][2]]);
            tri_grads[i].second += (phi[i][1] * values[tris[i][0]] + 
                                    phi[i][3] * values[tris[i][1]] +
                                    phi[i][5] * values[tris[i][2]]);
        }
        for(unsigned i = 0; i < mesh.n_faces(); ++i)
        {
            tri_grads[i].first *= det_E[i];
            tri_grads[i].second *= det_E[i];
        }

        std::vector<std::pair<real, real>> v_grads(mesh.n_vertices(), {0,0});
        std::vector<real> v_volumes(mesh.n_vertices(), 0.0);
        for(unsigned i = 0; i < tris.size(); ++i)
        {
            const std::array<unsigned, 3>& tri = tris[i];
            v_grads[tri[0]].first += tri_grads[i].first;
            v_grads[tri[1]].first += tri_grads[i].first;
            v_grads[tri[2]].first += tri_grads[i].first;
            v_grads[tri[0]].second += tri_grads[i].second;
            v_grads[tri[1]].second += tri_grads[i].second;
            v_grads[tri[2]].second += tri_grads[i].second;

            v_volumes[tri[0]] += det_E[i];
            v_volumes[tri[1]] += det_E[i];
            v_volumes[tri[2]] += det_E[i];
        }

        for(unsigned i = 0; i < mesh.n_vertices(); ++i)
        {
            v_grads[i].first /= v_volumes[i];
            v_grads[i].second /= v_volumes[i];
        }

        return v_grads;
    }

    std::vector<Matrix2d> calc_vertex_metric(const Trimesh2d& mesh, const std::vector<real>& values)
    {
        std::vector<std::pair<real, real>> v_grads = grad_recovery(mesh, values);

        std::vector<real> v_grads_x(mesh.n_vertices());
        std::vector<real> v_grads_y(mesh.n_vertices());

        for(unsigned i = 0; i < v_grads.size(); ++i)
        {
            v_grads_x[i] = v_grads[i].first;
            v_grads_y[i] = v_grads[i].second;
        }

        std::vector<std::pair<real, real>> v_hessian_x = grad_recovery(mesh, v_grads_x);
        std::vector<std::pair<real, real>> v_hessian_y = grad_recovery(mesh, v_grads_y);

        std::vector<Matrix2d> v_hessian(mesh.n_vertices());
        for(unsigned i = 0; i < mesh.n_vertices(); ++i)
        {
            if(!mesh.is_boundary(mesh.vertex_handle(i)))
            {
                Matrix2d H(v_hessian_x[i].first, v_hessian_x[i].second, v_hessian_y[i].first, v_hessian_y[i].second);
                v_hessian[i] = (H + H.transpose()) * 0.5;
            }
        }

        for(unsigned i = 0; i < mesh.n_vertices(); ++i)
        {
            v_hessian[i] = (Matrix2d(1,0,0,1) + v_hessian[i]).sqrt();
        }

        return v_hessian;
    }
}