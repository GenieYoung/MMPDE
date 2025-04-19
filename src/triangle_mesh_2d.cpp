// #include <triangle_mesh_2d.h>

// namespace MMPDE
// {
//     void build_triangle_mesh(TriangleMesh2d& mesh, const std::vector<real>& rows, const std::vector<real>& cols)
//     {
//         mesh.clear();

//         assert(rows.size() > 1 && cols.size() > 1);

//         std::vector<TriangleMesh2d::VertexHandle> vhs;
//         vhs.reserve(rows.size() * cols.size());
//         for(unsigned i = 0; i < rows.size(); ++i)
//         {
//             for(unsigned j = 0; j < cols.size(); ++j)
//             {
//                 vhs.push_back(mesh.add_vertex(CMTL::geo2d::Point<real>(rows[i], cols[j])));
//             }
//         }

//         for(unsigned i = 0; i < rows.size() - 1; ++i)
//         {
//             for(unsigned j = 0; j < cols.size() - 1; ++j)
//             {
//                 TriangleMesh2d::VertexHandle v0 = vhs[i * cols.size() + j];
//                 TriangleMesh2d::VertexHandle v1 = vhs[i * cols.size() + j + 1];
//                 TriangleMesh2d::VertexHandle v2 = vhs[(i + 1) * cols.size() + j + 1];
//                 TriangleMesh2d::VertexHandle v3 = vhs[(i + 1) * cols.size() + j];
//                 mesh.add_face(v0, v1, v2);
//                 mesh.add_face(v0, v2, v3);
//             }
//         }
//     }

//     void write_obj(const TriangleMesh2d& mesh, const std::string& file)
//     {
//         CMTL::io::write_obj(mesh, file);
//     }
// }