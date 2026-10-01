#include "ConvergenceMonitor.h"

#include <algorithm>
#include <cmath>

namespace
{
    // f(face, area) for every face, the area from inputVertexPositions (no geometry quantity required)
    template <typename F>
    void for_each_face_area(ManifoldSurfaceMesh &mesh, const VertexPositionGeometry &geometry, F &&f)
    {
        for (Face face : mesh.faces())
        {
            Halfedge he = face.halfedge();
            const Vector3 &a = geometry.inputVertexPositions[he.vertex()];
            const Vector3 &b = geometry.inputVertexPositions[he.next().vertex()];
            const Vector3 &c = geometry.inputVertexPositions[he.next().next().vertex()];
            f(face, 0.5 * norm(cross(b - a, c - a)));
        }
    }
} // namespace

double ConvergenceMonitor::window_dE_rel(double E) const
{
    return std::fabs(sum_dE) / std::max(std::fabs(E), p.E_floor);
}

double force_density_rms(ManifoldSurfaceMesh &mesh, const VertexPositionGeometry &geometry,
                         const VertexData<Vector3> &force, const VertexData<Vector3> *normals)
{
    VertexData<double> dual_area(mesh, 0.0);
    for_each_face_area(mesh, geometry, [&](Face face, double area)
                       {
        for (Vertex v : face.adjacentVertices())
            dual_area[v] += area / 3.0; });

    double sum = 0.0;
    size_t n = 0;
    for (Vertex v : mesh.vertices())
    {
        double f = normals ? dot(force[v], (*normals)[v]) : norm(force[v]);
        double p = f / dual_area[v];
        sum += p * p;
        n++;
    }
    return n ? std::sqrt(sum / n) : 0.0;
}

double force_density_scale(ManifoldSurfaceMesh &mesh, const VertexPositionGeometry &geometry, double bending_modulus)
{
    double A = 0.0;
    for_each_face_area(mesh, geometry, [&](Face, double area)
                       { A += area; });
    double R = std::sqrt(A / (4.0 * M_PI));
    double KB = bending_modulus > 0.0 ? bending_modulus : 1.0;
    return KB / (R * R * R);
}

double bending_modulus(const std::vector<std::string> &energies, const std::vector<std::vector<double>> &constants)
{
    for (size_t i = 0; i < energies.size(); i++)
    {
        const std::string &e = energies[i];
        if ((e == "Bending" || e == "Bending_tan" || e == "H1_Bending" || e == "H2_Bending") && !constants[i].empty() &&
            constants[i][0] >= 1e-5)
            return constants[i][0];
    }
    return 0.0;
}
