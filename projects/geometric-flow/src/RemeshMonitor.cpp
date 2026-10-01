#include "RemeshMonitor.h"

#include <algorithm>
#include <cmath>

namespace
{
    // isDelaunay_improv in remeshing.cpp: 2(cos a + cos b) >= -0.1, with a and b
    // the angles opposite the edge
    bool is_delaunay(const VertexData<Vector3> &P, Edge e)
    {
        Halfedge he = e.halfedge();
        Vector3 p0 = P[he.vertex()];
        Vector3 p1 = P[he.next().vertex()];
        Vector3 p2 = P[he.next().next().vertex()];
        Vector3 p3 = P[he.twin().next().next().vertex()];
        double la = norm(p0 - p1);
        double lb = norm(p1 - p2);
        double lc = norm(p2 - p0);
        double ld = norm(p3 - p1);
        double le = norm(p3 - p0);
        return (lb * lb + lc * lc - la * la) / (lb * lc) + (ld * ld + le * le - la * la) / (ld * le) >= -0.1;
    }

    // Metric part of shouldCollapse: collapsing e to its midpoint keeps every
    // edge of the link of e below 0.9 in the sizing metric
    bool metric_allows_collapse(const VertexData<Vector3> &P, const VertexData<double> &sizing, Edge e)
    {
        Vertex v1 = e.halfedge().vertex();
        Vertex v2 = e.halfedge().twin().vertex();
        Vector3 mid = 0.5 * (P[v1] + P[v2]);
        double new_sizing = std::max(sizing[v1], sizing[v2]);

        for (Vertex end : {v1, v2})
        {
            Vertex other = end == v1 ? v2 : v1;
            for (Halfedge he : end.outgoingHalfedges())
            {
                Halfedge link = he.next();
                if (link.tailVertex() == other || link.tipVertex() == other)
                    continue;
                Vertex b = link.tailVertex(), c = link.tipVertex();
                if (norm(P[c] - mid) * std::sqrt((new_sizing + sizing[c]) / 2.0) > 0.9)
                    return false;
                if (norm(P[b] - mid) * std::sqrt((new_sizing + sizing[b]) / 2.0) > 0.9)
                    return false;
                if (norm(P[b] - P[c]) * std::sqrt((sizing[b] + sizing[c]) / 2.0) > 0.9)
                    return false;
            }
        }
        return true;
    }
} // namespace

VertexData<double> remesh_vertex_sizing(ManifoldSurfaceMesh &mesh, const VertexPositionGeometry &geometry,
                                        const RemeshOptions &options)
{
    const VertexData<Vector3> &P = geometry.inputVertexPositions;

    EdgeData<double> length(mesh), dihedral(mesh);
    for (Edge e : mesh.edges())
    {
        length[e] = geometry.edgeLength(e);
        dihedral[e] = geometry.edgeDihedralAngle(e);
    }

    // Face sizing as in computeFaceSizing + remesh(): largest eigenvalue of
    // S^T S with S = (1/A) sum_edges -alpha L/2 n n^T (n the in-plane edge
    // normal), divided by refine_angle^2 and clamped. S is symmetric, so that
    // eigenvalue is the square of its spectral radius.
    const double inv_angle2 = 1.0 / (options.refine_angle * options.refine_angle);
    const double sizing_lo = 1.0 / (options.max_absolute_length * options.max_absolute_length);
    const double sizing_hi = 1.0 / (options.min_absolute_length * options.min_absolute_length);
    FaceData<double> area(mesh), face_sizing(mesh);
    for (Face f : mesh.faces())
    {
        Halfedge he = f.halfedge();
        Vector3 a = P[he.vertex()], b = P[he.next().vertex()], c = P[he.next().next().vertex()];
        Vector3 N = cross(b - a, c - a);
        double A = 0.5 * norm(N);
        area[f] = A;
        if (A <= 0.0)
        {
            face_sizing[f] = sizing_hi;
            continue;
        }
        N /= 2.0 * A;
        Vector3 t1 = unit(b - a), t2 = cross(N, t1);
        double s00 = 0.0, s01 = 0.0, s11 = 0.0;
        for (Halfedge h : f.adjacentHalfedges())
        {
            Vector3 n = unit(cross(N, P[h.tipVertex()] - P[h.tailVertex()]));
            double x = dot(n, t1), y = dot(n, t2);
            double w = -0.5 * dihedral[h.edge()] * length[h.edge()] / A;
            s00 += w * x * x;
            s01 += w * x * y;
            s11 += w * y * y;
        }
        double rho = std::fabs(0.5 * (s00 + s11)) + std::sqrt(0.25 * (s00 - s11) * (s00 - s11) + s01 * s01);
        face_sizing[f] = std::min(std::max(rho * rho * inv_angle2, sizing_lo), sizing_hi);
    }

    // computeVertexSizing: area weighted mean over the adjacent faces
    VertexData<double> sizing(mesh, 0.0);
    for (Vertex v : mesh.vertices())
    {
        double weighted = 0.0, total = 0.0;
        for (Face f : v.adjacentFaces())
        {
            weighted += area[f] * face_sizing[f];
            total += area[f];
        }
        sizing[v] = total > 0.0 ? weighted / total : sizing_hi;
    }
    return sizing;
}

MeshQuality measure_mesh_quality(ManifoldSurfaceMesh &mesh, const VertexPositionGeometry &geometry,
                                 const RemeshOptions &options)
{
    const VertexData<Vector3> &P = geometry.inputVertexPositions;
    VertexData<double> sizing = remesh_vertex_sizing(mesh, geometry, options);

    MeshQuality q;
    q.n_edges = mesh.nEdges();
    for (Edge e : mesh.edges())
    {
        Vertex v1 = e.halfedge().vertex();
        Vertex v2 = e.halfedge().twin().vertex();
        double L = geometry.edgeLength(e);

        bool is_long = L * std::sqrt((sizing[v1] + sizing[v2]) / 2.0) > 1.0 &&
                       L > 2.0 * options.min_absolute_length;

        bool interior = !e.isBoundary() && !v1.isBoundary() && !v2.isBoundary();
        bool is_short = interior && L < options.max_absolute_length &&
                        (v1.degree() == 3 || v1.degree() == 4 || v2.degree() == 3 || v2.degree() == 4 ||
                         metric_allows_collapse(P, sizing, e));

        bool is_flip = !e.isBoundary() && !is_delaunay(P, e);

        q.n_long += is_long;
        q.n_short += is_short;
        q.n_flip += is_flip;
        q.n_bad += is_long || is_short || is_flip;
    }
    return q;
}
