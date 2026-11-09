// Checks for the isotropic remesher (RemeshOptions::isotropic):
//   1. on a stretched, jittered sphere it brings the edges into [size_min, size_max]
//      and leaves (almost) no non-Delaunay edge
//   2. it does not use the curvature: an open flat mesh ends in the same range
//      (pass it as the first argument; boundary vertices are not collapsed)
//   3. nothing is NaN, the enclosed volume of a closed mesh is kept
//   4. a second remesh finds little to do
//
// Usage: Test_isotropic_remesh [mesh.obj] [size_min] [size_max]
//        (defaults ../../../input/Simplesphere.obj 0.05 0.2)

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/remeshing.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "RemeshMonitor.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace
{
    struct Stats
    {
        size_t n_vertices, n_edges;
        double l_min, l_max, l_mean;
        double frac_out;         // outside [size_min, size_max]
        double frac_non_delaunay; // 2(cos a + cos b) < 0, the strict test
        double frac_flip;         // what the remesher's own (-0.1) test would flip
        double min_angle, volume, area;
        size_t valence[3]; // < 5, 5-7, > 7
        bool nan;
    };

    Stats measure(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, const RemeshOptions &options)
    {
        geom.refreshQuantities();
        const VertexData<Vector3> &P = geom.inputVertexPositions;
        Stats s{};
        s.n_vertices = mesh.nVertices();
        s.n_edges = mesh.nEdges();
        s.l_min = 1e300;
        size_t out = 0, non_delaunay = 0;
        for (Edge e : mesh.edges())
        {
            double L = norm(P[e.firstVertex()] - P[e.secondVertex()]);
            s.l_min = std::min(s.l_min, L);
            s.l_max = std::max(s.l_max, L);
            s.l_mean += L / s.n_edges;
            out += L < options.min_absolute_length || L > options.max_absolute_length;
            if (!e.isBoundary())
            {
                Halfedge he = e.halfedge();
                Vector3 p0 = P[he.vertex()], p1 = P[he.next().vertex()];
                Vector3 p2 = P[he.next().next().vertex()], p3 = P[he.twin().next().next().vertex()];
                double la = norm(p0 - p1), lb = norm(p1 - p2), lc = norm(p2 - p0), ld = norm(p3 - p1),
                       le = norm(p3 - p0);
                double cc = (lb * lb + lc * lc - la * la) / (lb * lc) + (ld * ld + le * le - la * la) / (ld * le);
                non_delaunay += cc < 0.0;
            }
        }
        s.frac_out = double(out) / s.n_edges;
        s.frac_non_delaunay = double(non_delaunay) / s.n_edges;
        MeshQuality quality = measure_mesh_quality(mesh, geom, options);
        s.frac_flip = quality.fraction(quality.n_flip);
        s.min_angle = 1e300;
        for (Face f : mesh.faces())
        {
            Halfedge he = f.halfedge();
            Vector3 a = P[he.vertex()], b = P[he.next().vertex()], c = P[he.next().next().vertex()];
            s.area += 0.5 * norm(cross(b - a, c - a));
            s.volume += dot(a, cross(b, c)) / 6.0;
            double angles[3] = {angle(b - a, c - a), angle(a - b, c - b), angle(a - c, b - c)};
            for (double t : angles)
                s.min_angle = std::min(s.min_angle, t);
        }
        for (Vertex v : mesh.vertices())
        {
            size_t d = v.degree();
            s.valence[d < 5 ? 0 : (d <= 7 ? 1 : 2)]++;
            if (!(std::isfinite(P[v].x) && std::isfinite(P[v].y) && std::isfinite(P[v].z)))
                s.nan = true;
        }
        return s;
    }

    void print(const char *label, const Stats &s)
    {
        std::printf("%-20s V %5zu E %5zu  L [%.4f %.4f] mean %.4f  out %.4f  non-Delaunay %.4f  (remesher flip %.4f)\n"
                    "%-20s min angle %5.1f deg  valence <5 %zu  5-7 %zu  >7 %zu  area %.4f  volume %.4f\n",
                    label, s.n_vertices, s.n_edges, s.l_min, s.l_max, s.l_mean, s.frac_out, s.frac_non_delaunay,
                    s.frac_flip, "", s.min_angle * 180.0 / M_PI, s.valence[0], s.valence[1], s.valence[2], s.area,
                    s.volume);
    }

    void stretch_and_jitter(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geometry)
    {
        unsigned long long state = 12345;
        auto uniform = [&state]()
        {
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            return double(state >> 11) / double(1ULL << 53) - 0.5;
        };
        for (Vertex v : mesh.vertices())
        {
            Vector3 &p = geometry.inputVertexPositions[v];
            p *= 3.0;
            p.x *= 1.6;
            p += 0.02 * Vector3{uniform(), uniform(), uniform()};
        }
        geometry.refreshQuantities();
    }
} // namespace

int main(int argc, char **argv)
{
    std::string path = argc > 1 ? argv[1] : "../../../input/Simplesphere.obj";
    RemeshOptions options;
    options.isotropic = true;
    options.min_absolute_length = argc > 2 ? std::stod(argv[2]) : 0.05;
    options.max_absolute_length = argc > 3 ? std::stod(argv[3]) : 0.2;
    options.maxIterations = 5;
    options.aspect_min = 0.2;

    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(path);
    stretch_and_jitter(*mesh, *geometry);

    bool ok = true;
    Stats before = measure(*mesh, *geometry, options);
    print("stretched + jittered", before);

    int ops = remesh(*mesh, *geometry, options); // dispatches to isotropic_remesh
    Stats after = measure(*mesh, *geometry, options);
    print("after remesh", after);
    std::printf("  operations %d\n", ops);

    int ops2 = remesh(*mesh, *geometry, options);
    Stats again = measure(*mesh, *geometry, options);
    print("after second remesh", again);
    std::printf("  operations %d\n", ops2);

    // The same input through the curvature remesher, for comparison
    {
        std::unique_ptr<ManifoldSurfaceMesh> mesh_c;
        std::unique_ptr<VertexPositionGeometry> geometry_c;
        std::tie(mesh_c, geometry_c) = readManifoldSurfaceMesh(path);
        stretch_and_jitter(*mesh_c, *geometry_c);
        RemeshOptions curv = options;
        curv.isotropic = false;
        curv.refine_angle = 0.6;
        curv.maxIterations = 1;
        remesh(*mesh_c, *geometry_c, curv);
        print("curvature remesher", measure(*mesh_c, *geometry_c, options));
    }

    if (after.nan)
        ok = false, std::printf("FAIL: NaN position\n");
    if (!(after.frac_out < 0.02))
        ok = false, std::printf("FAIL: more than 2%% of the edges are outside [size_min, size_max]\n");
    if (!(after.frac_non_delaunay < 0.02))
        ok = false, std::printf("FAIL: more than 2%% of the edges are not Delaunay\n");
    // The smoothing after the flips moves the vertices, so a few edges can be flippable again
    if (!(after.frac_flip < 0.005))
        ok = false, std::printf("FAIL: more than 0.5%% of the edges are flippable for the remesher's flip test\n");
    if (!mesh->hasBoundary() && !(std::fabs(after.volume - before.volume) < 0.05 * std::fabs(before.volume)))
        ok = false, std::printf("FAIL: the volume changed by more than 5%%\n");
    if (!(after.min_angle > 10.0 * M_PI / 180.0))
        ok = false, std::printf("FAIL: a triangle angle below 10 degrees\n");
    if (!(ops2 < 0.5 * ops))
        ok = false, std::printf("FAIL: the second remesh still does half the operations of the first\n");

    std::printf(ok ? "PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
