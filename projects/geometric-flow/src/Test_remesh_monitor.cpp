// Checks for RemeshMonitor (the "quality" adaptive remeshing mode):
//   1. remesh_vertex_sizing matches geometry-central's faceSizing/vertexSizing
//      as remesh() sets them up
//   2. measure_mesh_quality does not move the mesh
//   3. on a stretched, jittered sphere remesh() removes most of the defects
//
// Usage: Test_remesh_monitor [mesh.obj]   (default ../../../input/Simplesphere.obj)

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
    void print(const char *label, const MeshQuality &q)
    {
        std::printf("%-22s edges %6zu  bad %.4f  long %.4f  short %.4f  flip %.4f\n", label, q.n_edges,
                    q.bad_fraction(), q.fraction(q.n_long), q.fraction(q.n_short), q.fraction(q.n_flip));
    }
} // namespace

int main(int argc, char **argv)
{
    std::string path = argc > 1 ? argv[1] : "../../../input/Simplesphere.obj";
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(path);

    // Same remesher settings as regression config b
    RemeshOptions options;
    options.max_absolute_length = 0.45;
    options.min_absolute_length = 0.05;
    options.refine_angle = 0.6;
    options.aspect_min = 0.2;
    options.maxIterations = 1;

    // Stretch along x and jitter, so there is something to fix
    unsigned long long state = 12345;
    auto uniform = [&state]()
    {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return double(state >> 11) / double(1ULL << 53) - 0.5;
    };
    for (Vertex v : mesh->vertices())
    {
        Vector3 &p = geometry->inputVertexPositions[v];
        p *= 3.0;
        p.x *= 1.6;
        p += 0.02 * Vector3{uniform(), uniform(), uniform()};
    }
    geometry->refreshQuantities();

    bool ok = true;

    // 1. Sizing against geometry-central (on a copy, so the caches of the
    //    geometry used below are not touched)
    {
        std::unique_ptr<VertexPositionGeometry> reference = geometry->copy();
        reference->requireFaceSizing();
        for (Face f : mesh->faces())
            reference->faceSizing[f] = clamp(reference->faceSizing[f] / (options.refine_angle * options.refine_angle),
                                             1.0 / (options.max_absolute_length * options.max_absolute_length),
                                             1.0 / (options.min_absolute_length * options.min_absolute_length));
        reference->requireVertexSizing();
        VertexData<double> mine = remesh_vertex_sizing(*mesh, *geometry, options);
        double worst = 0.0;
        for (Vertex v : mesh->vertices())
            worst = std::max(worst, std::fabs(mine[v] - reference->vertexSizing[v]) / reference->vertexSizing[v]);
        std::printf("sizing vs geometry-central: max rel diff %.3e\n", worst);
        if (worst > 1e-9)
        {
            std::printf("FAIL: sizing differs\n");
            ok = false;
        }
    }

    // 2. Read only
    VertexData<Vector3> positions = geometry->inputVertexPositions;
    MeshQuality before = measure_mesh_quality(*mesh, *geometry, options);
    for (Vertex v : mesh->vertices())
        if (positions[v] != geometry->inputVertexPositions[v])
        {
            std::printf("FAIL: measure_mesh_quality moved a vertex\n");
            ok = false;
            break;
        }

    // 3. A remesh removes the defects; a second one right after only finds
    //    what the first left behind
    print("stretched + jittered", before);
    int ops = remesh(*mesh, *geometry, options);
    geometry->refreshQuantities();
    MeshQuality after = measure_mesh_quality(*mesh, *geometry, options);
    print("after remesh", after);
    std::printf("  operations %d\n", ops);
    int ops2 = remesh(*mesh, *geometry, options);
    geometry->refreshQuantities();
    MeshQuality again = measure_mesh_quality(*mesh, *geometry, options);
    print("after second remesh", again);
    std::printf("  operations %d\n", ops2);

    if (!(after.bad_fraction() < 0.5 * before.bad_fraction()))
    {
        std::printf("FAIL: remesh did not halve the bad edge fraction\n");
        ok = false;
    }

    std::printf(ok ? "PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
