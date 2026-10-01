// Checks for Plane_Adhesion (E = -W * membrane area in contact with a plane):
//   1. membrane and plane forces against finite differences, on a sphere
//      touching a tilted plane
//   2. a flat square lying on a tilted plane gives exactly -W * area
//   3. moving the plane further than delta away gives E = 0
//
// Usage: Test_plane_adhesion [mesh.obj]   (default ../../../input/Simplesphere.obj)

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "Beads.h"
#include "Interaction.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace
{
    struct Setup
    {
        std::unique_ptr<ManifoldSurfaceMesh> mesh;
        std::unique_ptr<VertexPositionGeometry> geometry;
        std::unique_ptr<Plane_Adhesion> plane;
        Bead bead;

        void attach(double W, double delta, Vector3 normal, Vector3 point)
        {
            plane.reset(new Plane_Adhesion(mesh.get(), geometry.get(), {W, 1.0, delta, normal.x, normal.y, normal.z}));
            bead.mesh = mesh.get();
            bead.geometry = geometry.get();
            bead.Pos = point;
            bead.sigma = 1.0;
            bead.strength = W;
            bead.Bead_I = plane.get();
            plane->Bead_1 = &bead;
            bead.CoverageForce = Vector3{0.0, 0.0, 0.0};
        }
    };

    bool check(bool ok, const char *what)
    {
        std::printf("%s: %s\n", ok ? "ok  " : "FAIL", what);
        return ok;
    }
} // namespace

int main(int argc, char **argv)
{
    std::string path = argc > 1 ? argv[1] : "../../../input/Simplesphere.obj";
    const double W = 2.0, delta = 0.4;
    const Vector3 n = unit(Vector3{0.3, -0.2, 1.0});
    bool ok = true;

    // 1. Sphere of radius 3 whose bottom touches a tilted plane
    {
        Setup s;
        std::tie(s.mesh, s.geometry) = readManifoldSurfaceMesh(path);
        for (Vertex v : s.mesh->vertices())
            s.geometry->inputVertexPositions[v] *= 3.0;
        s.geometry->refreshQuantities();
        s.attach(W, delta, n, -2.9 * n + Vector3{0.05, 0.02, 0.0});

        double E0 = s.plane->Tot_Energy();
        VertexData<Vector3> F = s.plane->Gradient();
        Vector3 Fplane = s.bead.Total_force;
        int touching = 0;
        for (Face f : s.mesh->faces())
            touching += s.plane->Face_Energy(f) != 0.0;
        std::printf("sphere on tilted plane: E = %.10g, %d faces in contact\n", E0, touching);
        ok &= check(touching > 0 && E0 < 0.0, "some faces in contact, negative energy");

        const double h = 1e-6;
        double worst = 0.0, fmax = 0.0;
        for (Vertex v : s.mesh->vertices())
        {
            for (int k = 0; k < 3; k++)
            {
                Vector3 &p = s.geometry->inputVertexPositions[v];
                const double old = p[k];
                p[k] = old + h;
                double Ep = s.plane->Tot_Energy();
                p[k] = old - h;
                double Em = s.plane->Tot_Energy();
                p[k] = old;
                double fd = -(Ep - Em) / (2 * h);
                worst = std::max(worst, std::fabs(fd - F[v][k]));
                fmax = std::max(fmax, std::fabs(F[v][k]));
            }
        }
        std::printf("  membrane force: max |FD - analytic| %.3e (max |F| %.3e)\n", worst, fmax);
        ok &= check(worst < 1e-6 * std::max(1.0, fmax), "membrane force matches finite differences");

        double pworst = 0.0;
        for (int k = 0; k < 3; k++)
        {
            const double old = s.bead.Pos[k];
            s.bead.Pos[k] = old + h;
            double Ep = s.plane->Tot_Energy();
            s.bead.Pos[k] = old - h;
            double Em = s.plane->Tot_Energy();
            s.bead.Pos[k] = old;
            double fd = -(Ep - Em) / (2 * h);
            std::printf("  plane force %d: analytic %.8g FD %.8g\n", k, Fplane[k], fd);
            pworst = std::max(pworst, std::fabs(fd - Fplane[k]));
        }
        ok &= check(pworst < 1e-5 * std::max(1.0, norm(Fplane)), "plane force matches finite differences");
        ok &= check(norm(Fplane - dot(Fplane, n) * n) < 1e-10 * std::max(1.0, norm(Fplane)),
                    "plane force is along the normal");
    }

    // 2./3. A flat 2 x 1.5 rectangle lying in the tilted plane, facing it
    {
        Vector3 t1 = unit(cross(n, Vector3{1.0, 0.0, 0.0}));
        Vector3 t2 = cross(n, t1);
        Vector3 P{0.4, -1.0, 2.0};
        const double Lx = 2.0, Ly = 1.5;
        std::vector<Vector3> corners = {P, P + Lx * t1, P + Lx * t1 + Ly * t2, P + Ly * t2};
        // Clockwise seen from +n, so the outward normal points into the plane
        std::vector<std::vector<size_t>> faces = {{0, 2, 1}, {0, 3, 2}};

        Setup s;
        s.mesh.reset(new ManifoldSurfaceMesh(faces));
        s.geometry.reset(new VertexPositionGeometry(*s.mesh));
        for (Vertex v : s.mesh->vertices())
            s.geometry->inputVertexPositions[v] = corners[v.getIndex()];
        s.geometry->refreshQuantities();

        s.attach(W, delta, n, P);
        double E = s.plane->Tot_Energy();
        std::printf("flat patch on the plane: E = %.15g, -W * area = %.15g\n", E, -W * Lx * Ly);
        ok &= check(std::fabs(E + W * Lx * Ly) < 1e-12, "flat patch gives -W * area");

        s.bead.Pos = P - 1.01 * delta * n;
        ok &= check(s.plane->Tot_Energy() == 0.0, "no energy with the plane further than delta");
    }

    std::printf(ok ? "PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
