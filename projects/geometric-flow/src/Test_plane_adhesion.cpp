// Checks for Plane_Adhesion (E = -W * membrane area in contact with a plane):
//   1. membrane and plane forces against finite differences, on a sphere
//      touching a tilted plane
//   2. a flat square lying on a tilted plane gives exactly -W * area
//   3. moving the plane further than delta away gives E = 0
//   4. Frenkel, LJ and Morse profiles (the planar version of the sphere bead): forces against
//      finite differences, a flat patch at the minimum height gives -W * area, repulsion
//      for a patch closer to the plane, zero beyond the cutoff, C1 weight below the floor
//   5. back_repulsion R: faces that look away from the plane are repelled (E = +R W g |a|,
//      g = 1 - w below the radius): forces against finite differences on a flat vesicle
//      over the plane (front and back faces, edge-on ones included), flat patch value, zero
//      from the radius on, force always away from the plane, none for R = 0
//   6. the same with the plain cosine layer: the layer with the sign flipped on the back faces
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

        void attach(double W, double delta, Vector3 normal, Vector3 point, double sigma = 1.0,
                    bead_geometry::ShellKind kind = bead_geometry::ShellKind::Cosine, double alpha = 6.0, double back = 0.0)
        {
            plane.reset(new Plane_Adhesion(mesh.get(), geometry.get(), {W, sigma, delta, normal.x, normal.y, normal.z, back}));
            plane->shellKind = kind;
            plane->shellAlpha = alpha;
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

    // 4. Pair-potential profiles: the minimum is at height sigma above the plane, cutoff sigma (1 + width)
    {
        using bead_geometry::ShellKind;
        const double sigma = 0.3, width = 1.5, cutoff = sigma * (1.0 + width);
        struct Case
        {
            ShellKind kind;
            const char *name;
            double alpha;
        };
        const Case cases[] = {{ShellKind::Frenkel, "frenkel", 6.0}, {ShellKind::LJ, "lj", 6.0}, {ShellKind::Morse, "morse", 5.0}};

        for (const Case &c : cases)
        {
            // Sphere of radius 3 whose bottom is 0.55 sigma above the plane: part of it is in the
            // repulsive core, part in the well, so both regimes enter the finite differences
            Setup s;
            std::tie(s.mesh, s.geometry) = readManifoldSurfaceMesh(path);
            for (Vertex v : s.mesh->vertices())
                s.geometry->inputVertexPositions[v] *= 3.0;
            s.geometry->refreshQuantities();
            s.attach(W, cutoff, n, -(3.0 + 0.6 * sigma) * n + Vector3{0.05, 0.02, 0.0}, sigma, c.kind, c.alpha);

            double E0 = s.plane->Tot_Energy();
            VertexData<Vector3> F = s.plane->Gradient();
            Vector3 Fplane = s.bead.Total_force;
            int touching = 0;
            for (Face f : s.mesh->faces())
                touching += s.plane->Face_Energy(f) != 0.0;
            std::printf("%s: sphere above the plane, E = %.10g, %d faces in range\n", c.name, E0, touching);
            ok &= check(touching > 0 && std::isfinite(E0), "faces in range, finite energy");

            const double h = 1e-7;
            double worst = 0.0, fmax = 0.0;
            for (Vertex v : s.mesh->vertices())
                for (int k = 0; k < 3; k++)
                {
                    Vector3 &p = s.geometry->inputVertexPositions[v];
                    const double old = p[k];
                    p[k] = old + h;
                    double Ep = s.plane->Tot_Energy();
                    p[k] = old - h;
                    double Em = s.plane->Tot_Energy();
                    p[k] = old;
                    worst = std::max(worst, std::fabs(-(Ep - Em) / (2 * h) - F[v][k]));
                    fmax = std::max(fmax, std::fabs(F[v][k]));
                }
            std::printf("  membrane force: max |FD - analytic| %.3e (max |F| %.3e)\n", worst, fmax);
            ok &= check(worst < 1e-5 * std::max(1.0, fmax), "membrane force matches finite differences");

            double pworst = 0.0;
            for (int k = 0; k < 3; k++)
            {
                const double old = s.bead.Pos[k];
                s.bead.Pos[k] = old + h;
                double Ep = s.plane->Tot_Energy();
                s.bead.Pos[k] = old - h;
                double Em = s.plane->Tot_Energy();
                s.bead.Pos[k] = old;
                pworst = std::max(pworst, std::fabs(-(Ep - Em) / (2 * h) - Fplane[k]));
            }
            ok &= check(pworst < 1e-5 * std::max(1.0, norm(Fplane)), "plane force matches finite differences");
            ok &= check(norm(Fplane - dot(Fplane, n) * n) < 1e-10 * std::max(1.0, norm(Fplane)), "plane force is along the normal");
        }

        // Flat 2 x 1.5 rectangle facing the plane at height sigma, closer, and beyond the cutoff
        Vector3 t1 = unit(cross(n, Vector3{1.0, 0.0, 0.0}));
        Vector3 t2 = cross(n, t1);
        Vector3 P{0.4, -1.0, 2.0};
        const double Lx = 2.0, Ly = 1.5;
        std::vector<Vector3> corners = {P, P + Lx * t1, P + Lx * t1 + Ly * t2, P + Ly * t2};
        std::vector<std::vector<size_t>> faces = {{0, 2, 1}, {0, 3, 2}};
        for (const Case &c : cases)
        {
            Setup s;
            s.mesh.reset(new ManifoldSurfaceMesh(faces));
            s.geometry.reset(new VertexPositionGeometry(*s.mesh));
            for (Vertex v : s.mesh->vertices())
                s.geometry->inputVertexPositions[v] = corners[v.getIndex()];
            s.geometry->refreshQuantities();
            s.attach(W, cutoff, n, P - sigma * n, sigma, c.kind, c.alpha);

            // the patch sits at h = sigma above the plane: the minimum of the well, w = 1
            double E = s.plane->Tot_Energy();
            std::printf("%s: flat patch at the minimum, E = %.15g, -W * area = %.15g\n", c.name, E, -W * Lx * Ly);
            ok &= check(std::fabs(E + W * Lx * Ly) < 1e-10, "flat patch at h = sigma gives -W * area");

            s.bead.Pos = P - 0.6 * sigma * n;
            ok &= check(s.plane->Tot_Energy() > 0.0, "repulsive core: positive energy below the minimum");
            s.bead.Pos = P + 0.2 * sigma * n; // the plane is above the patch (behind it): below the floor
            double Eb = s.plane->Tot_Energy();
            ok &= check(std::isfinite(Eb) && Eb > 0.0, "finite repulsive energy even behind the plane");
            s.bead.Pos = P - 1.01 * cutoff * n;
            ok &= check(s.plane->Tot_Energy() == 0.0, "no energy beyond the cutoff");
        }

        // C1 at the floor: w and dw/dh just above and below
        bead_geometry::ShellProfile profile;
        profile.kind = ShellKind::Frenkel;
        profile.width = width;
        const double hf = bead_geometry::PLANE_PROFILE_FLOOR * sigma, eps = 1e-9;
        double dA, dB;
        const double wA = bead_geometry::planeProfileWeight(hf + eps, sigma, dA, profile);
        const double wB = bead_geometry::planeProfileWeight(hf - eps, sigma, dB, profile);
        ok &= check(std::fabs(wA - wB) < 1e-3 * std::fabs(wA) && std::fabs(dA - dB) < 1e-6 * std::fabs(dA), "weight is C1 at the floor");
    }

    // 5. Repulsion of the faces that look away from the plane (back_repulsion)
    {
        using bead_geometry::ShellKind;
        const double sigma = 1.0, width = 1.5, cutoff = sigma * (1.0 + width), R = 0.7;
        struct Case
        {
            ShellKind kind;
            const char *name;
            double alpha;
        };
        const Case cases[] = {{ShellKind::Frenkel, "frenkel", 6.0}, {ShellKind::LJ, "lj", 6.0}, {ShellKind::Morse, "morse", 5.0}};

        for (const Case &c : cases)
        {
            // A flat vesicle (a sphere of radius 3 squeezed to a thickness of 0.3) over a slightly tilted plane: its
            // lower sheet looks at the plane (front, in the core of the well), its upper sheet away (back, in the
            // repulsive range), and the rim is edge-on in between
            Setup s;
            std::tie(s.mesh, s.geometry) = readManifoldSurfaceMesh(path);
            for (Vertex v : s.mesh->vertices())
            {
                Vector3 &p = s.geometry->inputVertexPositions[v];
                p = Vector3{3.0 * p.x, 3.0 * p.y, 0.15 * p.z};
            }
            s.geometry->refreshQuantities();
            const Vector3 nt = unit(Vector3{0.05, -0.03, 1.0});
            s.attach(W, cutoff, nt, -(0.15 + 0.55 * sigma) * nt, sigma, c.kind, c.alpha, R);

            int front = 0, back = 0;
            for (Face f : s.mesh->faces())
            {
                Halfedge he = f.halfedge();
                const Vector3 p0 = s.geometry->inputVertexPositions[he.vertex()];
                const Vector3 p1 = s.geometry->inputVertexPositions[he.next().vertex()];
                const Vector3 p2 = s.geometry->inputVertexPositions[he.next().next().vertex()];
                const double e = s.plane->Face_Energy(f);
                const double a = -0.5 * dot(cross(p1 - p0, p2 - p0), n);
                if (e != 0.0 && a > 0.0)
                    front++;
                if (e != 0.0 && a < 0.0)
                    back++;
            }
            double E0 = s.plane->Tot_Energy();
            VertexData<Vector3> F = s.plane->Gradient();
            Vector3 Fplane = s.bead.Total_force;
            std::printf("%s: flat vesicle over the plane, E = %.10g, %d front and %d back faces in range\n", c.name, E0, front, back);
            ok &= check(front > 0 && back > 0 && std::isfinite(E0), "both orientations in range, finite energy");

            const double h = 1e-7;
            double worst = 0.0, fmax = 0.0;
            for (Vertex v : s.mesh->vertices())
                for (int k = 0; k < 3; k++)
                {
                    Vector3 &p = s.geometry->inputVertexPositions[v];
                    const double old = p[k];
                    p[k] = old + h;
                    double Ep = s.plane->Tot_Energy();
                    p[k] = old - h;
                    double Em = s.plane->Tot_Energy();
                    p[k] = old;
                    worst = std::max(worst, std::fabs(-(Ep - Em) / (2 * h) - F[v][k]));
                    fmax = std::max(fmax, std::fabs(F[v][k]));
                }
            std::printf("  membrane force: max |FD - analytic| %.3e (max |F| %.3e)\n", worst, fmax);
            ok &= check(worst < 1e-5 * std::max(1.0, fmax), "membrane force matches finite differences");

            double pworst = 0.0;
            for (int k = 0; k < 3; k++)
            {
                const double old = s.bead.Pos[k];
                s.bead.Pos[k] = old + h;
                double Ep = s.plane->Tot_Energy();
                s.bead.Pos[k] = old - h;
                double Em = s.plane->Tot_Energy();
                s.bead.Pos[k] = old;
                pworst = std::max(pworst, std::fabs(-(Ep - Em) / (2 * h) - Fplane[k]));
            }
            ok &= check(pworst < 1e-5 * std::max(1.0, norm(Fplane)), "plane force matches finite differences");
        }

        // Flat 2 x 1.5 rectangle that looks away from the plane (the faces are reversed), at heights around the radius
        Vector3 t1 = unit(cross(n, Vector3{1.0, 0.0, 0.0}));
        Vector3 t2 = cross(n, t1);
        Vector3 P{0.4, -1.0, 2.0};
        const double Lx = 2.0, Ly = 1.5, area = Lx * Ly;
        std::vector<Vector3> corners = {P, P + Lx * t1, P + Lx * t1 + Ly * t2, P + Ly * t2};
        std::vector<std::vector<size_t>> reversed = {{0, 1, 2}, {0, 2, 3}}; // counter-clockwise seen from +n: the normal points along +n
        for (const Case &c : cases)
        {
            Setup s;
            s.mesh.reset(new ManifoldSurfaceMesh(reversed));
            s.geometry.reset(new VertexPositionGeometry(*s.mesh));
            for (Vertex v : s.mesh->vertices())
                s.geometry->inputVertexPositions[v] = corners[v.getIndex()];
            s.geometry->refreshQuantities();

            s.attach(W, cutoff, n, P - sigma * n, sigma, c.kind, c.alpha, R);
            ok &= check(std::fabs(s.plane->Tot_Energy()) < 1e-12, "back patch at h = radius: no energy (g = 0 there)");
            s.bead.Pos = P - 1.5 * sigma * n;
            ok &= check(s.plane->Tot_Energy() == 0.0, "back patch above the radius: no energy");

            // Below the radius the energy is R W g |a| > 0 and falls monotonically with the height, so the force pushes away from the plane
            s.bead.Pos = P - 0.6 * sigma * n;
            const double E = s.plane->Tot_Energy();
            bead_geometry::ShellProfile profile;
            profile.kind = c.kind;
            profile.width = width;
            profile.alpha = c.alpha;
            double dg;
            const double g = bead_geometry::planeProfileRepulsion(0.6 * sigma, sigma, dg, profile);
            std::printf("%s: back patch at 0.6 radius, E = %.10g, R W g area = %.10g\n", c.name, E, R * W * g * area);
            ok &= check(E > 0.0 && std::fabs(E - R * W * g * area) < 1e-9 * E, "back patch below the radius: E = R W g area");
            VertexData<Vector3> Fv = s.plane->Gradient();
            double push = 0.0;
            for (Vertex v : s.mesh->vertices())
                push += dot(Fv[v], n);
            ok &= check(push > 0.0 && dot(s.bead.Total_force, n) < 0.0, "the patch is pushed along +n, away from the plane (the plane the other way)");
            double Eprev = 1e300;
            bool decreasing = true;
            for (double hh = 0.3 * sigma; hh < sigma; hh += 0.05 * sigma)
            {
                s.bead.Pos = P - hh * n;
                const double Eh = s.plane->Tot_Energy();
                decreasing &= Eh < Eprev;
                Eprev = Eh;
            }
            ok &= check(decreasing, "energy falls monotonically with the height up to the radius");

            // R = 0 switches it off
            s.attach(W, cutoff, n, P - 0.6 * sigma * n, sigma, c.kind, c.alpha, 0.0);
            ok &= check(s.plane->Tot_Energy() == 0.0, "back_repulsion = 0: back faces feel nothing");
        }

        // Tilting the patch towards edge-on takes the energy to 0 continuously
        {
            Setup s;
            double previous = -1.0;
            bool vanishes = true;
            for (double tilt : {0.5, 0.1, 0.01, 1e-4})
            {
                // the back patch rotated by a small angle from the plane around t1: its normal turns to be nearly along t2
                const double ang = 1.5707963267948966 - tilt, ca = std::cos(ang), sa = std::sin(ang);
                std::vector<Vector3> tilted;
                for (const Vector3 &q : corners)
                {
                    const Vector3 d = 0.1 * (q - P); // a small patch, so that all of it stays within the radius
                    const double dn = dot(d, n), d2 = dot(d, t2), d1 = dot(d, t1);
                    tilted.push_back(P + d1 * t1 + (ca * d2 - sa * dn) * t2 + (sa * d2 + ca * dn) * n);
                }
                s.mesh.reset(new ManifoldSurfaceMesh(reversed));
                s.geometry.reset(new VertexPositionGeometry(*s.mesh));
                for (Vertex v : s.mesh->vertices())
                    s.geometry->inputVertexPositions[v] = tilted[v.getIndex()];
                s.geometry->refreshQuantities();
                s.attach(W, cutoff, n, P - 0.3 * sigma * n, sigma, ShellKind::Frenkel, 6.0, R);
                const double E = std::fabs(s.plane->Tot_Energy());
                vanishes &= (previous < 0.0 || E < previous);
                previous = E;
            }
            std::printf("tilting a back patch towards edge-on: |E| ends at %.3e\n", previous);
            ok &= check(vanishes, "back energy decreases as the patch turns edge-on");
        }
    }

    // 6. back_repulsion with the plain cosine layer: the same layer, sign flipped, on the faces that look away
    {
        const double delta = 0.4, R = 0.7;

        // A flat vesicle straddling a slightly tilted plane: the lower sheet looks at it (adhesion), the upper one away
        // (repulsion), both within the layer
        Setup s;
        std::tie(s.mesh, s.geometry) = readManifoldSurfaceMesh(path);
        for (Vertex v : s.mesh->vertices())
        {
            Vector3 &p = s.geometry->inputVertexPositions[v];
            p = Vector3{3.0 * p.x, 3.0 * p.y, 0.15 * p.z};
        }
        s.geometry->refreshQuantities();
        const Vector3 nt = unit(Vector3{0.05, -0.03, 1.0});
        s.attach(W, delta, nt, Vector3{0.0, 0.0, 0.0}, 1.0, bead_geometry::ShellKind::Cosine, 6.0, R);

        int front = 0, back = 0;
        for (Face f : s.mesh->faces())
        {
            Halfedge he = f.halfedge();
            const Vector3 p0 = s.geometry->inputVertexPositions[he.vertex()];
            const Vector3 p1 = s.geometry->inputVertexPositions[he.next().vertex()];
            const Vector3 p2 = s.geometry->inputVertexPositions[he.next().next().vertex()];
            const double e = s.plane->Face_Energy(f);
            const double a = -0.5 * dot(cross(p1 - p0, p2 - p0), nt);
            front += e != 0.0 && a > 0.0;
            back += e != 0.0 && a < 0.0;
        }
        double E0 = s.plane->Tot_Energy();
        VertexData<Vector3> F = s.plane->Gradient();
        Vector3 Fplane = s.bead.Total_force;
        std::printf("cosine: flat vesicle across the plane, E = %.10g, %d front and %d back faces in range\n", E0, front, back);
        ok &= check(front > 0 && back > 0 && std::isfinite(E0), "both orientations in range, finite energy");

        const double h = 1e-7;
        double worst = 0.0, fmax = 0.0;
        for (Vertex v : s.mesh->vertices())
            for (int k = 0; k < 3; k++)
            {
                Vector3 &p = s.geometry->inputVertexPositions[v];
                const double old = p[k];
                p[k] = old + h;
                double Ep = s.plane->Tot_Energy();
                p[k] = old - h;
                double Em = s.plane->Tot_Energy();
                p[k] = old;
                worst = std::max(worst, std::fabs(-(Ep - Em) / (2 * h) - F[v][k]));
                fmax = std::max(fmax, std::fabs(F[v][k]));
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
            pworst = std::max(pworst, std::fabs(-(Ep - Em) / (2 * h) - Fplane[k]));
        }
        ok &= check(pworst < 1e-6 * std::max(1.0, norm(Fplane)), "plane force matches finite differences");

        // Flat rectangle that looks away from the plane: +R W area on the plane, zero beyond the layer, none for R = 0
        Vector3 t1 = unit(cross(n, Vector3{1.0, 0.0, 0.0}));
        Vector3 t2 = cross(n, t1);
        Vector3 P{0.4, -1.0, 2.0};
        const double Lx = 2.0, Ly = 1.5;
        std::vector<Vector3> corners = {P, P + Lx * t1, P + Lx * t1 + Ly * t2, P + Ly * t2};
        std::vector<std::vector<size_t>> reversed = {{0, 1, 2}, {0, 2, 3}};
        Setup r;
        r.mesh.reset(new ManifoldSurfaceMesh(reversed));
        r.geometry.reset(new VertexPositionGeometry(*r.mesh));
        for (Vertex v : r.mesh->vertices())
            r.geometry->inputVertexPositions[v] = corners[v.getIndex()];
        r.geometry->refreshQuantities();
        r.attach(W, delta, n, P, 1.0, bead_geometry::ShellKind::Cosine, 6.0, R);
        const double E = r.plane->Tot_Energy();
        std::printf("cosine: back patch on the plane, E = %.15g, R W area = %.15g\n", E, R * W * Lx * Ly);
        ok &= check(std::fabs(E - R * W * Lx * Ly) < 1e-12, "back patch on the plane: E = +R W area (front would give -W area)");
        r.bead.Pos = P - 1.01 * delta * n;
        ok &= check(r.plane->Tot_Energy() == 0.0, "no energy beyond the layer");
        r.bead.Pos = P - 0.5 * delta * n; // the patch is 0.5 delta above the plane: pushed away from it
        VertexData<Vector3> Fv = r.plane->Gradient();
        double push = 0.0;
        for (Vertex v : r.mesh->vertices())
            push += dot(Fv[v], n);
        ok &= check(push > 0.0, "the patch above the plane is pushed along +n");
        r.attach(W, delta, n, P, 1.0, bead_geometry::ShellKind::Cosine, 6.0, 0.0);
        ok &= check(r.plane->Tot_Energy() == 0.0, "back_repulsion = 0: back faces feel nothing");
    }

    std::printf(ok ? "PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
