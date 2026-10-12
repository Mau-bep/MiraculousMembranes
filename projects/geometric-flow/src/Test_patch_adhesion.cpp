// Checks for the angular patch of the Adhesion interaction ("patch_angle", "patch_width", "patch_axis"):
//   E = -W sigma^2 sum_f w(r_f) w_patch(t_f) |omega_f|,   t_f = dot(u_f, axis),
//   w_patch = 0 for t <= cos(theta_Wr), 1 for t >= cos(theta_Wr - width), a C1 cosine ramp in between
//
//   1. the patch weight: 0 below the edge, 1 above the ramp, continuous with zero slope at both ends of the ramp,
//      monotone, d w / d t against central differences
//   2. Adhesion::Add_Face_Force against central finite differences of Adhesion::Face_Energy, membrane vertices AND
//      bead position, on random triangles whose centroid sits at polar angles below, on and above the ramp (and in the
//      shell at x = 0, +-0.5, +-0.9)
//   3. the same for a whole sphere of faces around the bead (Tot_Energy / Gradient of the interaction), with the
//      bead and the membrane shifted at random
//   4. the energy of a sphere of radius sigma around the bead, normals pointing to the bead, is -W 4 pi sigma^2 f with f
//      the area fraction of the patch (a quadrature of the ramp), within the mesh error
//   5. an Adhesion without the patch constants is not changed (same energy as the explicit formula)
//   6. the free patch axis ("patch_free"): started 12 degrees off, the Newton iteration finds the axis of the covered cap, the
//      energy is the one at that axis (and below the one at the start), and the force with the axis at its maximum matches
//      finite differences of the profiled energy (vertices and bead), because dE/dx = dE/dx at fixed axis there
//
// Usage: Test_patch_adhesion [mesh.obj]   (default ../../../input/Simplesphere.obj)

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "Beads.h"
#include "BeadGeometry.h"
#include "Interaction.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace
{
    bool check(bool ok, const char *what)
    {
        std::printf("%s: %s\n", ok ? "ok  " : "FAIL", what);
        return ok;
    }

    const double PI_ = bead_geometry::PI_VALUE;

    std::vector<double> patchConstants(double W, double a, double s, double thetaDeg, double widthDeg, Vector3 axis)
    {
        const double th = thetaDeg * PI_ / 180.0, w = widthDeg * PI_ / 180.0;
        axis = axis.normalize();
        return {W, a, 0.0, s, 1.0, std::cos(th), std::cos(th - w) - std::cos(th), axis.x, axis.y, axis.z};
    }

    std::vector<double> patchConstantsFree(double W, double a, double s, double thetaDeg, double widthDeg, Vector3 axis)
    {
        std::vector<double> c = patchConstants(W, a, s, thetaDeg, widthDeg, axis);
        c.push_back(1.0);
        return c;
    }

    struct Setup
    {
        std::unique_ptr<ManifoldSurfaceMesh> mesh;
        std::unique_ptr<VertexPositionGeometry> geometry;
        std::unique_ptr<Adhesion> adhesion;
        Bead bead;

        void attach(const std::vector<double> &constants)
        {
            adhesion.reset(new Adhesion(mesh.get(), geometry.get(), constants));
            bead.mesh = mesh.get();
            bead.geometry = geometry.get();
            bead.Pos = Vector3{0.0, 0.0, 0.0};
            bead.sigma = constants[1];
            bead.strength = constants[0];
            bead.Bead_I = adhesion.get();
            adhesion->Bead_1 = &bead;
            bead.CoverageForce = Vector3{0.0, 0.0, 0.0};
        }
    };
} // namespace

int main(int argc, char **argv)
{
    std::string path = argc > 1 ? argv[1] : "../../../input/Simplesphere.obj";
    bool ok = true;

    // ---- 1. the weight ---------------------------------------------------------------------------
    {
        const double thetaDeg = 143.0, widthDeg = 5.0;
        const double th = thetaDeg * PI_ / 180.0, w = widthDeg * PI_ / 180.0;
        const double cosEdge = std::cos(th), dCos = std::cos(th - w) - std::cos(th);
        double d;
        ok &= check(dCos > 0.0, "dCos > 0");
        ok &= check(bead_geometry::patchWeight(cosEdge, cosEdge, dCos, d) == 0.0 && d == 0.0, "weight 0 at the edge, slope 0");
        ok &= check(bead_geometry::patchWeight(cosEdge - 0.3, cosEdge, dCos, d) == 0.0, "weight 0 beyond the edge (inactive cap)");
        ok &= check(bead_geometry::patchWeight(cosEdge + dCos, cosEdge, dCos, d) == 1.0 && d == 0.0, "weight 1 at the end of the ramp, slope 0");
        ok &= check(bead_geometry::patchWeight(1.0, cosEdge, dCos, d) == 1.0, "weight 1 at the pole");
        double dw, wprev = 0.0, worst = 0.0;
        bool monotone = true;
        for (int i = 1; i < 1000; i++)
        {
            const double t = cosEdge + dCos * i / 1000.0;
            const double wv = bead_geometry::patchWeight(t, cosEdge, dCos, dw);
            monotone &= wv >= wprev && dw >= 0.0;
            wprev = wv;
            const double h = 1e-7 * dCos;
            double d1, d2;
            const double wp = bead_geometry::patchWeight(t + h, cosEdge, dCos, d1);
            const double wm = bead_geometry::patchWeight(t - h, cosEdge, dCos, d2);
            worst = std::max(worst, std::abs((wp - wm) / (2 * h) - dw) * dCos);
        }
        ok &= check(monotone, "weight is monotone in t on the ramp");
        char msg[100];
        std::snprintf(msg, sizeof msg, "dw/dt against central differences, worst |error| * dCos = %.2e", worst);
        ok &= check(worst < 1e-6, msg);
    }

    // ---- 2. single triangles ---------------------------------------------------------------------
    std::mt19937 gen(4321);
    std::normal_distribution<double> gauss(0.0, 1.0);
    {
        const double a = 1.0, s = 0.25, W = 1.7;
        const double thetaDeg = 130.0, widthDeg = 8.0;
        const Vector3 axis = unit(Vector3{0.3, -0.5, 0.8});
        std::printf("\nsingle triangles around the patch edge (theta_Wr = %g deg, ramp %g deg), polar angle theta of the centroid\n", thetaDeg, widthDeg);
        std::printf("  theta[deg]   x       E            max|F|      rel.err(membrane)  rel.err(bead)\n");
        double worstAll = 0.0;
        // polar angles: pole, inside the patch, in the ramp (122..130), at the edge, beyond it
        // (a centroid exactly on the edge is left out: the energy is C1 there, a central difference of the quadratic onset is O(h) off)
        for (double thetaPos : {10.0, 60.0, 110.0, 121.0, 123.0, 125.0, 127.0, 129.0, 129.9, 130.1, 131.0, 150.0})
        {
            for (double x : {0.0, -0.5, 0.5, -0.9, 0.9})
            {
                for (int trial = 0; trial < 3; trial++)
                {
                    // a direction at polar angle thetaPos from the axis, random azimuth
                    Vector3 e1 = unit(cross(axis, Vector3{gauss(gen), gauss(gen), gauss(gen)}));
                    Vector3 e2 = cross(axis, e1);
                    const double phi = 2.0 * PI_ * (0.5 + 0.1 * gauss(gen));
                    const double thr = thetaPos * PI_ / 180.0 + 0.002 * trial;
                    Vector3 u = std::cos(thr) * axis + std::sin(thr) * (std::cos(phi) * e1 + std::sin(phi) * e2);
                    u = unit(u);
                    Vector3 t1 = unit(cross(u, Vector3{gauss(gen), gauss(gen), gauss(gen)}));
                    Vector3 t2 = cross(u, t1);
                    const double r = a * (1.0 + s * x);
                    std::vector<Vector3> pos(3);
                    Vector3 mean{0, 0, 0};
                    for (int i = 0; i < 3; i++)
                    {
                        pos[i] = 0.1 * (gauss(gen) * t1 + gauss(gen) * t2) + 0.02 * gauss(gen) * u;
                        mean += pos[i] / 3.0;
                    }
                    Vector3 centre = r * u;
                    for (int i = 0; i < 3; i++)
                        pos[i] = pos[i] - mean + centre;
                    Vector3 n = cross(pos[1] - pos[0], pos[2] - pos[0]);
                    std::vector<std::vector<size_t>> faces = {{0, 1, 2}};
                    if (dot(n, u) > 0)
                        faces = {{0, 2, 1}};
                    const double h = 1e-7;
                    Setup S;
                    S.mesh.reset(new ManifoldSurfaceMesh(faces));
                    S.geometry.reset(new VertexPositionGeometry(*S.mesh));
                    for (Vertex v : S.mesh->vertices())
                        S.geometry->inputVertexPositions[v] = pos[v.getIndex()];
                    S.geometry->refreshQuantities();
                    S.attach(patchConstants(W, a, s, thetaDeg, widthDeg, axis));
                    Face f = S.mesh->face(0);

                    const double E0 = S.adhesion->Face_Energy(f);
                    VertexData<Vector3> F(*S.mesh, Vector3{0, 0, 0});
                    Vector3 Fb{0, 0, 0};
                    S.adhesion->Add_Face_Force(f, F, Fb);

                    double worstM = 0.0, fmax = 0.0;
                    for (Vertex v : S.mesh->vertices())
                        for (int k = 0; k < 3; k++)
                        {
                            Vector3 &q = S.geometry->inputVertexPositions[v];
                            const double old = q[k];
                            q[k] = old + h;
                            double Ep = S.adhesion->Face_Energy(f);
                            q[k] = old - h;
                            double Em = S.adhesion->Face_Energy(f);
                            q[k] = old;
                            double fd = -(Ep - Em) / (2 * h);
                            worstM = std::max(worstM, std::abs(fd - F[v][k]));
                            fmax = std::max(fmax, std::abs(F[v][k]));
                        }
                    double worstB = 0.0;
                    for (int k = 0; k < 3; k++)
                    {
                        const double old = S.bead.Pos[k];
                        S.bead.Pos[k] = old + h;
                        double Ep = S.adhesion->Face_Energy(f);
                        S.bead.Pos[k] = old - h;
                        double Em = S.adhesion->Face_Energy(f);
                        S.bead.Pos[k] = old;
                        double fd = -(Ep - Em) / (2 * h);
                        worstB = std::max(worstB, std::abs(fd - Fb[k]));
                        fmax = std::max(fmax, std::abs(Fb[k]));
                    }
                    const double scale = std::max(fmax, 1e-8);
                    if (trial == 0)
                        std::printf("  %-10g %-7g %-12.5e %-11.4e %-18.3e %.3e\n", thetaPos, x, E0, fmax, worstM / scale, worstB / scale);
                    worstAll = std::max(worstAll, std::max(worstM, worstB) / scale);
                    if (thetaPos >= 130.1)
                        ok &= check(E0 == 0.0 && fmax == 0.0, "beyond the patch edge: no energy and no force");
                }
            }
        }
        std::printf("worst relative error over all single-triangle checks: %.3e\n", worstAll);
        ok &= check(worstAll < 1e-6, "single triangles: force matches finite differences (rel. < 1e-6)");
    }

    // ---- 3. a sphere of faces around the bead, shifted at random ----------------------------------
    {
        const double a = 1.0, s = 0.25, W = 2.3;
        std::printf("\nwhole mesh (radius-1 sphere around the bead, normals towards the bead), patch 120 deg / ramp 10 deg\n");
        Setup S;
        std::tie(S.mesh, S.geometry) = readManifoldSurfaceMesh(path);
        // reverse the orientation: faces looking at the bead
        {
            std::vector<std::vector<size_t>> faces;
            for (Face f : S.mesh->faces())
            {
                Halfedge he = f.halfedge();
                faces.push_back({he.vertex().getIndex(), he.next().next().vertex().getIndex(), he.next().vertex().getIndex()});
            }
            std::vector<Vector3> pos;
            for (Vertex v : S.mesh->vertices())
                pos.push_back(S.geometry->inputVertexPositions[v]);
            S.mesh.reset(new ManifoldSurfaceMesh(faces));
            S.geometry.reset(new VertexPositionGeometry(*S.mesh));
            for (Vertex v : S.mesh->vertices())
                S.geometry->inputVertexPositions[v] = pos[v.getIndex()];
        }
        for (Vertex v : S.mesh->vertices())
            S.geometry->inputVertexPositions[v] = S.geometry->inputVertexPositions[v] * (1.0 + 0.03 * gauss(gen));
        S.geometry->refreshQuantities();
        const Vector3 axis = unit(Vector3{0.2, 0.1, -1.0});
        S.attach(patchConstants(W, a, s, 120.0, 10.0, axis));
        S.bead.Pos = Vector3{0.03, -0.02, 0.01};
        double E0 = S.adhesion->Tot_Energy();
        VertexData<Vector3> F = S.adhesion->Gradient();
        Vector3 Fb = S.bead.Total_force;
        const double h = 1e-6;
        double worst = 0.0, fmax = 0.0;
        for (Vertex v : S.mesh->vertices())
            for (int k = 0; k < 3; k++)
            {
                Vector3 &q = S.geometry->inputVertexPositions[v];
                const double old = q[k];
                q[k] = old + h;
                double Ep = S.adhesion->Tot_Energy();
                q[k] = old - h;
                double Em = S.adhesion->Tot_Energy();
                q[k] = old;
                worst = std::max(worst, std::abs(-(Ep - Em) / (2 * h) - F[v][k]));
                fmax = std::max(fmax, std::abs(F[v][k]));
            }
        double wb = 0.0;
        for (int k = 0; k < 3; k++)
        {
            const double old = S.bead.Pos[k];
            S.bead.Pos[k] = old + h;
            double Ep = S.adhesion->Tot_Energy();
            S.bead.Pos[k] = old - h;
            double Em = S.adhesion->Tot_Energy();
            S.bead.Pos[k] = old;
            wb = std::max(wb, std::abs(-(Ep - Em) / (2 * h) - Fb[k]));
            fmax = std::max(fmax, std::abs(Fb[k]));
        }
        std::printf("  E = %.8g, max|F| %.3e, rel.err membrane %.2e, bead %.2e\n", E0, fmax, worst / fmax, wb / fmax);
        ok &= check(E0 < 0.0 && worst < 1e-6 * fmax && wb < 1e-6 * fmax, "whole mesh: forces match finite differences");
    }

    // ---- 4. the area of the patch and 5. no patch: unchanged --------------------------------------
    {
        const double a = 1.0, s = 0.25, W = 1.0;
        std::printf("\nsphere of radius sigma around the bead (the shell weight is 1 on the sphere, the faces' centroids sit slightly inside)\n");
        Setup S;
        std::tie(S.mesh, S.geometry) = readManifoldSurfaceMesh(path);
        {
            std::vector<std::vector<size_t>> faces;
            for (Face f : S.mesh->faces())
            {
                Halfedge he = f.halfedge();
                faces.push_back({he.vertex().getIndex(), he.next().next().vertex().getIndex(), he.next().vertex().getIndex()});
            }
            std::vector<Vector3> pos;
            for (Vertex v : S.mesh->vertices())
                pos.push_back(unit(S.geometry->inputVertexPositions[v]));
            S.mesh.reset(new ManifoldSurfaceMesh(faces));
            S.geometry.reset(new VertexPositionGeometry(*S.mesh));
            for (Vertex v : S.mesh->vertices())
                S.geometry->inputVertexPositions[v] = pos[v.getIndex()];
        }
        S.geometry->refreshQuantities();
        const Vector3 axis{0.0, 0.0, 1.0};
        // the reference without patch: the faces' share of the sphere
        S.attach({W, a, 0.0, s, 1.0});
        const double Efull = S.adhesion->Tot_Energy();
        for (double thetaDeg : {60.0, 90.0, 143.13, 165.0, 180.0})
        {
            const double widthDeg = 6.0;
            S.attach(patchConstants(W, a, s, thetaDeg, widthDeg, axis));
            const double E = S.adhesion->Tot_Energy();
            // quadrature of the area fraction of the weight on the sphere
            const double th = thetaDeg * PI_ / 180.0, wd = widthDeg * PI_ / 180.0;
            double frac = 0.0;
            const int N = 20000;
            for (int i = 0; i < N; i++)
            {
                const double theta = PI_ * (i + 0.5) / N; // polar angle from the axis
                const double t = std::cos(theta);
                double d;
                frac += bead_geometry::patchWeight(t, std::cos(th), std::cos(th - wd) - std::cos(th), d) * std::sin(theta) * PI_ / N / 2.0;
            }
            const double ratio = E / Efull, expected = frac;
            std::printf("  theta_Wr = %-7g E/E(no patch) = %.5f, quadrature %.5f (area fraction of the cap %.5f)\n", thetaDeg, ratio, expected,
                        (1.0 - std::cos(th)) / 2.0);
            char msg[120];
            std::snprintf(msg, sizeof msg, "theta_Wr = %g: energy ratio matches the quadrature of the weight within 2e-2", thetaDeg);
            ok &= check(std::abs(ratio - expected) < 2e-2, msg);
        }
        // 5. the old formula, written out: E = -W sigma^2 sum w |omega|
        S.attach({W, a, 0.0, s, 1.0});
        double Eold = 0.0;
        for (Face f : S.mesh->faces())
        {
            double dwdr;
            const bead_geometry::FaceCoverageData face = bead_geometry::evaluateFaceCoverage(*S.geometry, f, S.bead.Pos, a, s, 1.0);
            if (face.selected)
                Eold += -W * a * a * face.weight * face.omegaUnsigned;
        }
        ok &= check(std::abs(S.adhesion->Tot_Energy() - Eold) <= 1e-12 * std::abs(Eold), "no patch constants: the energy is the old formula");
    }

    // ---- 6. the free patch axis -------------------------------------------------------------------
    {
        const double a = 1.0, s = 0.25, W = 2.3;
        std::printf("\nfree patch axis: sphere around the bead that lifts off beyond 100 deg of the true axis (unique maximum)\n");
        Setup S;
        std::tie(S.mesh, S.geometry) = readManifoldSurfaceMesh(path);
        {
            std::vector<std::vector<size_t>> faces;
            for (Face f : S.mesh->faces())
            {
                Halfedge he = f.halfedge();
                faces.push_back({he.vertex().getIndex(), he.next().next().vertex().getIndex(), he.next().vertex().getIndex()});
            }
            std::vector<Vector3> pos;
            for (Vertex v : S.mesh->vertices())
                pos.push_back(unit(S.geometry->inputVertexPositions[v]));
            S.mesh.reset(new ManifoldSurfaceMesh(faces));
            S.geometry.reset(new VertexPositionGeometry(*S.mesh));
            const Vector3 trueAxis = unit(Vector3{0.3, -0.2, 0.9});
            for (Vertex v : S.mesh->vertices())
            {
                Vector3 p = pos[v.getIndex()];
                // on the sphere up to 100 deg of the axis, then the membrane lifts off smoothly (the weight falls with the radius),
                // so the maximum of the covered weight over the axis is unique (a cap that fits inside the patch would be a plateau)
                const double polar = std::acos(std::max(-1.0, std::min(1.0, dot(p, trueAxis))));
                const double y = std::min(1.0, std::max(0.0, (polar - 100.0 * PI_ / 180.0) / (60.0 * PI_ / 180.0)));
                S.geometry->inputVertexPositions[v] = (1.0 + 0.35 * y * y * (3.0 - 2.0 * y)) * p;
            }
            S.geometry->refreshQuantities();
            // start the axis 12 degrees away from the true one
            Vector3 off = unit(cross(trueAxis, Vector3{0.2, 1.0, -0.3}));
            const double tilt = 12.0 * PI_ / 180.0;
            const Vector3 start = std::cos(tilt) * trueAxis + std::sin(tilt) * off;
            S.attach(patchConstantsFree(W, a, s, 120.0, 10.0, start));
            const double Efixed = [&]
            {
                Setup T;
                T.mesh.reset(new ManifoldSurfaceMesh(S.mesh->getFaceVertexList()));
                T.geometry.reset(new VertexPositionGeometry(*T.mesh));
                for (Vertex v : T.mesh->vertices())
                    T.geometry->inputVertexPositions[v] = S.geometry->inputVertexPositions[v];
                T.geometry->refreshQuantities();
                T.attach(patchConstants(W, a, s, 120.0, 10.0, start));
                return T.adhesion->Tot_Energy();
            }();
            const double Efree = S.adhesion->Tot_Energy();
            const Vector3 found = S.adhesion->PatchAxis();
            const double err = std::acos(std::max(-1.0, std::min(1.0, dot(found, trueAxis))));
            std::printf("  start tilt 12 deg, found axis %.3e rad from the true one, E(start axis) = %.8f, E(free) = %.8f\n", err, Efixed, Efree);
            // the mesh is not symmetric about the axis, so the maximum sits a little off the true axis: check it is near, and that
            // it is a stationary point (a first difference of the energy under tilts of the axis vanishes, the second is positive)
            ok &= check(err < 0.05, "the free axis is near the axis of the covered cap");
            ok &= check(Efree <= Efixed, "the energy at the free axis is below the one at the start axis");
            {
                const Vector3 e1 = unit(cross(found, Vector3{0.0, 1.0, 0.3}));
                const Vector3 e2 = cross(found, e1);
                auto energyAt = [&](const Vector3 &axisHere)
                {
                    Setup T;
                    T.mesh.reset(new ManifoldSurfaceMesh(S.mesh->getFaceVertexList()));
                    T.geometry.reset(new VertexPositionGeometry(*T.mesh));
                    for (Vertex v : T.mesh->vertices())
                        T.geometry->inputVertexPositions[v] = S.geometry->inputVertexPositions[v];
                    T.geometry->refreshQuantities();
                    T.attach(patchConstants(W, a, s, 120.0, 10.0, axisHere));
                    return T.adhesion->Tot_Energy();
                };
                const double d = 1e-4, E0 = energyAt(found);
                double firstWorst = 0.0, secondMin = 1e30;
                for (const Vector3 &e : {e1, e2})
                {
                    const double Ep = energyAt(unit(found + d * e)), Em = energyAt(unit(found - d * e));
                    firstWorst = std::max(firstWorst, std::abs(Ep - Em) / (2 * d));
                    secondMin = std::min(secondMin, (Ep - 2 * E0 + Em) / (d * d));
                }
                std::printf("  at the found axis: |dE/dtilt| = %.2e, d2E/dtilt2 = %.3f (E0 = %.8f)\n", firstWorst, secondMin, E0);
                ok &= check(firstWorst < 1e-7 && secondMin > 0.0, "the free axis is a stationary minimum of the energy over the axis");
            }

            // forces of the profiled energy against finite differences
            S.bead.Pos = Vector3{0.03, -0.02, 0.01};
            for (Vertex v : S.mesh->vertices())
                S.geometry->inputVertexPositions[v] = S.geometry->inputVertexPositions[v] +
                                                      0.02 * Vector3{gauss(gen), gauss(gen), gauss(gen)};
            S.geometry->refreshQuantities();
            VertexData<Vector3> F = S.adhesion->Gradient();
            Vector3 Fb = S.bead.Total_force;
            const double h = 1e-6;
            double worst = 0.0, fmax = 0.0;
            for (Vertex v : S.mesh->vertices())
                for (int k = 0; k < 3; k++)
                {
                    Vector3 &q = S.geometry->inputVertexPositions[v];
                    const double old = q[k];
                    q[k] = old + h;
                    double Ep = S.adhesion->Tot_Energy();
                    q[k] = old - h;
                    double Em = S.adhesion->Tot_Energy();
                    q[k] = old;
                    worst = std::max(worst, std::abs(-(Ep - Em) / (2 * h) - F[v][k]));
                    fmax = std::max(fmax, std::abs(F[v][k]));
                }
            double wb = 0.0;
            for (int k = 0; k < 3; k++)
            {
                const double old = S.bead.Pos[k];
                S.bead.Pos[k] = old + h;
                double Ep = S.adhesion->Tot_Energy();
                S.bead.Pos[k] = old - h;
                double Em = S.adhesion->Tot_Energy();
                S.bead.Pos[k] = old;
                wb = std::max(wb, std::abs(-(Ep - Em) / (2 * h) - Fb[k]));
                fmax = std::max(fmax, std::abs(Fb[k]));
            }
            std::printf("  profiled energy: max|F| %.3e, rel.err membrane %.2e, bead %.2e\n", fmax, worst / fmax, wb / fmax);
            ok &= check(worst < 1e-6 * fmax && wb < 1e-6 * fmax, "free axis: forces match finite differences of the profiled energy");
        }
    }

    std::printf(ok ? "PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
