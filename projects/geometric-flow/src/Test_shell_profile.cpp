// Checks for the radial profiles of the Adhesion weight: Frenkel, Lennard-Jones and Morse
// (bead_geometry::ShellProfile, "shell_profile" in the input file).
//   1. the weight: w(a) = 1 and w'(a) = 0 (the minimum of the well, a = bead radius), w = 0 from rc = a (1 + s) on,
//      w < 1 elsewhere, the repulsive core is negative, the potential formulas against independent evaluations,
//      Frenkel has zero slope at rc, dw/dr against central differences, the Cosine kind is bit-identical to
//      coverageShellWeight, invalid widths / alpha are rejected
//   2. Adhesion::Add_Face_Force against central finite differences of Adhesion::Face_Energy (vertices and bead),
//      random triangles at several r / a from deep in the core to the cutoff
//   3. the same for a whole patch of a sphere (Tot_Energy / Gradient)
//
// Usage: Test_shell_profile [mesh.obj]   (default ../../../input/Simplesphere.obj)

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
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
using bead_geometry::ShellKind;
using bead_geometry::ShellProfile;

namespace
{
    bool check(bool ok, const char *what)
    {
        std::printf("%s: %s\n", ok ? "ok  " : "FAIL", what);
        return ok;
    }

    ShellProfile make(ShellKind kind, double width, double alpha = 6.0)
    {
        ShellProfile p;
        p.kind = kind;
        p.width = width;
        p.alpha = alpha;
        return p;
    }

    // Independent evaluation of the potentials (depth 1 at r = a), without the shift
    double frenkelV(double r, double a, double rc)
    {
        const double s2 = 2.0 * rc * rc * a * a / (3.0 * rc * rc - a * a);
        const double alpha = 2.0 * (rc * rc / s2) * std::pow(3.0 / (2.0 * (rc * rc / s2 - 1.0)), 3.0);
        return alpha * (s2 / (r * r) - 1.0) * std::pow(rc * rc / (r * r) - 1.0, 2.0);
    }
    double ljV(double r, double a) { return std::pow(a / r, 12) - 2.0 * std::pow(a / r, 6); }
    double morseV(double r, double a, double alpha) { return std::pow(1.0 - std::exp(-alpha * (r - a) / a), 2.0) - 1.0; }

    struct Setup
    {
        std::unique_ptr<ManifoldSurfaceMesh> mesh;
        std::unique_ptr<VertexPositionGeometry> geometry;
        std::unique_ptr<Adhesion> adhesion;
        Bead bead;

        void attach(double W, double a, const ShellProfile &profile)
        {
            adhesion.reset(new Adhesion(mesh.get(), geometry.get(), {W, a, 0.0, profile.width}));
            adhesion->shellKind = profile.kind;
            adhesion->shellAlpha = profile.alpha;
            bead.mesh = mesh.get();
            bead.geometry = geometry.get();
            bead.Pos = Vector3{0.0, 0.0, 0.0};
            bead.sigma = a;
            bead.strength = W;
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

    struct Case
    {
        ShellKind kind;
        double width;
        double alpha;
        const char *name;
    };
    const std::vector<Case> cases = {
        {ShellKind::Frenkel, 0.25, 6.0, "frenkel s=0.25"}, {ShellKind::Frenkel, 1.0, 6.0, "frenkel s=1"},
        {ShellKind::LJ, 0.25, 6.0, "lj s=0.25"},           {ShellKind::LJ, 1.0, 6.0, "lj s=1"},
        {ShellKind::Morse, 0.25, 6.0, "morse s=0.25 alpha=6"}, {ShellKind::Morse, 1.0, 3.0, "morse s=1 alpha=3"},
        // rc = 2.5 and 2 sigma_LJ (sigma_LJ = a / 2^(1/6)) and the Morse wells of the planar potential scan (a = 2, 10 are the ones used)
        {ShellKind::LJ, 1.2272467953508484, 6.0, "lj rc=2.5 sigma"}, {ShellKind::LJ, 0.7817974362806785, 6.0, "lj rc=2 sigma"},
        {ShellKind::Morse, 1.2272467953508484, 0.25, "morse a=0.25 rc=2.5 sigma"}, {ShellKind::Morse, 1.2272467953508484, 0.1, "morse a=0.1 rc=2.5 sigma"},
        {ShellKind::Morse, 1.2272467953508484, 2.0, "morse a=2 rc=2.5 sigma"}, {ShellKind::Morse, 1.2272467953508484, 10.0, "morse a=10 rc=2.5 sigma"}};

    // ---- 1. the weight ---------------------------------------------------------------------------
    {
        const double a = 1.3;
        for (const Case &c : cases)
        {
            char msg[200];
            const ShellProfile p = make(c.kind, c.width, c.alpha);
            const double rc = a * (1.0 + c.width);
            double dw;
            const double w0 = bead_geometry::shellProfileWeight(a, a, dw, p);
            std::snprintf(msg, sizeof msg, "%s: w(a) = %.15g, w'(a) = %.2e", c.name, w0, dw);
            ok &= check(std::abs(w0 - 1.0) < 1e-12 && std::abs(dw) < 1e-9, msg);

            double wc = bead_geometry::shellProfileWeight(rc, a, dw, p);
            std::snprintf(msg, sizeof msg, "%s: w(rc) = %g exactly zero, w = 0 beyond", c.name, wc);
            ok &= check(wc == 0.0 && bead_geometry::shellProfileWeight(rc * 1.01, a, dw, p) == 0.0 && dw == 0.0, msg);

            double dIn;
            const double wIn = bead_geometry::shellProfileWeight(rc * (1.0 - 1e-9), a, dIn, p);
            std::snprintf(msg, sizeof msg, "%s: just inside rc w = %.2e, w' = %.3e%s", c.name, wIn, dIn,
                          c.kind == ShellKind::Frenkel ? " (zero slope)" : " (shifted: the slope is finite)");
            ok &= check(std::abs(wIn) < 1e-6 && (c.kind != ShellKind::Frenkel || std::abs(dIn) < 1e-6), msg);

            // w <= 1 everywhere (the minimum of the well), attractive tail, repulsive core
            double wmax = -1e300;
            bool tailPositive = true, corePositive = false;
            for (int i = 0; i < 4000; i++)
            {
                const double r = a * (0.3 + i * (1.0 + c.width - 0.3) / 4000.0);
                if (r >= rc)
                    break;
                const double w = bead_geometry::shellProfileWeight(r, a, dw, p);
                wmax = std::max(wmax, w);
                if (r > a && w < 0.0)
                    tailPositive = false;
                if (r < 0.8 * a && w < 0.0)
                    corePositive = true;
            }
            // a Morse well that wide (alpha < 1) has no repulsive core at r >= 0.3 a, only a weaker attraction
            const bool softMorse = c.kind == ShellKind::Morse && c.alpha < 1.0;
            std::snprintf(msg, sizeof msg, "%s: max w = %.12g on a grid, w > 0 for a < r < rc%s", c.name, wmax,
                          softMorse ? " (no core: alpha < 1)" : ", repulsive (w < 0) at r < 0.8 a");
            ok &= check(wmax <= 1.0 + 1e-12 && tailPositive && (corePositive || softMorse), msg);

            // against the independent formulas
            double worstF = 0.0;
            for (double x : {0.7, 0.85, 0.95, 1.0, 1.01, 1.05, 1.1, 1.2})
            {
                const double r = a * std::min(x, 1.0 + 0.9 * c.width);
                double V, Vc, depth;
                switch (c.kind)
                {
                case ShellKind::Frenkel:
                    V = frenkelV(r, a, rc), Vc = 0.0;
                    break;
                case ShellKind::LJ:
                    V = ljV(r, a), Vc = ljV(rc, a);
                    break;
                default:
                    V = morseV(r, a, c.alpha), Vc = morseV(rc, a, c.alpha);
                    break;
                }
                depth = 1.0 + Vc;
                const double ref = (Vc - V) / depth;
                const double w = bead_geometry::shellProfileWeight(r, a, dw, p);
                worstF = std::max(worstF, std::abs(w - ref) / std::max(1.0, std::abs(ref)));
            }
            std::snprintf(msg, sizeof msg, "%s: w against the independent potential formula, worst error %.2e", c.name, worstF);
            ok &= check(worstF < 1e-12, msg);

            // derivative against central differences
            double worst = 0.0;
            for (double x : {0.5, 0.7, 0.85, 0.95, 0.99, 1.0, 1.01, 1.05, 1.5 })
            {
                const double r = a * std::min(x, 1.0 + 0.999 * c.width), h = 1e-6 * a;
                if (r + h >= rc)
                    continue;
                double dd, d1, d2;
                bead_geometry::shellProfileWeight(r, a, dd, p);
                const double wp = bead_geometry::shellProfileWeight(r + h, a, d1, p);
                const double wm = bead_geometry::shellProfileWeight(r - h, a, d2, p);
                worst = std::max(worst, std::abs((wp - wm) / (2 * h) - dd) / std::max(1e-3, std::abs(dd)));
            }
            std::snprintf(msg, sizeof msg, "%s: dw/dr against central differences, worst relative error %.2e", c.name, worst);
            ok &= check(worst < 1e-6, msg);
        }

        // the Cosine kind is the old shell, bit for bit
        {
            ShellProfile p;
            p.width = 0.3;
            p.power = 2.0;
            long diff = 0;
            for (int i = -1500; i <= 1500; i++)
            {
                const double r = a * (1.0 + 0.0004 * i);
                double d1, d2;
                const double w1 = bead_geometry::shellProfileWeight(r, a, d1, p);
                const double w2 = bead_geometry::coverageShellWeight(r, a, d2, 0.3, 2.0);
                diff += !(w1 == w2 && d1 == d2);
            }
            ok &= check(diff == 0, "Cosine kind: bit-identical to coverageShellWeight");
        }

        // names and invalid parameters
        {
            ok &= check(bead_geometry::shellKindFromName("lj") == ShellKind::LJ && bead_geometry::shellKindFromName("cosine") == ShellKind::Cosine &&
                            bead_geometry::shellKindFromName("morse") == ShellKind::Morse && bead_geometry::shellKindFromName("frenkel") == ShellKind::Frenkel,
                        "profile names");
            bool threw = false;
            try
            {
                bead_geometry::shellKindFromName("gauss");
            }
            catch (const std::invalid_argument &)
            {
                threw = true;
            }
            ok &= check(threw, "unknown profile name is rejected");
            for (double badWidth : {0.0, -0.1, 4.5, std::numeric_limits<double>::quiet_NaN()})
            {
                threw = false;
                try
                {
                    double d;
                    bead_geometry::shellProfileWeight(a, a, d, make(ShellKind::LJ, badWidth));
                }
                catch (const std::invalid_argument &)
                {
                    threw = true;
                }
                char msg[80];
                std::snprintf(msg, sizeof msg, "width = %g is rejected", badWidth);
                ok &= check(threw, msg);
            }
            for (double badAlpha : {0.0, -2.0, std::numeric_limits<double>::infinity()})
            {
                threw = false;
                try
                {
                    double d;
                    bead_geometry::shellProfileWeight(a, a, d, make(ShellKind::Morse, 0.25, badAlpha));
                }
                catch (const std::invalid_argument &)
                {
                    threw = true;
                }
                char msg[80];
                std::snprintf(msg, sizeof msg, "Morse alpha = %g is rejected", badAlpha);
                ok &= check(threw, msg);
            }
        }
    }

    // ---- 2. single triangles ---------------------------------------------------------------------
    std::mt19937 gen(4242);
    std::normal_distribution<double> gauss(0.0, 1.0);
    {
        const double a = 1.0, W = 1.7;
        std::printf("\nsingle triangles, relative error = max|FD - analytic| / max|analytic|\n");
        for (const Case &c : cases)
        {
            const ShellProfile p = make(c.kind, c.width, c.alpha);
            double worstAll = 0.0;
            int nsel = 0;
            for (double rr : {0.6, 0.8, 0.95, 1.0, 1.05, 1.1, 1.2, 0.9 + 0.0, 1.0 + 0.5 * c.width, 1.0 + 0.9 * c.width, 1.0 + 0.995 * c.width})
            {
                for (int trial = 0; trial < 3; trial++)
                {
                    Vector3 u = unit(Vector3{gauss(gen), gauss(gen), gauss(gen)});
                    Vector3 t1 = unit(cross(u, Vector3{gauss(gen), gauss(gen), gauss(gen)}));
                    Vector3 t2 = cross(u, t1);
                    const double r = a * rr;
                    if (r >= a * (1.0 + c.width))
                        continue;
                    std::vector<Vector3> pos(3);
                    Vector3 mean{0, 0, 0};
                    for (int i = 0; i < 3; i++)
                    {
                        pos[i] = 0.2 * (gauss(gen) * t1 + gauss(gen) * t2) + 0.02 * gauss(gen) * u;
                        mean += pos[i] / 3.0;
                    }
                    Vector3 centre = r * u;
                    for (int i = 0; i < 3; i++)
                        pos[i] = pos[i] - mean + centre;
                    Vector3 n = cross(pos[1] - pos[0], pos[2] - pos[0]);
                    std::vector<std::vector<size_t>> faces = {{0, 1, 2}};
                    if (dot(n, u) > 0)
                        faces = {{0, 2, 1}};
                    // FD step: far from the cutoff, where the shifted potentials have a kink in the slope
                    const double h = std::min(1e-6, 1e-3 * (a * (1.0 + c.width) - r));
                    Setup S;
                    S.mesh.reset(new ManifoldSurfaceMesh(faces));
                    S.geometry.reset(new VertexPositionGeometry(*S.mesh));
                    for (Vertex v : S.mesh->vertices())
                        S.geometry->inputVertexPositions[v] = pos[v.getIndex()];
                    S.geometry->refreshQuantities();
                    S.attach(W, a, p);
                    Face f = S.mesh->face(0);

                    const double E0 = S.adhesion->Face_Energy(f);
                    if (E0 != 0.0)
                        nsel++;
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
                            worstM = std::max(worstM, std::abs(-(Ep - Em) / (2 * h) - F[v][k]));
                            fmax = std::max(fmax, std::abs(F[v][k]));
                        }
                    for (int k = 0; k < 3; k++)
                    {
                        const double old = S.bead.Pos[k];
                        S.bead.Pos[k] = old + h;
                        double Ep = S.adhesion->Face_Energy(f);
                        S.bead.Pos[k] = old - h;
                        double Em = S.adhesion->Face_Energy(f);
                        S.bead.Pos[k] = old;
                        worstM = std::max(worstM, std::abs(-(Ep - Em) / (2 * h) - Fb[k]));
                        fmax = std::max(fmax, std::abs(Fb[k]));
                    }
                    // near rc the force of a stiff Morse is ~1e-6 and the finite differences carry ~1e-12 of roundoff: the floor
                    // keeps that from being read as a relative error (the triangles at smaller r check the formulas at full scale)
                    worstAll = std::max(worstAll, worstM / std::max(fmax, 1e-4));
                }
            }
            char msg[160];
            std::snprintf(msg, sizeof msg, "%s: %d faces selected, force matches finite differences (worst rel. %.2e)", c.name, nsel, worstAll);
            ok &= check(nsel > 20 && worstAll < 1e-6, msg);
        }

        // the face energy has the sign of the weight: negative weight (core) -> positive energy
        {
            const ShellProfile p = make(ShellKind::LJ, 0.5);
            Vector3 u = unit(Vector3{0.2, 0.3, 0.9});
            Vector3 t1 = unit(cross(u, Vector3{1, 0, 0})), t2 = cross(u, t1);
            double Ein = 0.0, Eout = 0.0;
            for (double r : {0.8, 1.0})
            {
                std::vector<std::vector<size_t>> faces = {{0, 2, 1}};
                Setup S;
                S.mesh.reset(new ManifoldSurfaceMesh(faces));
                S.geometry.reset(new VertexPositionGeometry(*S.mesh));
                Vector3 c = r * u;
                std::vector<Vector3> pos = {c + 0.2 * t1, c - 0.2 * t1 + 0.2 * t2, c - 0.2 * t1 - 0.2 * t2};
                for (Vertex v : S.mesh->vertices())
                    S.geometry->inputVertexPositions[v] = pos[v.getIndex()];
                S.geometry->refreshQuantities();
                S.attach(1.0, 1.0, p);
                (r < 0.9 ? Ein : Eout) = S.adhesion->Face_Energy(S.mesh->face(0));
            }
            char msg[120];
            std::snprintf(msg, sizeof msg, "LJ core is repulsive (E(0.8 a) = %.3g > 0) and the well attractive (E(a) = %.3g < 0)", Ein, Eout);
            ok &= check(Ein > 0.0 && Eout < 0.0, msg);
        }
    }

    // ---- 3. a sphere patch -----------------------------------------------------------------------
    {
        const double a = 1.0, W = 1.0;
        std::printf("\nsphere patch (radius 3 sphere, bead outside at r = a from the pole)\n");
        for (const Case &c : cases)
        {
            const ShellProfile p = make(c.kind, c.width, c.alpha);
            Setup S;
            std::tie(S.mesh, S.geometry) = readManifoldSurfaceMesh(path);
            for (Vertex v : S.mesh->vertices())
                S.geometry->inputVertexPositions[v] = 3.0 * S.geometry->inputVertexPositions[v] +
                                                      0.04 * Vector3{gauss(gen), gauss(gen), gauss(gen)};
            S.geometry->refreshQuantities();
            S.attach(W, a, p);
            S.bead.Pos = Vector3{0.15, -0.1, 3.0 + 1.0 * a};
            const double E0 = S.adhesion->Tot_Energy();
            int nsel = 0;
            for (Face f : S.mesh->faces())
                if (S.adhesion->Face_Energy(f) != 0.0)
                    nsel++;
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
            std::printf("  %-22s E = %.8g, %d faces selected, max|F| %.3e, rel.err membrane %.2e, bead %.2e\n", c.name, E0, nsel, fmax,
                        worst / fmax, wb / fmax);
            char msg[100];
            std::snprintf(msg, sizeof msg, "%s patch: forces match finite differences", c.name);
            ok &= check(nsel > 5 && worst < 1e-6 * fmax && wb < 1e-6 * fmax, msg);
        }
    }

    std::printf(ok ? "PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
