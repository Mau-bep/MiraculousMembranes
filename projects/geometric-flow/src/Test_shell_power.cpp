// Checks for the shell exponent p of the Adhesion shell,
//   w(r) = [ 1/2 (1 + cos(pi x)) ]^p,   x = (r - a) / (s a)   (p = 1: the plain shell)
//   1. the weight: w(0) = 1, w(+-1) = 0, w'(+-1) = 0, p = 1 is bit-identical to the old formula,
//      dw/dr against central differences, p < 1 and non-finite p are rejected
//   2. Adhesion::Add_Face_Force against central finite differences of Adhesion::Face_Energy,
//      vertices AND bead position, on random triangles whose centroid sits at x = 0, +-0.3, +-0.6, +-0.9,
//      +-0.99, +-0.999 inside the shell and also on a triangle almost edge-on to the bead
//      (face-selection boundary); p = 1, 1.5, 2, 4
//   3. the same for a whole patch of a sphere (Tot_Energy / Gradient of the interaction)
//
// Usage: Test_shell_power [mesh.obj]   (default ../../../input/Simplesphere.obj)

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

namespace
{
    bool check(bool ok, const char *what)
    {
        std::printf("%s: %s\n", ok ? "ok  " : "FAIL", what);
        return ok;
    }

    // The shell as it was before the exponent existed
    double oldWeight(double r, double sigma, double &dw, double s)
    {
        const double delta = s * sigma;
        const double x = (r - sigma) / delta;
        dw = 0.0;
        if (std::abs(x) >= 1.0)
            return 0.0;
        const double w = 0.5 * (1.0 + std::cos(bead_geometry::PI_VALUE * x));
        dw = -0.5 * bead_geometry::PI_VALUE * std::sin(bead_geometry::PI_VALUE * x) / delta;
        return w;
    }

    struct Setup
    {
        std::unique_ptr<ManifoldSurfaceMesh> mesh;
        std::unique_ptr<VertexPositionGeometry> geometry;
        std::unique_ptr<Adhesion> adhesion;
        Bead bead;

        void attach(double W, double a, double s, double p)
        {
            adhesion.reset(new Adhesion(mesh.get(), geometry.get(), {W, a, 0.0, s, p}));
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
    const double PI = bead_geometry::PI_VALUE;

    // ---- 1. the weight ---------------------------------------------------------------------------
    {
        const double a = 1.3, s = 0.25, delta = s * a;
        double dw;
        for (double p : {1.0, 1.5, 2.0, 4.0})
        {
            char msg[160];
            double w0 = bead_geometry::coverageShellWeight(a, a, dw, s, p);
            double d0 = dw;
            std::snprintf(msg, sizeof msg, "p=%g: w(0) = 1 exactly (got %.17g), w'(0) = %.3g", p, w0, d0);
            ok &= check(w0 == 1.0 && std::abs(d0) < 1e-14, msg);
            for (double sgn : {-1.0, 1.0})
            {
                double wE = bead_geometry::coverageShellWeight(a + sgn * delta, a, dw, s, p);
                std::snprintf(msg, sizeof msg, "p=%g: w(%+d) = %g, w'(%+d) = %g", p, (int)sgn, wE, (int)sgn, dw);
                // r = a +- delta can land one ulp inside the shell (x = 1 - 1e-16): then only roundoff is left
                ok &= check(wE < 1e-30 && std::abs(dw) < 1e-12, msg);
            }
            // slope just inside the edge is small and goes to zero (p >= 1)
            double dwin, win;
            win = bead_geometry::coverageShellWeight(a + delta * (1.0 - 1e-6), a, dwin, s, p);
            std::snprintf(msg, sizeof msg, "p=%g: just inside the edge w = %.3e, w' = %.3e", p, win, dwin);
            ok &= check(win >= 0.0 && std::abs(dwin) < 1e-3, msg);

            // derivative against central differences
            double worst = 0.0;
            for (double x : {-0.999, -0.9, -0.5, -0.2, 0.0, 0.1, 0.4, 0.7, 0.95, 0.999})
            {
                const double r = a + x * delta, h = 1e-7 * delta;
                double dd, d1, d2;
                bead_geometry::coverageShellWeight(r, a, dd, s, p);
                double wp = bead_geometry::coverageShellWeight(r + h, a, d1, s, p);
                double wm = bead_geometry::coverageShellWeight(r - h, a, d2, s, p);
                worst = std::max(worst, std::abs((wp - wm) / (2 * h) - dd) / std::max(1e-3 / delta, std::abs(dd)));
            }
            std::snprintf(msg, sizeof msg, "p=%g: dw/dr vs central differences, worst relative error %.2e", p, worst);
            ok &= check(worst < 1e-6, msg);
        }
        // p = 1 default and explicit are bit-identical to the old formula, also on a fine r grid
        long diff = 0, n = 0;
        for (double ss : {0.1, 0.25, 0.5})
            for (int i = -1500; i <= 1500; i++)
            {
                const double r = a * (1.0 + 0.001 * i * ss * 1.3);
                double dold, dnew, dnew2;
                double wold = oldWeight(r, a, dold, ss);
                double wnew = bead_geometry::coverageShellWeight(r, a, dnew, ss);
                double wnew2 = bead_geometry::coverageShellWeight(r, a, dnew2, ss, 1.0);
                diff += !(wold == wnew && dold == dnew && wold == wnew2 && dold == dnew2);
                n++;
            }
        std::printf("      p=1 against the old formula on %ld points: %ld differences\n", n, diff);
        ok &= check(diff == 0, "p = 1 is bit-identical to the old weight and slope");
        // weight is cos^(2p)(pi x / 2)
        {
            double d;
            double w = bead_geometry::coverageShellWeight(a + 0.37 * delta, a, d, s, 3.0);
            double ref = std::pow(std::cos(PI * 0.37 / 2.0), 6.0);
            ok &= check(std::abs(w - ref) < 1e-14, "w = cos^(2p)(pi x / 2) for p = 3");
        }
        // invalid exponents
        for (double bad : {0.5, 0.0, -1.0, 0.999999, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        {
            bool threw = false;
            try
            {
                double d;
                bead_geometry::coverageShellWeight(a, a, d, s, bad);
            }
            catch (const std::invalid_argument &)
            {
                threw = true;
            }
            char msg[100];
            std::snprintf(msg, sizeof msg, "shellPower = %g is rejected", bad);
            ok &= check(threw, msg);
        }
    }

    // ---- 2. single triangles ---------------------------------------------------------------------
    std::mt19937 gen(12345);
    std::normal_distribution<double> gauss(0.0, 1.0);
    {
        const double a = 1.0, s = 0.25, W = 1.7;
        std::printf("\nsingle triangles, FD step min(1e-6, 1e-3 * distance to the shell edge), relative error = max|FD - analytic| / max|analytic|\n");
        std::printf("  p     x        E            max|F|      rel.err(membrane)  rel.err(bead)\n");
        double worstAll = 0.0;
        for (double p : {1.0, 1.5, 2.0, 4.0})
        {
            for (double x : {0.0, -0.3, 0.3, -0.6, 0.6, -0.9, 0.9, -0.99, 0.99, -0.999, 0.999, -0.9999, 0.9999})
            {
                for (int trial = 0; trial < 3; trial++)
                {
                    Vector3 u = unit(Vector3{gauss(gen), gauss(gen), gauss(gen)});
                    Vector3 t1 = unit(cross(u, Vector3{gauss(gen), gauss(gen), gauss(gen)}));
                    Vector3 t2 = cross(u, t1);
                    const double r = a * (1.0 + s * x);
                    std::vector<Vector3> pos(3);
                    Vector3 mean{0, 0, 0};
                    for (int i = 0; i < 3; i++)
                    {
                        pos[i] = 0.25 * (gauss(gen) * t1 + gauss(gen) * t2) + 0.03 * gauss(gen) * u;
                        mean += pos[i] / 3.0;
                    }
                    Vector3 centre = r * u;
                    for (int i = 0; i < 3; i++)
                        pos[i] = pos[i] - mean + centre;
                    // outward normal must point to the bead (antiparallel to u)
                    Vector3 n = cross(pos[1] - pos[0], pos[2] - pos[0]);
                    std::vector<std::vector<size_t>> faces = {{0, 1, 2}};
                    if (dot(n, u) > 0)
                        faces = {{0, 2, 1}};
                    const double h = std::min(1e-6, 1e-3 * a * s * (1.0 - std::abs(x)));
                    Setup S;
                    S.mesh.reset(new ManifoldSurfaceMesh(faces));
                    S.geometry.reset(new VertexPositionGeometry(*S.mesh));
                    for (Vertex v : S.mesh->vertices())
                        S.geometry->inputVertexPositions[v] = pos[v.getIndex()];
                    S.geometry->refreshQuantities();
                    S.attach(W, a, s, p);
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
                    // x = +-0.9999 is within one FD step of the edge on one side for large steps; fine
                    const double scale = std::max(fmax, 1e-8);
                    if (trial == 0)
                        std::printf("  %-4g %-8g %-12.5e %-11.4e %-18.3e %.3e\n", p, x, E0, fmax, worstM / scale, worstB / scale);
                    worstAll = std::max(worstAll, std::max(worstM, worstB) / scale);
                    if (!(E0 < 0.0) && std::abs(x) < 1.0)
                    {
                        ok &= check(false, "face was not selected / zero energy inside the shell");
                    }
                }
            }
        }
        std::printf("worst relative error over all single-triangle checks: %.3e\n", worstAll);
        ok &= check(worstAll < 1e-6, "single triangles: force matches finite differences (rel. < 1e-6)");

        // outside the shell: exactly zero
        {
            Vector3 u = unit(Vector3{0.2, 0.3, 0.9});
            std::vector<std::vector<size_t>> faces = {{0, 2, 1}};
            Setup S;
            S.mesh.reset(new ManifoldSurfaceMesh(faces));
            S.geometry.reset(new VertexPositionGeometry(*S.mesh));
            Vector3 t1 = unit(cross(u, Vector3{1, 0, 0})), t2 = cross(u, t1);
            Vector3 c = 1.26 * u;
            std::vector<Vector3> pos = {c + 0.2 * t1, c - 0.2 * t1 + 0.2 * t2, c - 0.2 * t1 - 0.2 * t2};
            for (Vertex v : S.mesh->vertices())
                S.geometry->inputVertexPositions[v] = pos[v.getIndex()];
            S.geometry->refreshQuantities();
            S.attach(W, a, s, 4.0);
            ok &= check(S.adhesion->Face_Energy(S.mesh->face(0)) == 0.0, "no energy outside the shell (r = 1.26)");
        }

        // almost edge-on (face-selection boundary): the energy goes to zero continuously with omega
        {
            std::printf("\nface-selection boundary: triangle at x = 0.2 tilted so that dot(radial, normal) -> 0-\n");
            for (double p : {1.0, 4.0})
            {
                double worstAll2 = 0.0;
                for (double tilt : {1e-1, 1e-2, 1e-3, 1e-4})
                {
                    Vector3 u = unit(Vector3{0.2, 0.3, 0.9});
                    Vector3 t1 = unit(cross(u, Vector3{1, 0, 0})), t2 = cross(u, t1);
                    // normal tilted away from -u by angle (pi/2 - tilt): dot(u, n) = -sin(tilt)
                    Vector3 nrm = -std::sin(tilt) * u + std::cos(tilt) * t1;
                    Vector3 e1 = unit(cross(nrm, u));
                    Vector3 e2 = cross(nrm, e1);
                    Vector3 c = a * (1.0 + s * 0.2) * u;
                    std::vector<Vector3> pos = {c + 0.2 * e1, c - 0.15 * e1 + 0.2 * e2, c - 0.15 * e1 - 0.2 * e2};
                    Vector3 n = cross(pos[1] - pos[0], pos[2] - pos[0]);
                    std::vector<std::vector<size_t>> faces = {{0, 1, 2}};
                    if (dot(n, nrm) < 0)
                        faces = {{0, 2, 1}};
                    Setup S;
                    S.mesh.reset(new ManifoldSurfaceMesh(faces));
                    S.geometry.reset(new VertexPositionGeometry(*S.mesh));
                    for (Vertex v : S.mesh->vertices())
                        S.geometry->inputVertexPositions[v] = pos[v.getIndex()];
                    S.geometry->refreshQuantities();
                    S.attach(W, a, s, p);
                    Face f = S.mesh->face(0);
                    const double E0 = S.adhesion->Face_Energy(f);
                    VertexData<Vector3> F(*S.mesh, Vector3{0, 0, 0});
                    Vector3 Fb{0, 0, 0};
                    S.adhesion->Add_Face_Force(f, F, Fb);
                    double worst = 0.0, fmax = 0.0;
                    const double hh = 1e-8;
                    for (Vertex v : S.mesh->vertices())
                        for (int k = 0; k < 3; k++)
                        {
                            Vector3 &q = S.geometry->inputVertexPositions[v];
                            const double old = q[k];
                            q[k] = old + hh;
                            double Ep = S.adhesion->Face_Energy(f);
                            q[k] = old - hh;
                            double Em = S.adhesion->Face_Energy(f);
                            q[k] = old;
                            worst = std::max(worst, std::abs(-(Ep - Em) / (2 * hh) - F[v][k]));
                            fmax = std::max(fmax, std::abs(F[v][k]));
                        }
                    std::printf("  p=%g tilt=%.0e  dot=%.2e  E=%.4e  max|F|=%.3e  max|FD-F|=%.3e\n", p, tilt,
                                dot(S.geometry->faceNormal(f) / norm(S.geometry->faceNormal(f)), unit(pos[0] + pos[1] + pos[2])),
                                E0, fmax, worst);
                    worstAll2 = std::max(worstAll2, worst / std::max(fmax, 1e-12));
                }
                char msg[100];
                std::snprintf(msg, sizeof msg, "p=%g near-edge-on triangles: FD matches (worst rel %.2e)", p, worstAll2);
                ok &= check(worstAll2 < 1e-5, msg);
            }
        }
    }

    // ---- 3. a sphere patch -----------------------------------------------------------------------
    {
        const double a = 1.0, s = 0.6, W = 1.0; // wide shell so that many faces are inside
        std::printf("\nsphere patch (radius 3 sphere, bead outside), many faces at different x\n");
        for (double p : {1.0, 2.0, 4.0})
        {
            Setup S;
            std::tie(S.mesh, S.geometry) = readManifoldSurfaceMesh(path);
            for (Vertex v : S.mesh->vertices())
                S.geometry->inputVertexPositions[v] = 3.0 * S.geometry->inputVertexPositions[v] +
                                                      0.04 * Vector3{gauss(gen), gauss(gen), gauss(gen)};
            S.geometry->refreshQuantities();
            S.attach(W, a, s, p);
            // sphere top at z = 3: put the bead so that the pole is at r = a and the neighbours fan out
            S.bead.Pos = Vector3{0.15, -0.1, 3.0 + 1.0 * a};
            double E0 = S.adhesion->Tot_Energy();
            int nsel = 0;
            double xmin = 9, xmax = -9;
            for (Face f : S.mesh->faces())
                if (S.adhesion->Face_Energy(f) != 0.0)
                {
                    nsel++;
                    Halfedge he = f.halfedge();
                    Vector3 c = (S.geometry->inputVertexPositions[he.vertex()] + S.geometry->inputVertexPositions[he.next().vertex()] +
                                 S.geometry->inputVertexPositions[he.next().next().vertex()]) /
                                3.0;
                    double x = (norm(c - S.bead.Pos) - a) / (s * a);
                    xmin = std::min(xmin, x);
                    xmax = std::max(xmax, x);
                }
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
            std::printf("  p=%g: E = %.8g, %d faces selected, x in [%.3f, %.3f], max|F| %.3e, rel.err membrane %.2e, bead %.2e\n", p, E0,
                        nsel, xmin, xmax, fmax, worst / fmax, wb / fmax);
            char msg[100];
            std::snprintf(msg, sizeof msg, "p=%g patch: forces match finite differences", p);
            ok &= check(nsel > 5 && worst < 1e-6 * fmax && wb < 1e-6 * fmax, msg);
        }
    }

    std::printf(ok ? "PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
