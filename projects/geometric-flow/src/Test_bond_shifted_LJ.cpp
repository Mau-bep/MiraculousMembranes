// Checks for the "Shifted_LJ" bond between two beads (soft repulsion: LJ cut at
// its minimum rc = 2^(1/6) sigma and shifted up by epsilon):
//   1. the bond energy has the expected values: epsilon at r = sigma, zero at
//      and beyond rc, positive and decreasing inside
//   2. the bond force against central finite differences of the bond energy,
//      for several distances (inside, close to rc, outside) and directions
//   3. the bond Hessian against central finite differences of the bond force
//   4. Newton's third law and symmetry of the Hessian
//
// Usage: Test_bond_shifted_LJ [results.txt]   (default Bond_shifted_LJ_fd.txt)
// The same lines go to the screen and to the results file.

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "Beads.h"
#include "Interaction.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace
{
    FILE *results = nullptr;

    void out(const char *fmt, ...)
    {
        va_list args;
        va_start(args, fmt);
        std::vfprintf(stdout, fmt, args);
        va_end(args);
        if (results)
        {
            va_start(args, fmt);
            std::vfprintf(results, fmt, args);
            va_end(args);
        }
    }

    bool check(bool ok, const char *what)
    {
        out("%s: %s\n", ok ? "ok  " : "FAIL", what);
        return ok;
    }

    // Two beads bonded to each other (both list the bond), on a one-triangle mesh
    struct Pair
    {
        std::unique_ptr<ManifoldSurfaceMesh> mesh;
        std::unique_ptr<VertexPositionGeometry> geometry;
        No_mem_Inter inter[2];
        Bead bead[2];

        Pair(double epsilon, double sigma)
        {
            mesh.reset(new ManifoldSurfaceMesh(std::vector<std::vector<size_t>>{{0, 1, 2}}));
            geometry.reset(new VertexPositionGeometry(*mesh));
            for (int i = 0; i < 2; i++)
            {
                inter[i].mesh = mesh.get();
                inter[i].geometry = geometry.get();
                inter[i].Bead_1 = &bead[i];
                bead[i].mesh = mesh.get();
                bead[i].geometry = geometry.get();
                bead[i].Bead_I = &inter[i];
                bead[i].Bead_id = i;
                bead[i].Total_beads = 2;
            }
            bead[0].Add_bead(&bead[1], "Shifted_LJ", {epsilon, sigma});
            bead[1].Add_bead(&bead[0], "Shifted_LJ", {epsilon, sigma});
        }

        double energy() { return inter[0].Bond_energy() + inter[1].Bond_energy(); }

        // bead 1 sits at the origin, bead 2 at r * dir
        void place(double r, Vector3 dir)
        {
            bead[0].Pos = Vector3{0.0, 0.0, 0.0};
            bead[1].Pos = r * dir;
        }

        Eigen::Matrix<double, 6, 1> force()
        {
            Vector3 f0 = inter[0].Bond_force(), f1 = inter[1].Bond_force();
            Eigen::Matrix<double, 6, 1> f;
            f << f0.x, f0.y, f0.z, f1.x, f1.y, f1.z;
            return f;
        }

        Eigen::Matrix<double, 6, 6> hessian()
        {
            Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();
            const int offset = 3 * int(mesh->nVertices());
            for (int i = 0; i < 2; i++)
                for (const Eigen::Triplet<double> &t : inter[i].Hessian_bonds_triplet())
                    H(t.row() - offset, t.col() - offset) += t.value();
            return H;
        }

        double &coordinate(int c) { return bead[c / 3].Pos[c % 3]; }
    };
} // namespace

int main(int argc, char **argv)
{
    const char *path = argc > 1 ? argv[1] : "Bond_shifted_LJ_fd.txt";
    results = std::fopen(path, "w");
    if (!results)
        std::fprintf(stderr, "Cannot write %s, printing to the screen only\n", path);

    const double epsilon = 1.7, sigma = 0.8;
    const double rc = std::pow(2.0, 1.0 / 6.0) * sigma;
    const double h = 1e-6;
    bool ok = true;

    Pair pair(epsilon, sigma);
    out("Shifted_LJ bond: epsilon = %g, sigma = %g, rc = 2^(1/6) sigma = %.12g, FD step h = %g\n", epsilon, sigma, rc, h);
    out("Both beads list the bond: E = sum of both listings, force on each bead = -dE/dx\n\n");

    // 1. Energy values
    {
        const Vector3 x{1.0, 0.0, 0.0};
        pair.place(sigma, x);
        double E_sigma = pair.energy();
        pair.place(rc, x);
        double E_rc = pair.energy();
        double F_rc = norm(pair.inter[0].Bond_force());
        pair.place(rc * (1.0 + 1e-9), x);
        double E_out = pair.energy();
        pair.place(2.5 * sigma, x);
        double E_far = pair.energy();
        double F_far = norm(pair.inter[0].Bond_force());
        out("E(sigma) = %.15g (epsilon = %g)\n", E_sigma, epsilon);
        out("E(rc) = %.3e, |F|(rc) = %.3e, E(rc(1+1e-9)) = %.3e, E(2.5 sigma) = %.3e, |F|(2.5 sigma) = %.3e\n",
            E_rc, F_rc, E_out, E_far, F_far);
        ok &= check(std::fabs(E_sigma - epsilon) < 1e-12, "E(sigma) = epsilon");
        ok &= check(std::fabs(E_rc) < 1e-12 && F_rc < 1e-10, "E and force vanish at the cutoff");
        ok &= check(E_out == 0.0 && E_far == 0.0 && F_far == 0.0, "no energy or force beyond the cutoff");

        bool monotone = true, positive = true;
        double previous = 1e300;
        for (double r = 0.6 * sigma; r < rc; r += 1e-3)
        {
            pair.place(r, x);
            double E = pair.energy();
            monotone &= E < previous;
            positive &= E > 0.0;
            previous = E;
        }
        ok &= check(monotone && positive, "E > 0 and decreasing for r < rc (purely repulsive)");
    }

    // 2./3. Force and Hessian against finite differences
    out("\n%8s %12s %14s %14s %14s %14s %14s\n", "r/sigma", "E", "|F|", "max|F-FD(E)|", "max|dF|/max", "max|H-FD(F)|", "max|H-H^T|");
    std::mt19937 rng(7);
    std::normal_distribution<double> gauss(0.0, 1.0);
    const double distances[] = {0.85, 0.9, 0.95, 1.0, 1.05, 1.10, 1.115, 1.12, 1.13, 1.3, 2.0};
    double worst_force = 0.0, worst_hess = 0.0, worst_sym = 0.0, worst_newton = 0.0;
    for (double ratio : distances)
    {
        const double r = ratio * sigma;
        // Skip a step that would straddle the cutoff, where E'' jumps
        const bool straddles = std::fabs(r - rc) < 2 * h;
        double row_force = 0.0, row_hess = 0.0, row_sym = 0.0, row_F = 0.0, row_E = 0.0;
        for (int d = 0; d < 3; d++)
        {
            Vector3 dir = unit(Vector3{gauss(rng), gauss(rng), gauss(rng)});
            pair.place(r, dir);
            row_E = pair.energy();
            Eigen::Matrix<double, 6, 1> F = pair.force();
            Eigen::Matrix<double, 6, 6> H = pair.hessian();
            row_F = std::max(row_F, F.head<3>().norm());
            worst_newton = std::max(worst_newton, (F.head<3>() + F.tail<3>()).norm());
            row_sym = std::max(row_sym, (H - H.transpose()).cwiseAbs().maxCoeff());

            if (straddles)
                continue;
            Eigen::Matrix<double, 6, 6> Hfd;
            for (int c = 0; c < 6; c++)
            {
                const double x0 = pair.coordinate(c);
                pair.coordinate(c) = x0 + h;
                double Ep = pair.energy();
                Eigen::Matrix<double, 6, 1> Fp = pair.force();
                pair.coordinate(c) = x0 - h;
                double Em = pair.energy();
                Eigen::Matrix<double, 6, 1> Fm = pair.force();
                pair.coordinate(c) = x0;
                row_force = std::max(row_force, std::fabs(-(Ep - Em) / (2 * h) - F(c)));
                Hfd.col(c) = -(Fp - Fm) / (2 * h);
            }
            row_hess = std::max(row_hess, (H - Hfd).cwiseAbs().maxCoeff());
        }
        const double scale_F = std::max(1.0, row_F);
        const double scale_H = std::max(1.0, row_F / r);
        out("%8.4f %12.6g %14.6e %14.3e %14.3e %14.3e %14.3e%s\n", ratio, row_E, row_F, row_force, row_force / scale_F,
            row_hess, row_sym, straddles ? "   (FD skipped: within 2h of rc)" : "");
        worst_force = std::max(worst_force, row_force / scale_F);
        worst_hess = std::max(worst_hess, row_hess / scale_H);
        worst_sym = std::max(worst_sym, row_sym);
    }
    out("\nworst relative force error  %.3e\nworst relative Hessian error %.3e\n"
        "worst Hessian asymmetry      %.3e\nworst |F1 + F2|              %.3e\n",
        worst_force, worst_hess, worst_sym, worst_newton);
    ok &= check(worst_force < 1e-6, "bond force matches finite differences of the bond energy");
    ok &= check(worst_hess < 1e-5, "bond Hessian matches finite differences of the bond force");
    ok &= check(worst_sym < 1e-9, "bond Hessian is symmetric");
    ok &= check(worst_newton < 1e-12, "the two beads feel opposite forces");

    out(ok ? "PASS\n" : "FAILED\n");
    if (results)
        std::fclose(results);
    return ok ? 0 : 1;
}
