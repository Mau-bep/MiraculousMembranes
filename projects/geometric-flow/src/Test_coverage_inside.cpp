// Checks for the Coverage energy of a bead INSIDE the closed membrane ("outside": -1 of the bead, Bead::side = -1),
//   1. the input file loads (it used to be refused) and sets Bead::side; the default is +1
//   2. the faces it counts: a bead of radius 1 inside a sphere of radius 3, tangent to the wall, counts faces on the wall
//      (E_Coverage equals K (Omega - target)^2 with Omega summed by bead_geometry::evaluateFaceCoverage with side = -1),
//      the same bead with side +1 counts none, and a bead outside the sphere counts none with side -1
//   3. mirror: side -1 on the outward sphere is side +1 on the sphere with the faces reversed: the same energy and the same
//      forces on the vertices and on the bead
//   4. the forces are -dE/dx: central finite differences of the energy against F_Coverage, on the vertices with the largest
//      force and on the bead, for a bead inside and (control) outside, on a perturbed mesh
//   5. Debug_Coverage runs for a bead inside
//
// Usage: Test_coverage_inside [mesh.obj [workdir]]   (default ../../../input/Simplesphere.obj, the current folder; the
//        temporary input files are written to workdir and removed)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "geometrycentral/surface/meshio.h"

#include "BeadGeometry.h"
#include "SimConfig.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace
{
    bool check(bool ok, const char *what)
    {
        std::printf("  [%s] %s\n", ok ? "ok" : "FAIL", what);
        return ok;
    }

    const double K_COVERAGE = 600.0;
    const double TARGET = 0.5;

    // The same xorshift64* as Test_exactness, so the perturbation is the same on every machine
    struct Noise
    {
        uint64_t state;
        double next() // uniform in [-1, 1)
        {
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            return double((state * 2685821657736338717ULL) >> 11) * (2.0 / 9007199254740992.0) - 1.0;
        }
        Vector3 vec() { return Vector3{next(), next(), next()}; }
    };

    std::string replace_all(std::string text, const std::string &from, const std::string &to)
    {
        for (size_t k = text.find(from); k != std::string::npos; k = text.find(from, k + to.size()))
            text.replace(k, from.size(), to);
        return text;
    }

    // A sphere (rescale 3: radius 3) and one bead of radius 1 with only the Coverage energy
    void write_config(const std::string &path, const std::string &mesh, int outside, const std::string &pos)
    {
        std::string text = R"({
    "init_file": "%MESH%", "first_dir": "OVERRIDDEN/", "save_interval": 50, "timesteps": 10,
    "remeshing": false, "remesh_every": 5, "adapt_remesh": false, "Count_remesh": false,
    "Integration": "Gradient_descent", "BFGS_saved_states": 15, "backtrack": true, "Timestep": 1.0,
    "boundary": false, "recentering": false, "saving_states": false, "momentum": false, "rescale": 3.0,
    "Switches": ["BFGS"], "Switch_times": [150],
    "Energies": [ {"Name": "Coverage", "constants": [%K%, %TARGET%]} ],
    "Beads": [ {"gradient_order": "Bead", "radius": 1.0, "rc": 1.25, "state": "froze", "outside": %OUTSIDE%,
                "Pos": [%POS%], "mem_inter": "None", "inter_str": 0.0, "bonds": [], "bonds_constants": [], "Beads": []} ],
    "remesher": {"aspect_min": 0.2, "refine_angle": 0.6, "size_max": 0.45, "size_min": 0.05}
})";
        text = replace_all(text, "%MESH%", mesh);
        text = replace_all(text, "%K%", std::to_string(K_COVERAGE));
        text = replace_all(text, "%TARGET%", std::to_string(TARGET));
        text = replace_all(text, "%OUTSIDE%", std::to_string(outside));
        text = replace_all(text, "%POS%", pos);
        std::ofstream(path.c_str()) << text;
    }

    // The same mesh with every face turned around (the vertices keep their order)
    void write_reversed_mesh(const std::string &from, const std::string &to)
    {
        std::ifstream in(from.c_str());
        std::ofstream out(to.c_str());
        std::string line;
        while (std::getline(in, line))
        {
            if (line.compare(0, 2, "f ") == 0)
            {
                std::istringstream words(line.substr(2));
                std::vector<std::string> corner;
                std::string word;
                while (words >> word)
                    corner.push_back(word);
                if (corner.size() == 3)
                    line = "f " + corner[0] + " " + corner[2] + " " + corner[1];
            }
            out << line << "\n";
        }
    }

    std::unique_ptr<Simulation> make(const std::string &config)
    {
        std::unique_ptr<Simulation> sim(new Simulation());
        build_simulation(load_config(config), *sim);
        sim->Sim_handler.update_vertex_normals();
        return sim;
    }

    double energy(Simulation &sim)
    {
        double E = 0.0;
        sim.Sim_handler.Calculate_energies(&E);
        return E;
    }

    // Omega of bead i, summed by the face selection of bead_geometry with the side of the bead
    double omega(Simulation &sim, size_t i, int &selected)
    {
        const Bead &bead = sim.Beads[i];
        double sum = 0.0;
        selected = 0;
        for (Face f : sim.mesh->faces())
        {
            const bead_geometry::FaceCoverageData d =
                bead_geometry::evaluateFaceCoverage(*sim.geometry, f, bead.Pos, bead.sigma, 0.25, 1.0, bead.side);
            if (!d.selected)
                continue;
            selected++;
            sum += d.weight * d.omegaUnsigned / (4.0 * bead_geometry::PI_VALUE);
        }
        return sum;
    }

    double relative(double a, double b)
    {
        return std::abs(a - b) / std::max(std::abs(a), std::abs(b));
    }

    // Energy and forces of the current state, the vertex forces in `forces` and the force of bead 0 in `bead`
    double evaluate(Simulation &sim, std::vector<Vector3> &forces, Vector3 &bead)
    {
        const double E = energy(sim);
        sim.Sim_handler.Calculate_gradient();
        forces.assign(sim.mesh->nVertices(), Vector3{0.0, 0.0, 0.0});
        for (Vertex v : sim.mesh->vertices())
            forces[v.getIndex()] = sim.Sim_handler.Current_grad[v];
        bead = sim.Beads[0].CoverageForce;
        return E;
    }

    // Forces against central differences of the energy; returns the worst error relative to the largest force
    double finite_difference_error(Simulation &sim, double &fmax, size_t &checked)
    {
        std::vector<Vector3> forces;
        Vector3 beadForce;
        evaluate(sim, forces, beadForce);

        std::vector<std::pair<double, size_t>> order;
        for (size_t i = 0; i < forces.size(); i++)
            order.push_back(std::make_pair(norm(forces[i]), i));
        std::sort(order.rbegin(), order.rend());
        fmax = std::max(order[0].first, norm(beadForce));

        const double h = 1e-6;
        double worst = 0.0;
        checked = 0;
        const size_t nVertices = std::min<size_t>(10, order.size());
        for (size_t k = 0; k < nVertices; k++)
        {
            const Vertex v = sim.mesh->vertex(order[k].second);
            const Vector3 x0 = sim.geometry->inputVertexPositions[v];
            for (int axis = 0; axis < 3; axis++)
            {
                Vector3 step{0.0, 0.0, 0.0};
                step[axis] = h;
                sim.geometry->inputVertexPositions[v] = x0 + step;
                const double Ep = energy(sim);
                sim.geometry->inputVertexPositions[v] = x0 - step;
                const double Em = energy(sim);
                sim.geometry->inputVertexPositions[v] = x0;
                const double fd = -(Ep - Em) / (2.0 * h);
                worst = std::max(worst, std::abs(fd - forces[order[k].second][axis]));
                checked++;
            }
        }
        const Vector3 p0 = sim.Beads[0].Pos;
        for (int axis = 0; axis < 3; axis++)
        {
            Vector3 step{0.0, 0.0, 0.0};
            step[axis] = h;
            sim.Beads[0].Pos = p0 + step;
            const double Ep = energy(sim);
            sim.Beads[0].Pos = p0 - step;
            const double Em = energy(sim);
            sim.Beads[0].Pos = p0;
            worst = std::max(worst, std::abs(-(Ep - Em) / (2.0 * h) - beadForce[axis]));
            checked++;
        }
        return worst / fmax;
    }

    void perturb(Simulation &sim, uint64_t seed)
    {
        Noise noise{seed};
        const double amplitude = 0.02 * sim.geometry->meanEdgeLength();
        for (Vertex v : sim.mesh->vertices())
            sim.geometry->inputVertexPositions[v] += amplitude * noise.vec();
        sim.Beads[0].Pos += amplitude * noise.vec();
    }
} // namespace

int main(int argc, char **argv)
{
    const std::string mesh = argc > 1 ? argv[1] : "../../../input/Simplesphere.obj";
    const std::string dir = argc > 2 ? std::string(argv[2]) + "/" : "";
    const std::string reversed_mesh = dir + "Test_coverage_inside_reversed.obj";
    const std::string inside_cfg = dir + "Test_coverage_inside_in.json";
    const std::string outside_cfg = dir + "Test_coverage_inside_out.json";
    const std::string mirror_cfg = dir + "Test_coverage_inside_mirror.json";
    const std::string far_cfg = dir + "Test_coverage_inside_far.json";
    const std::string control_cfg = dir + "Test_coverage_inside_control.json";

    bool ok = true;
    try
    {
        write_reversed_mesh(mesh, reversed_mesh);
        // inside: the wall is at distance 3 from the centre, the bead of radius 1 is tangent to it from the inside
        write_config(inside_cfg, mesh, -1, "2.0, 0.0, 0.0");
        // the same bead and "outside": 1 (the old convention), and a bead outside the sphere with "outside": -1
        write_config(outside_cfg, mesh, 1, "2.0, 0.0, 0.0");
        write_config(mirror_cfg, reversed_mesh, 1, "2.0, 0.0, 0.0");
        write_config(far_cfg, mesh, -1, "4.0, 0.0, 0.0");

        std::cout << "1. the input file\n";
        {
            std::unique_ptr<Simulation> in = make(inside_cfg);
            std::unique_ptr<Simulation> out = make(outside_cfg);
            ok &= check(in->Beads[0].side == -1.0, "\"outside\": -1 with the Coverage energy loads and sets side = -1");
            ok &= check(out->Beads[0].side == 1.0, "\"outside\": 1 gives side = +1");
            ok &= check(Bead().side == 1.0, "a default Bead has side = +1");
        }

        std::cout << "2. the faces it counts\n";
        {
            std::unique_ptr<Simulation> in = make(inside_cfg);
            std::unique_ptr<Simulation> out = make(outside_cfg);
            std::unique_ptr<Simulation> far = make(far_cfg);
            int nIn = 0, nOut = 0, nFar = 0;
            const double omegaIn = omega(*in, 0, nIn);
            const double omegaOut = omega(*out, 0, nOut);
            const double omegaFar = omega(*far, 0, nFar);
            std::printf("  inside: %d faces, Omega = %.6f | same bead, side +1: %d faces | bead outside, side -1: %d faces\n", nIn, omegaIn, nOut, nFar);
            ok &= check(nIn > 0 && omegaIn > 0.0, "the bead inside counts faces of the wall");
            ok &= check(nOut == 0 && omegaOut == 0.0, "the same bead with side +1 counts none (the normals point away from it)");
            ok &= check(nFar == 0 && omegaFar == 0.0, "a bead outside the sphere with side -1 counts none");
            const double Ein = energy(*in), Eout = energy(*out), Efar = energy(*far);
            const double expectedIn = K_COVERAGE * (omegaIn - TARGET) * (omegaIn - TARGET);
            ok &= check(relative(Ein, expectedIn) < 1e-12, "E_Coverage of the bead inside = K (Omega - target)^2 with Omega of side -1");
            ok &= check(Eout == K_COVERAGE * TARGET * TARGET && Efar == K_COVERAGE * TARGET * TARGET, "the beads that count nothing have E = K target^2");
            ok &= check(Ein < Eout, "the bead inside lowers the energy (it covers part of the target)");
        }

        std::cout << "3. mirror: side -1 on the sphere = side +1 on the reversed sphere\n";
        {
            std::unique_ptr<Simulation> a = make(inside_cfg);
            std::unique_ptr<Simulation> b = make(mirror_cfg);
            perturb(*a, 0x9E3779B97F4A7C15ULL);
            perturb(*b, 0x9E3779B97F4A7C15ULL);
            std::vector<Vector3> fa, fb;
            Vector3 ba, bb;
            const double Ea = evaluate(*a, fa, ba);
            const double Eb = evaluate(*b, fb, bb);
            double fmax = 0.0, diff = 0.0;
            for (size_t i = 0; i < fa.size(); i++)
            {
                fmax = std::max(fmax, norm(fa[i]));
                diff = std::max(diff, norm(fa[i] - fb[i]));
            }
            std::printf("  E = %.10f vs %.10f, largest force %.4g, largest difference %.3g, bead force difference %.3g\n", Ea, Eb, fmax, diff, norm(ba - bb));
            ok &= check(Ea < K_COVERAGE * TARGET * TARGET, "the perturbed bead inside still counts faces");
            ok &= check(relative(Ea, Eb) < 1e-12, "the same energy");
            ok &= check(diff <= 1e-10 * fmax && norm(ba - bb) <= 1e-10 * fmax, "the same forces on the vertices and on the bead");
        }

        std::cout << "4. the forces against finite differences of the energy\n";
        {
            std::unique_ptr<Simulation> in = make(inside_cfg);
            // the control: the bead outside the sphere with the usual convention, same sphere
            write_config(control_cfg, mesh, 1, "4.0, 0.0, 0.0");
            std::unique_ptr<Simulation> ctrl = make(control_cfg);
            perturb(*in, 0x1234567ULL);
            perturb(*ctrl, 0x7654321ULL);
            double fmaxIn = 0.0, fmaxCtrl = 0.0;
            size_t nIn = 0, nCtrl = 0;
            const double errIn = finite_difference_error(*in, fmaxIn, nIn);
            const double errCtrl = finite_difference_error(*ctrl, fmaxCtrl, nCtrl);
            std::printf("  bead inside: largest force %.4g, worst error %.3g of it (%zu derivatives) | bead outside: %.4g, %.3g (%zu)\n",
                        fmaxIn, errIn, nIn, fmaxCtrl, errCtrl, nCtrl);
            ok &= check(fmaxIn > 0.0 && errIn < 1e-6, "bead inside: F_Coverage matches finite differences (rel. < 1e-6)");
            ok &= check(fmaxCtrl > 0.0 && errCtrl < 1e-6, "bead outside (control): the same");
        }

        std::cout << "5. Debug_Coverage\n";
        {
            std::unique_ptr<Simulation> in = make(inside_cfg);
            in->Sim_handler.Debug_Coverage(in->Sim_handler.Energy_constants[0]);
            ok &= check(true, "Debug_Coverage ran for a bead inside");
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        ok = false;
    }

    std::remove(reversed_mesh.c_str());
    std::remove(inside_cfg.c_str());
    std::remove(outside_cfg.c_str());
    std::remove(mirror_cfg.c_str());
    std::remove(far_cfg.c_str());
    std::remove(control_cfg.c_str());
    std::cout << (ok ? "\nAll checks passed\n" : "\nSome checks FAILED\n");
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
