// Checks for the Bilayer energy, E = K/2 (da - da0)^2 with the reduced area difference da = M / (4 pi R0) of the
// bilayer-coupling model (Seifert, Berndl, Lipowsky, PRA 44, 1182),
//   1. da of a sphere is 1, and it does not depend on the size of the vesicle
//   2. da of spheroids (the sphere stretched along z) against the quadrature of int H dA over the surface of revolution,
//      and the sign: the sphere with the faces reversed has da = -1
//   3. the input file: the constants K and da0 are read (with the config completed to two constants), the energy is
//      exactly K/2 (da - da0)^2, and K = 0 switches the term off
//   4. the forces are -dE/dx: central finite differences of the energy against F_Bilayer on a perturbed sphere and on a
//      perturbed spheroid, and Calculate_gradient adds the force of the term
//
// Usage: Test_bilayer [mesh.obj [workdir]]   (default ../../../input/Simplesphere.obj, the current folder; the temporary
//        input files are written to workdir and removed)

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

    const double K_BILAYER = 7.0;
    const double DA_TARGET = 0.8;

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

    // The mesh rescaled by `rescale` with only the Bilayer energy; `constants` is the text after "constants":
    void write_config(const std::string &path, const std::string &mesh, double rescale, const std::string &constants)
    {
        std::string text = R"({
    "init_file": "%MESH%", "first_dir": "OVERRIDDEN/", "save_interval": 50, "timesteps": 10,
    "remeshing": false, "remesh_every": 5, "adapt_remesh": false, "Count_remesh": false,
    "Integration": "Gradient_descent", "BFGS_saved_states": 15, "backtrack": true, "Timestep": 1.0,
    "boundary": false, "recentering": false, "saving_states": false, "momentum": false, "rescale": %RESCALE%,
    "Switches": ["BFGS"], "Switch_times": [150],
    "Energies": [ {"Name": "Bilayer", "constants": %CONSTANTS%} ],
    "Beads": [],
    "remesher": {"aspect_min": 0.2, "refine_angle": 0.6, "size_max": 0.45, "size_min": 0.05}
})";
        text = replace_all(text, "%MESH%", mesh);
        text = replace_all(text, "%RESCALE%", std::to_string(rescale));
        text = replace_all(text, "%CONSTANTS%", constants);
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

    // The mean radius of the vertices, from the centroid
    double mean_radius(Simulation &sim)
    {
        Vector3 c{0.0, 0.0, 0.0};
        for (Vertex v : sim.mesh->vertices())
            c += sim.geometry->inputVertexPositions[v];
        c /= double(sim.mesh->nVertices());
        double r = 0.0;
        for (Vertex v : sim.mesh->vertices())
            r += norm(sim.geometry->inputVertexPositions[v] - c);
        return r / double(sim.mesh->nVertices());
    }

    // The sphere stretched along z by a factor c (its points stay on the spheroid with semi-axes R, R, c R)
    void stretch(Simulation &sim, double c)
    {
        for (Vertex v : sim.mesh->vertices())
            sim.geometry->inputVertexPositions[v].z *= c;
    }

    // da of the spheroid with semi-axes a, a, c by quadrature: the surface r = a sin t, z = c cos t has the principal
    // curvatures k1 = a c / w^3 (meridian) and k2 = c / (a w) (parallel), w = sqrt(a^2 cos^2 t + c^2 sin^2 t), and
    // dA = 2 pi a sin t w dt; M = int (k1 + k2) / 2 dA and da = M / sqrt(4 pi A)
    double spheroid_da(double a, double c)
    {
        const int n = 200000;
        double M = 0.0, A = 0.0;
        for (int i = 0; i < n; i++)
        {
            const double t = (i + 0.5) * M_PI / n;
            const double w = std::sqrt(a * a * std::cos(t) * std::cos(t) + c * c * std::sin(t) * std::sin(t));
            const double dA = 2.0 * M_PI * a * std::sin(t) * w * (M_PI / n);
            M += 0.5 * (a * c / (w * w * w) + c / (a * w)) * dA;
            A += dA;
        }
        return M / std::sqrt(4.0 * M_PI * A);
    }

    double relative(double a, double b)
    {
        return std::abs(a - b) / std::max(std::abs(a), std::abs(b));
    }

    void perturb(Simulation &sim, uint64_t seed)
    {
        Noise noise{seed};
        const double amplitude = 0.02 * sim.geometry->meanEdgeLength();
        for (Vertex v : sim.mesh->vertices())
            sim.geometry->inputVertexPositions[v] += amplitude * noise.vec();
    }

    // The forces of the term against central differences of the energy on the vertices with the largest force; returns
    // the worst error relative to the largest force
    double finite_difference_error(Simulation &sim, const std::vector<double> &constants, double &fmax, size_t &checked)
    {
        const VertexData<Vector3> force = sim.Sim_handler.F_Bilayer(constants);
        std::vector<std::pair<double, size_t>> order;
        for (Vertex v : sim.mesh->vertices())
            order.push_back(std::make_pair(norm(force[v]), v.getIndex()));
        std::sort(order.rbegin(), order.rend());
        fmax = order[0].first;

        const double h = 1e-6;
        double worst = 0.0;
        checked = 0;
        const size_t n = std::min<size_t>(15, order.size());
        for (size_t k = 0; k < n; k++)
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
                worst = std::max(worst, std::abs(fd - force[v][axis]));
                checked++;
            }
        }
        return worst / fmax;
    }
} // namespace

int main(int argc, char **argv)
{
    const std::string mesh = argc > 1 ? argv[1] : "../../../input/Simplesphere.obj";
    const std::string dir = argc > 2 ? std::string(argv[2]) + "/" : "";
    const std::string reversed_mesh = dir + "Test_bilayer_reversed.obj";
    const std::string cfg = dir + "Test_bilayer.json";
    const std::string big_cfg = dir + "Test_bilayer_big.json";
    const std::string reversed_cfg = dir + "Test_bilayer_reversed.json";
    const std::string short_cfg = dir + "Test_bilayer_short.json";
    const std::string off_cfg = dir + "Test_bilayer_off.json";
    const std::string constants = "[" + std::to_string(K_BILAYER) + ", " + std::to_string(DA_TARGET) + "]";

    bool ok = true;
    try
    {
        write_reversed_mesh(mesh, reversed_mesh);
        write_config(cfg, mesh, 1.0, constants);
        write_config(big_cfg, mesh, 3.0, constants);
        write_config(reversed_cfg, reversed_mesh, 1.0, constants);
        write_config(short_cfg, mesh, 1.0, "[" + std::to_string(K_BILAYER) + "]");
        write_config(off_cfg, mesh, 1.0, "[0.0, " + std::to_string(DA_TARGET) + "]");

        std::cout << "1. the sphere\n";
        {
            std::unique_ptr<Simulation> sim = make(cfg);
            std::unique_ptr<Simulation> big = make(big_cfg);
            const double da = sim->Sim_handler.Reduced_area_difference();
            const double da_big = big->Sim_handler.Reduced_area_difference();
            std::printf("  %zu vertices, da = %.6f, rescaled by 3: %.6f\n", sim->mesh->nVertices(), da, da_big);
            ok &= check(std::abs(da - 1.0) < 0.03, "da of the sphere is 1 (within the discretisation error)");
            ok &= check(relative(da, da_big) < 1e-12, "da does not depend on the size of the vesicle");
        }

        std::cout << "2. spheroids and the orientation\n";
        {
            for (double c : {0.6, 1.5})
            {
                std::unique_ptr<Simulation> sim = make(cfg);
                const double R = mean_radius(*sim);
                stretch(*sim, c);
                const double da = sim->Sim_handler.Reduced_area_difference();
                const double exact = spheroid_da(R, c * R);
                std::printf("  stretched by %.1f: da = %.6f, quadrature %.6f\n", c, da, exact);
                ok &= check(std::abs(da - exact) < 0.03, "da of the spheroid agrees with int H dA of the surface");
            }
            std::unique_ptr<Simulation> sim = make(cfg);
            std::unique_ptr<Simulation> rev = make(reversed_cfg);
            const double da = sim->Sim_handler.Reduced_area_difference();
            const double da_rev = rev->Sim_handler.Reduced_area_difference();
            ok &= check(relative(da, -da_rev) < 1e-12, "the sphere with the faces reversed has da = -da");
        }

        std::cout << "3. the input file\n";
        {
            std::unique_ptr<Simulation> sim = make(cfg);
            const double da = sim->Sim_handler.Reduced_area_difference();
            const double expected = 0.5 * K_BILAYER * (da - DA_TARGET) * (da - DA_TARGET);
            ok &= check(sim->Sim_handler.Energy_constants[0][0] == K_BILAYER && sim->Sim_handler.Energy_constants[0][1] == DA_TARGET,
                        "the constants are K and da0");
            ok &= check(relative(energy(*sim), expected) < 1e-12, "E = K/2 (da - da0)^2");
            ok &= check(sim->Sim_handler.Energy_values[0] == energy(*sim), "the value of the term is the whole energy");
            ok &= check(energy(*make(off_cfg)) == 0.0, "K = 0 switches the term off");
            bool refused = false;
            try
            {
                make(short_cfg);
            }
            catch (const std::exception &)
            {
                refused = true;
            }
            ok &= check(refused, "a single constant is refused (da0 is needed)");
        }

        std::cout << "4. forces against finite differences\n";
        for (double c : {1.0, 1.4})
        {
            std::unique_ptr<Simulation> sim = make(cfg);
            stretch(*sim, c);
            perturb(*sim, 0x9E3779B97F4A7C15ULL);
            const std::vector<double> k = sim->Sim_handler.Energy_constants[0];
            double fmax = 0.0;
            size_t checked = 0;
            const double error = finite_difference_error(*sim, k, fmax, checked);
            std::printf("  stretched by %.1f: da = %.4f, largest force %.4e, %zu components, worst error / largest = %.2e\n", c,
                        sim->Sim_handler.Reduced_area_difference(), fmax, checked, error);
            ok &= check(fmax > 1e-3 && error < 1e-6, "F_Bilayer = -dE/dx");

            sim->Sim_handler.Calculate_gradient();
            const VertexData<Vector3> force = sim->Sim_handler.F_Bilayer(k);
            double diff = 0.0;
            for (Vertex v : sim->mesh->vertices())
                diff = std::max(diff, norm(sim->Sim_handler.Current_grad[v] - force[v]));
            ok &= check(diff == 0.0, "Calculate_gradient adds the force of the term");
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        ok = false;
    }

    std::remove(reversed_mesh.c_str());
    std::remove(cfg.c_str());
    std::remove(big_cfg.c_str());
    std::remove(reversed_cfg.c_str());
    std::remove(short_cfg.c_str());
    std::remove(off_cfg.c_str());
    std::cout << (ok ? "\nAll checks passed\n" : "\nSome checks FAILED\n");
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
