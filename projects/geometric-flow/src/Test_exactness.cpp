// Exactness harness for the energy and force code:
//   Test_exactness <dump_file> <input.json> [<input.json> ...]
//
// Builds each simulation with load_config/build_simulation, moves the vertices
// and beads to a few fixed pseudo-random perturbations of the initial state,
// then once more after a remesh, and writes Calculate_energies and
// Calculate_gradient at every state in hexfloat. A refactor that keeps the
// results bit-identical gives a byte-identical dump (compare with cmp).

#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "geometrycentral/surface/remeshing.h"

#include "SimConfig.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace
{
    // Portable generator (xorshift64*), so the states are the same on every machine
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

    void dump_evaluation(std::ofstream &out, Simulation &sim, const std::string &label)
    {
        E_Handler &handler = sim.Sim_handler;
        double E = 0.0;
        handler.Calculate_energies(&E);
        handler.Calculate_gradient();

        out << "state " << label << " vertices " << sim.mesh->nVertices() << " faces " << sim.mesh->nFaces() << "\n";
        out << "total_energy " << E << "\nenergies";
        for (size_t i = 0; i < handler.Energies.size(); i++)
            out << " " << handler.Energies[i] << " " << handler.Energy_values[i];
        out << "\ngradient_norms";
        for (double g : handler.Gradient_norms)
            out << " " << g;
        out << "\n";
        for (size_t i = 0; i < sim.Beads.size(); i++)
        {
            const Bead &b = sim.Beads[i];
            out << "bead " << i << " force " << b.Total_force.x << " " << b.Total_force.y << " " << b.Total_force.z
                << " coverage_force " << b.CoverageForce.x << " " << b.CoverageForce.y << " " << b.CoverageForce.z << "\n";
        }
        for (Vertex v : sim.mesh->vertices())
        {
            const Vector3 &f = handler.Current_grad[v];
            out << f.x << " " << f.y << " " << f.z << "\n";
        }
    }

    void run_config(std::ofstream &out, const std::string &path)
    {
        Simulation sim;
        build_simulation(load_config(path), sim);
        if (!sim.M3DG.boundary)
            sim.Sim_handler.Constraints.push_back("Volume_constraint");
        sim.Sim_handler.update_vertex_normals();

        out << std::hexfloat << "config " << path << "\n";
        std::vector<Vector3> x0(sim.mesh->nVertices());
        for (Vertex v : sim.mesh->vertices())
            x0[v.getIndex()] = sim.geometry->inputVertexPositions[v];
        std::vector<Vector3> b0;
        for (const Bead &b : sim.Beads)
            b0.push_back(b.Pos);

        const double amplitude = 0.02 * sim.geometry->meanEdgeLength();
        Noise noise{0x9E3779B97F4A7C15ULL};
        const int n_states = 4;
        for (int k = 0; k < n_states; k++)
        {
            // State 0 is the input mesh itself
            const double a = k == 0 ? 0.0 : amplitude;
            for (Vertex v : sim.mesh->vertices())
                sim.geometry->inputVertexPositions[v] = x0[v.getIndex()] + a * noise.vec();
            for (size_t i = 0; i < sim.Beads.size(); i++)
                sim.Beads[i].Pos = b0[i] + a * noise.vec();
            dump_evaluation(out, sim, "perturbed_" + std::to_string(k));
        }

        // The last perturbed state, remeshed as main_cluster does. The refresh
        // first, because the remesher starts from the cached vertexPositions
        if (sim.cfg.remeshing)
        {
            sim.geometry->refreshQuantities();
            remesh(*sim.mesh, *sim.geometry, sim.Options);
            sim.geometry->refreshQuantities();
            sim.Sim_handler.update_vertex_normals();
            dump_evaluation(out, sim, "remeshed");
        }
    }
} // namespace

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <dump_file> <input.json> [<input.json> ...]\n";
        return EXIT_FAILURE;
    }
    std::ofstream out(argv[1]);
    for (int i = 2; i < argc; i++)
    {
        try
        {
            run_config(out, argv[i]);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error with " << argv[i] << ": " << e.what() << "\n";
            return EXIT_FAILURE;
        }
    }
    std::cout << "Wrote " << argv[1] << "\n";
    return EXIT_SUCCESS;
}
