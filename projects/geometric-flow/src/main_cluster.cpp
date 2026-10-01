// Batch driver: main_cluster <input.json> <Nsim>
//
// Loads the simulation with load_config/build_simulation (SimConfig.h), then
// runs the time loop: switches -> remeshing -> saving -> one integrator step.

#include <sys/stat.h>

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>

#include "geometrycentral/surface/remeshing.h"

#include <EigenRand/EigenRand>

#include "SimConfig.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace
{

    // Mutable run state that the switches and integrators change over time.
    struct RunState
    {
        std::string Integration;
        bool remesher;
        int remesh_every;
        bool adapt_remesh;
        int save_interval;
        std::unordered_map<std::string, int> Switch_times_map;

        size_t last_remesh = 0;
        int remesh_op = 0;
        int remesh_op_last = 0;
        int trgt_remesh_op = 100;
        double integral_error = 0;

        double time = 0.0;
        double dt_sim = 0.0;
        std::vector<std::string> Constraints; // used by the Newton integrators

        std::string basic_name;
        std::string output_file;
        std::vector<std::string> Bead_filenames;
    };

    int index_of(const std::vector<std::string> &names, const std::string &name)
    {
        for (size_t i = 0; i < names.size(); i++)
            if (names[i] == name)
                return i;
        return -1;
    }

    void copy_file(const std::string &from, const std::string &to)
    {
        std::ifstream src(from, std::ios::binary);
        std::ofstream dst(to, std::ios::binary);
        dst << src.rdbuf();
    }

    void print_mesh_diagnostics(Simulation &sim)
    {
        ManifoldSurfaceMesh *mesh = sim.mesh;
        VertexPositionGeometry *geometry = sim.geometry;
        double avg_dih = 0, max_dih = 0, min_dih = 0.1;
        for (Edge e : mesh->edges())
        {
            double dih = fabs(geometry->dihedralAngle(e.halfedge()));
            avg_dih += dih;
            max_dih = std::max(max_dih, dih);
            min_dih = std::min(min_dih, dih);
        }
        std::cout << "Dihedral angles: avg " << avg_dih / mesh->nEdges() << " min " << min_dih << " max " << max_dih << "\n";

        // The sizing fields are only printed, but computing them registers
        // geometry quantities, so this stays to keep runs bit-identical.
        FaceData<double> F_sizings = sim.M3DG.Face_sizings();
        VertexData<double> Sizings = sim.M3DG.Vert_sizing(F_sizings);

        double min_edge_l = 1e4, max_edge_l = -1, avg_edge_l = 0.0;
        for (Edge e : mesh->edges())
        {
            double edge_l = geometry->edgeLength(e);
            max_edge_l = std::max(max_edge_l, edge_l);
            min_edge_l = std::min(min_edge_l, edge_l);
            avg_edge_l += edge_l;
        }
        std::cout << "Edge lengths: min " << min_edge_l << " max " << max_edge_l << " avg " << avg_edge_l / mesh->nEdges() << "\n";
    }

    // Initial perturbations requested in the input file
    void apply_initial_perturbations(Simulation &sim, const RunState &run)
    {
        ManifoldSurfaceMesh *mesh = sim.mesh;
        VertexPositionGeometry *geometry = sim.geometry;
        if (sim.cfg.has_initial_noise)
        {
            double noise_amp = sim.cfg.initial_noise;
            Eigen::Rand::P8_mt19937_64 urng{42};
            Eigen::VectorXd noise = Eigen::Rand::normal<Eigen::VectorXd>(mesh->nVertices(), 0, urng, noise_amp, noise_amp);
            VertexData<Vector3> V_Normals = sim.Sim_handler.F_Volume(std::vector<double>{1.0});
            for (Vertex v : mesh->vertices())
                geometry->inputVertexPositions[v] = geometry->inputVertexPositions[v] + noise(v.getIndex()) * V_Normals[v];
            Save_mesh(mesh, geometry, run.basic_name, 1);
        }
        if (sim.cfg.harmonic)
        {
            double radius = 2.0;
            double constant = 0.5 * sqrt(3 / 3.1415926535);
            for (Vertex v : mesh->vertices())
            {
                double displacement = fabs(geometry->inputVertexPositions[v].z) / (radius);
                geometry->inputVertexPositions[v] = geometry->inputVertexPositions[v] + displacement * geometry->inputVertexPositions[v] * constant;
            }
        }
    }

    void apply_switches(Simulation &sim, RunState &run, size_t current_t)
    {
        const std::vector<std::string> &Switches = sim.cfg.switches;
        for (size_t sw = 0; sw < Switches.size(); sw++)
        {
            const std::string &Switch = Switches[sw];
            int Switch_t = run.Switch_times_map[Switch];
            if (Switch_t < 0 || current_t != size_t(Switch_t))
                continue;

            if (Switch == "Newton" || Switch == "Newton-Normal")
            {
                run.Integration = Switch;
                run.remesh_every = -1;
            }
            else if (Switch == "BFGS")
            {
                run.Integration = "BFGS";
                sim.M3DG.BFGS_iter = 0;
            }
            else if (Switch == "IpOpt" || Switch == "IpOpt-Normal")
            {
                run.Integration = Switch;
                run.Switch_times_map[Switch] = current_t + 100;
            }
            else if (Switch == "BFGS-Normal")
            {
                run.Integration = "BFGS-Normal";
                sim.M3DG.BFGS_iter = 0;
                run.remesh_every = -1;
            }
            else if (Switch == "Freeze_beads" || Switch == "Free_beads")
            {
                std::cout << "Switching the beads to " << (Switch == "Freeze_beads" ? "froze" : "default") << "\n";
                for (Bead &bead : sim.Beads)
                    bead.state = Switch == "Freeze_beads" ? "froze" : "default";
                run.Switch_times_map[Switch] = -1;
            }
            else if (Switch == "No_remesh")
            {
                run.remesher = false;
                run.Switch_times_map[Switch] = -1;
                for (size_t i = 0; i < sim.Energies.size(); i++)
                {
                    if (sim.Energies[i] == "Edge_reg")
                    {
                        sim.Energy_constants[i][0] = sim.Energy_constants[i][1];
                        sim.Sim_handler.Energy_constants = sim.Energy_constants;
                        std::cout << "Setting edge reg constant to" << sim.Energy_constants[i][1] << "\n";
                    }
                }
            }
            else if (Switch == "Remesh_always")
            {
                run.remesher = true;
                run.adapt_remesh = false;
                run.remesh_every = 1;
                run.Switch_times_map[Switch] = -1;
                std::cout << "Remeshing every step activated\n";
            }
            else if (Switch == "Restore_remeshing")
            {
                run.remesher = true;
                run.remesh_every = sim.cfg.remesh_every;
                run.adapt_remesh = sim.cfg.adapt_remesh;
                run.Switch_times_map[Switch] = -1;
                std::cout << "Restored remeshing frequency to " << run.remesh_every << "\n";
            }
            else if (Switch == "Volume_constraint")
            {
                sim.Energies.push_back("Volume_constraint");
                sim.Energy_constants.push_back({10000, sim.V_bar});
                sim.Sim_handler.Energies = sim.Energies;
                sim.Sim_handler.Energy_constants = sim.Energy_constants;
                std::cout << "Added the volume constraint\n";
                run.Switch_times_map[Switch] = -1;
            }
            else if (Switch == "Save_all")
            {
                std::cout << "\t\tSaving last states activated\n";
                run.save_interval = 1;
            }
            else if (Switch == "Finer_mesh")
            {
                std::cout << "\t\tSwitching to finer mesh\n";
                sim.Options.max_absolute_length = sim.Options.max_absolute_length / 2;
                sim.Options.min_absolute_length = sim.Options.min_absolute_length / 2;
                run.Switch_times_map[Switch] = -1;
            }
            else if (Switch == "Adapt_remesh")
            {
                run.adapt_remesh = true;
                run.remesh_every = 1;
            }
            else if (Switch == "Break_bonds")
            {
                for (Bead *bead : sim.Sim_handler.Beads)
                    bead->Beads.resize(0);
                std::cout << "Bonds broke\n";
            }
        }
    }

    // Area_constraint with nu > 0: move the target area by dA each step
    void update_area_target(Simulation &sim)
    {
        if (fabs(sim.dA) <= 1e-15)
            return;
        int i = index_of(sim.Energies, "Area_constraint");
        if (i < 0)
            return;
        double dA_effective = (sim.dA / (fabs(sim.dA))) * std::min(fabs(sim.dA), fabs(sim.A_bar - sim.Energy_constants[i][1]));
        sim.Energy_constants[i][1] += dA_effective;
        sim.Sim_handler.Energy_constants[i][1] += dA_effective;
        if (fabs(sim.A_bar - sim.Energy_constants[i][1]) < fabs(sim.dA))
        {
            std::cout << "\t\t\tReached target area\n";
            sim.dA = 0.0;
            sim.Energy_constants[i][1] = sim.A_bar;
            sim.Sim_handler.Energy_constants[i][1] = sim.A_bar;
        }
    }

    bool has_small_angle(Simulation &sim, double threshold)
    {
        bool found = false;
        sim.geometry->requireCornerAngles();
        for (Corner c : sim.mesh->corners())
        {
            if (sim.geometry->cornerAngles[c] < threshold)
            {
                found = true;
                break;
            }
        }
        sim.geometry->unrequireCornerAngles();
        return found;
    }

    void fix_small_angles(Simulation &sim)
    {
        remeshSmallAngles(*sim.mesh, *sim.geometry, sim.Options);
        sim.M3DG.BFGS_iter = 0;
        sim.mesh->compress();
        sim.geometry->refreshQuantities();
    }

    void maybe_remesh(Simulation &sim, RunState &run, size_t current_t)
    {
        if (!run.remesher)
            return;

        // Collapse slivers every step
        bool flagSmallAngle = has_small_angle(sim, 0.2);
        if (flagSmallAngle)
            fix_small_angles(sim);

        bool due = (current_t - run.last_remesh) > size_t(run.remesh_every) && run.remesh_every > 0;
        if (!(due || run.dt_sim == 0.0 || (flagSmallAngle && run.remesh_every < 0)))
            return;

        if (has_small_angle(sim, sim.Options.angleThresh))
            fix_small_angles(sim);

        run.last_remesh = current_t;
        run.remesh_op = remesh(*sim.mesh, *sim.geometry, sim.Options);
        sim.geometry->refreshQuantities();
        sim.M3DG.BFGS_iter = 0;
        sim.Sim_handler.update_vertex_normals();

        double output = 0.0;
        if (run.adapt_remesh)
        {
            // PID-like value, only logged; remesh_every follows the simple rule below
            int error = run.remesh_op - run.trgt_remesh_op;
            run.integral_error = run.integral_error + error * run.remesh_every;
            double derivative = (error - (run.remesh_op_last - run.trgt_remesh_op)) / run.remesh_every;
            output = 0.6 * error + 1.2 * run.integral_error / 10000 + (3.0 / 4.0) * derivative / 10000;
            run.remesh_op_last = run.remesh_op;

            if (run.remesh_op > 50)
                run.remesh_every = std::max(run.remesh_every / 2, 1);
            else if (run.remesh_op < 20)
                run.remesh_every = run.remesh_every + 10;
            else
                run.remesh_every = run.remesh_every + 1;
            run.remesh_every = std::min(run.remesh_every, 100);
        }
        if (sim.cfg.count_remesh)
        {
            std::ofstream Remeshing_count(run.basic_name + "Remeshing_count.txt", std::ios_base::app);
            Remeshing_count << current_t << " " << run.remesh_op << " " << run.remesh_every << " " << output << "\n";
        }
    }

    // Constraint energies are replaced by Lagrange multipliers in the Newton
    // methods: zero their constants in the handler (restore_constraint_energies undoes it).
    void disable_constraint_energies(Simulation &sim)
    {
        for (size_t i = 0; i < sim.Energies.size(); i++)
            if (sim.Energies[i] == "Area_constraint" || sim.Energies[i] == "Volume_constraint")
                sim.Sim_handler.Energy_constants[i][0] = 0.0;
    }

    void restore_constraint_energies(Simulation &sim)
    {
        for (size_t i = 0; i < sim.Energies.size(); i++)
        {
            if (sim.Energies[i] == "Volume_constraint" || sim.Energies[i] == "Area_constraint")
                sim.Sim_handler.Energy_constants[i][0] = sim.Energy_constants[i][0];
            if (sim.Energies[i] == "Edge_reg")
                sim.Sim_handler.Energy_constants[i][0] = 0.0;
        }
    }

    void push_rigid_constraints(std::vector<std::string> &Constraints)
    {
        for (const char *name : {"CMx", "CMy", "CMz", "Rx", "Ry", "Rz"})
            Constraints.push_back(name);
    }

    int switch_time(const RunState &run, const std::string &name)
    {
        auto it = run.Switch_times_map.find(name);
        return it == run.Switch_times_map.end() ? -1 : it->second;
    }

    void step_newton(Simulation &sim, RunState &run, size_t current_t, std::ofstream &Sim_data, bool save)
    {
        E_Handler &Sim_handler = sim.Sim_handler;
        int Switch_t = switch_time(run, "Newton");
        if (current_t == 0 || int(current_t) == Switch_t)
        {
            if (!sim.M3DG.boundary)
            {
                int ia = index_of(sim.Energies, "Area_constraint");
                Eigen::VectorXd Lagrange_mults = Eigen::VectorXd::Zero(ia >= 0 ? 2 : 1);
                Sim_handler.Lagrange_mult = Lagrange_mults;
                Sim_handler.Trgt_vol = sim.V_bar;
                Sim_handler.Trgt_area = ia >= 0 ? sim.Energy_constants[ia][1] : 0.0;
            }
            else
            {
                Sim_handler.Lagrange_mult = Eigen::VectorXd(0);
            }
            run.Switch_times_map["Newton"] = -1;
        }

        // Rebuilt every step (with a boundary it used to grow by six entries per step)
        std::vector<std::string> &Constraints = run.Constraints;
        Constraints.clear();
        if (!sim.M3DG.boundary)
            Constraints.push_back("Volume");
        for (const std::string &name : sim.Energies)
            if (name == "Area_constraint")
                Constraints.push_back("Area");
        disable_constraint_energies(sim);
        push_rigid_constraints(Constraints);
        Sim_handler.Constraints = Constraints;

        std::vector<std::string> Data_filenames(0);
        run.dt_sim = sim.M3DG.integrate_Newton(Sim_data, run.time, run.Bead_filenames, save, Constraints, Data_filenames);

        if (sim.M3DG.small_TS)
        {
            run.Integration = "Gradient_descent";
            run.Switch_times_map["Newton"] = current_t + 1;
            run.remesh_every = 1;
            run.remesher = true;
            run.Switch_times_map["No_remesh"] = current_t + 2;
            restore_constraint_energies(sim);
        }
    }

    void step_newton_normal(Simulation &sim, RunState &run, size_t current_t, std::ofstream &Sim_data, bool save)
    {
        E_Handler &Sim_handler = sim.Sim_handler;
        std::vector<std::string> &Constraints = run.Constraints;
        int Switch_t = switch_time(run, "Newton-Normal");
        if (current_t == 1 || int(current_t) == Switch_t)
        {
            Constraints.resize(0);
            Sim_handler.update_vertex_normals();
            if (!sim.M3DG.boundary)
            {
                Constraints.push_back("Volume");
                int ia = index_of(sim.Energies, "Area_constraint");
                if (ia >= 0)
                    Constraints.push_back("Area");
                Eigen::VectorXd Lagrange_mults = Eigen::VectorXd::Zero(ia >= 0 ? 8 : 7);
                Sim_handler.Trgt_area = ia >= 0 ? sim.Energy_constants[ia][1] : 0.0;
                push_rigid_constraints(Constraints);
                Sim_handler.Constraints = Constraints;

                // Initial volume multiplier from the least squares fit of the normal gradient
                Sim_handler.Calculate_Jacobian_Normal();
                Sim_handler.Calculate_gradient();
                size_t nV = Sim_handler.mesh->nVertices();
                Eigen::MatrixXd Jt = Sim_handler.Jacobian_constraints.transpose();
                Jt = Jt.block(0, 0, nV, 1).eval(); // eval: the block aliases Jt
                Eigen::VectorXd df(nV);
                for (size_t i = 0; i < nV; i++)
                    df(i) = dot(Sim_handler.Current_grad[i], Sim_handler.Vertex_normals[i]);
                Eigen::MatrixXd J = Jt.transpose();
                Eigen::VectorXd delta_lambdas = (J * Jt).colPivHouseholderQr().solve(J * df);
                Lagrange_mults(0) = delta_lambdas(0);
                Sim_handler.Lagrange_mult = Lagrange_mults;
                Sim_handler.Trgt_vol = sim.V_bar;
                std::cout << "Initial Lagrange multipliers " << Lagrange_mults.transpose() << "\n";
            }
            else
            {
                Sim_handler.Lagrange_mult = Eigen::VectorXd::Zero(6);
            }
            run.Switch_times_map["Newton-Normal"] = -1;
        }

        disable_constraint_energies(sim);
        if (Constraints.size() > size_t(Sim_handler.Lagrange_mult.size()))
            Sim_handler.Lagrange_mult = Eigen::VectorXd::Zero(Constraints.size());
        Sim_handler.Constraints = Constraints;

        std::vector<std::string> Data_filenames(0);
        run.dt_sim = sim.M3DG.integrate_Newton_Normal(Sim_data, run.time, run.Bead_filenames, save, Constraints, Data_filenames);

        if (sim.M3DG.small_TS)
        {
            std::cout << "We will do GD for the next 100 steps\n";
            run.Integration = "Gradient_descent";
            run.Switch_times_map["Newton-Normal"] = current_t + 100;
            run.remesh_every = 1;
            run.remesher = true;
            run.Switch_times_map["No_remesh"] = current_t + 50;
            restore_constraint_energies(sim);
        }
    }

    // The IpOpt solvers are disabled; these only keep the bookkeeping the switch did.
    void step_ipopt(Simulation &sim, RunState &run, size_t current_t, bool normal)
    {
        const std::string name = normal ? "IpOpt-Normal" : "IpOpt";
        int Switch_t = switch_time(run, name);
        if (current_t == 0 || int(current_t) == Switch_t)
            run.Switch_times_map[name] = -1;

        std::vector<std::string> Constraints(0);
        if (!sim.M3DG.boundary)
            Constraints = std::vector<std::string>{"Volume"};
        double A_target = 0.0;
        for (size_t i = 0; i < sim.Energies.size(); i++)
        {
            if (sim.Energies[i] == "Area_constraint")
            {
                Constraints.push_back("Area");
                A_target = sim.Energy_constants[i][1];
            }
        }
        disable_constraint_energies(sim);
        sim.Sim_handler.Constraints = Constraints;
        if (normal)
        {
            sim.Sim_handler.Trgt_area = A_target * 1.05;
            run.Integration = "Gradient_descent";
        }
        else
        {
            sim.Sim_handler.Trgt_area = sim.geometry->totalArea();
        }
        std::cout << "IpOpt is disabled in this build, no step taken\n";
    }

    void step_integrator(Simulation &sim, RunState &run, size_t current_t, std::ofstream &Sim_data, bool save)
    {
        Mem3DG &M3DG = sim.M3DG;
        const std::string &Integration = run.Integration;
        if (Integration == "Gradient_descent")
        {
            run.dt_sim = M3DG.integrate(Sim_data, run.time, run.Bead_filenames, save);
        }
        else if (Integration == "BFGS")
        {
            M3DG.m = sim.cfg.bfgs_saved_states;
            run.dt_sim = M3DG.integrate_BFGS(Sim_data, run.time, run.Bead_filenames, save);
            if (M3DG.remesh_flag)
            {
                run.remesh_every = 1;
                M3DG.remesh_flag = false;
                std::cout << "We need to remesh before we can continue with the iteration\n";
            }
        }
        else if (Integration == "Newton")
        {
            step_newton(sim, run, current_t, Sim_data, save);
        }
        else if (Integration == "Newton-Normal")
        {
            step_newton_normal(sim, run, current_t, Sim_data, save);
        }
        else if (Integration == "IpOpt" || Integration == "IpOpt-Normal")
        {
            step_ipopt(sim, run, current_t, Integration == "IpOpt-Normal");
        }
        else if (Integration == "BFGS-Normal")
        {
            if (M3DG.BFGS_iter == 0)
                sim.Sim_handler.update_vertex_normals();
            M3DG.m = sim.cfg.bfgs_saved_states;
            run.dt_sim = M3DG.integrate_BFGS_Normal(Sim_data, run.time, run.Bead_filenames, save);
            if (M3DG.remesh_flag)
            {
                run.remesh_every = 1;
                M3DG.remesh_flag = false;
                std::cout << "We need to remesh before we can continue with the iteration\n";
            }
            else
            {
                run.remesh_every = -1;
            }
        }
    }

    void print_status(Simulation &sim, const RunState &run, size_t current_t, double elapsed_ms)
    {
        double Volume = sim.geometry->totalVolume();
        double Area = sim.geometry->totalArea();
        double nu_obs = 3 * Volume / (4 * PI * pow(Area / (4 * PI), 1.5));
        std::cout << "t = " << current_t << " system time " << run.time << " vertices " << sim.mesh->nVertices()
                  << " mean edge length " << std::fixed << std::setprecision(10) << sim.geometry->meanEdgeLength()
                  << " reduced volume " << nu_obs << "\n";
        std::cout << "A thousand iterations took " << elapsed_ms << " miliseconds\n";
    }

    void save_final_state(Simulation &sim, const RunState &run)
    {
        if (sim.cfg.saving_states)
        {
            // Bead positions are not tracked any more; the file keeps its old layout
            std::ofstream beads_saved(run.basic_name + "Saved_bead_info.txt");
            beads_saved << std::setprecision(std::numeric_limits<double>::max_digits10);
            for (int k = 0; k < 6; k++)
            {
                for (size_t i = 0; i < sim.Beads.size(); i++)
                    beads_saved << 0 << " " << 0 << " " << 0 << " ";
                beads_saved << " \n";
            }
        }
        std::ofstream o(run.basic_name + "Final_state.obj");
        o << "#This is a meshfile from a saved state\n";
        for (Vertex v : sim.mesh->vertices())
        {
            Vector3 Pos = sim.geometry->inputVertexPositions[v];
            o << "v " << Pos.x << " " << Pos.y << " " << Pos.z << "\n";
        }
        for (Face f : sim.mesh->faces())
        {
            o << "f";
            for (Vertex v : f.adjacentVertices())
                o << " " << v.getIndex() + 1;
            o << "\n";
        }
    }

} // namespace

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <input.json> <Nsim> [--check]\n";
        return EXIT_FAILURE;
    }
    int Nsim = std::stoi(argv[2]);

    if (argc > 3 && std::string(argv[3]) == "--check")
    {
        // Only validate the input file
        try
        {
            SimConfig cfg = load_config(argv[1]);
            std::cout << "OK: " << cfg.energies.size() << " energies, " << cfg.beads.size() << " beads, "
                      << cfg.switches.size() << " switches\n";
            return EXIT_SUCCESS;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Invalid input file: " << e.what() << "\n";
            return EXIT_FAILURE;
        }
    }

    Simulation sim;
    try
    {
        build_simulation(load_config(argv[1]), sim);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error setting up the simulation: " << e.what() << "\n";
        return EXIT_FAILURE;
    }
    const SimConfig &cfg = sim.cfg;
    Mem3DG &M3DG = sim.M3DG;

    RunState run;
    run.Integration = cfg.integration;
    run.remesher = cfg.remeshing;
    run.remesh_every = cfg.remesh_every;
    run.adapt_remesh = cfg.adapt_remesh;
    run.save_interval = cfg.save_interval;
    for (size_t i = 0; i < cfg.switches.size(); i++)
    {
        run.Switch_times_map[cfg.switches[i]] = cfg.switch_times[i];
        std::cout << "The switch " << cfg.switches[i] << " happens at time " << cfg.switch_times[i] << "\n";
    }

    print_mesh_diagnostics(sim);

    try
    {
        run.basic_name = make_numbered_dir(cfg.first_dir);
        record_run(cfg.first_dir, run.basic_name, build_descriptive_name(sim, Nsim), cfg.source_path);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << "\n";
        return EXIT_FAILURE;
    }
    std::cout << "The output directory is " << run.basic_name << "\n";
    copy_file(cfg.source_path, run.basic_name + "Input_file.json");
    run.output_file = run.basic_name + "Output_data.txt";
    run.Bead_filenames = open_output_files(sim, run.basic_name, false);

    Save_mesh(sim.mesh, sim.geometry, run.basic_name, 0);
    M3DG.basic_name = run.basic_name;
    M3DG.BFGS_iter = 0;
    M3DG.Newton_iter = 0;
    apply_initial_perturbations(sim, run);

    if (!M3DG.boundary)
        sim.Sim_handler.Constraints.push_back("Volume_constraint");

    int coverage_index = index_of(sim.Energies, "Coverage");

    using clock = std::chrono::steady_clock;
    auto ms_since = [](clock::time_point t0)
    { return double(std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - t0).count()); };
    double remeshing_elapsed_time = 0, integrate_elapsed_time = 0, saving_mesh_time = 0;
    auto start = clock::now();

    std::ofstream Sim_data;
    const size_t Final_t = cfg.timesteps;
    for (size_t current_t = 0; current_t <= Final_t; current_t++)
    {
        M3DG.discreteTs = current_t;
        apply_switches(sim, run, current_t);
        update_area_target(sim);

        auto t0 = clock::now();
        maybe_remesh(sim, run, current_t);
        remeshing_elapsed_time += ms_since(t0);

        bool save = current_t % run.save_interval == 0;
        if (save)
        {
            t0 = clock::now();
            Save_mesh(sim.mesh, sim.geometry, run.basic_name, current_t);
            saving_mesh_time += ms_since(t0);
            Sim_data = std::ofstream(run.output_file, std::ios_base::app);
        }
        if (current_t % 1000 == 0)
        {
            print_status(sim, run, current_t, ms_since(start));
            start = clock::now();
        }

        t0 = clock::now();
        if (current_t % 100 == 0 && coverage_index >= 0)
            sim.Sim_handler.Debug_Coverage(sim.Sim_handler.Energy_constants[coverage_index]);

        step_integrator(sim, run, current_t, Sim_data, save);

        for (Bead *b : M3DG.Beads)
            if (b->state == "manual")
                b->update_state();

        if (M3DG.small_TS && current_t > Final_t * 0.2 && cfg.finish_sim)
        {
            std::cout << "Ending sim due to small TS at t = " << current_t << "\n";
            break;
        }
        if (run.dt_sim < 0)
        {
            std::cout << "Sim broke or timestep very small at timestep " << current_t << "\n";
            if (run.Integration != "BFGS")
                break;
            run.Integration = "BFGS-Normal";
            std::cout << "Switching to BFGS-Normal\n";
            M3DG.BFGS_iter = 0;
        }
        else
        {
            run.time += run.dt_sim;
            M3DG.system_time += 1;
        }
        Sim_data.close();
        integrate_elapsed_time += ms_since(t0);

        if (current_t % 1000 == 0)
        {
            std::cout << "Remeshing has taken " << remeshing_elapsed_time << " ms, saving meshes " << saving_mesh_time
                      << " ms, integrating " << integrate_elapsed_time << " ms\n";
            std::ifstream statm("/proc/self/statm");
            long pages;
            statm >> pages;
            std::cout << "RSS at t=" << current_t << ": " << pages * 4 / 1024 << " MB\n\n";
        }
    }
    std::cout << "The simulation is finished\n";
    Sim_data.close();

    save_final_state(sim, run);
    return EXIT_SUCCESS;
}
