// Batch driver: main_cluster <input.json> [--check]
//
// Loads the simulation with load_config/build_simulation (SimConfig.h), then
// runs the time loop: switches -> remeshing -> saving -> one integrator step.
//
// Resuming: with "continue_sim": true the run in first_dir/Subfolder goes on from
// its last recorded step N (the last row of Output_data.txt that has a mesh),
// appending to the same files, for "extra_steps" more steps (else "timesteps").
// The saved Input_file.json is the base config; the file passed in only gives
// extra_steps/timesteps and, if every saved switch is done, new Switches whose
// times count from N. See Config_files/Resume_example.json.

#include <sys/stat.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "geometrycentral/surface/remeshing.h"

#include <EigenRand/EigenRand>

#include "RemeshMonitor.h"
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
        double quality_baseline = 0.0; // bad edge fraction right after the last remesh
        size_t n_remesh = 0;
        std::string monitored; // integrator the convergence windows belong to
        int tiny_steps = 0;    // consecutive L-BFGS line searches that ended below 1e-10
        bool convergence_log = false;
        int polish_left = 0; // gradient descent steps left in the current polish
        int n_polish = 0;
        size_t remesh_pause_until = 0; // no remeshing before this step (after a rollback)
        size_t n_rollback = 0;

        double time = 0.0;
        double dt_sim = 0.0;
        bool resumed_first_step = false; // the Newton methods set up their multipliers on it
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
                std::cout << (Switch == "Freeze_beads" ? "Freezing" : "Freeing") << " the beads\n";
                // Freed beads with rigid bonds keep them
                for (Bead &bead : sim.Beads)
                    bead.state = Switch == "Freeze_beads" ? "froze" : (bead.Has_rigid_bond() ? "rigid" : "default");
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

    // Total energy for Remesh_log; the handler's per-term values are put back
    double total_energy(Simulation &sim)
    {
        double E = 0.0;
        energy_terms(sim, E);
        return E;
    }

    void log_polish(Simulation &sim, const RunState &run, size_t current_t, const std::string &event)
    {
        double E = total_energy(sim);
        std::cout << "Polish " << run.n_polish << ": " << event << " at t = " << current_t << ", E = " << E << "\n";
        const std::string path = run.basic_name + "Polish_log.txt";
        const bool fresh = !std::ifstream(path).good();
        std::ofstream log(path, std::ios_base::app);
        if (fresh)
            log << "# timestep polish event E (slivers: before the remesh, relaxed: after the gradient descent)\n";
        log << std::setprecision(12) << current_t << " " << run.n_polish << " " << event << " " << E << "\n";
    }

    // Slivers in BFGS-Normal: remesh (the caller does it) and relax with
    // gradient descent before going back (PolishParams). Returns true if started.
    bool start_polish(Simulation &sim, RunState &run, size_t current_t)
    {
        const PolishParams &p = sim.cfg.polish;
        if (run.Integration != "BFGS-Normal" || p.gd_steps == 0 || (p.max_cycles >= 0 && run.n_polish >= p.max_cycles))
            return false;
        run.n_polish++;
        log_polish(sim, run, current_t, "slivers");
        run.Integration = "Gradient_descent";
        run.polish_left = p.gd_steps;
        run.remesh_every = sim.cfg.remesh_every;
        return true;
    }

    // Called after every successful step
    void advance_polish(Simulation &sim, RunState &run, size_t current_t)
    {
        if (run.polish_left == 0)
            return;
        if (run.Integration != "Gradient_descent")
        {
            run.polish_left = 0; // a switch changed the integrator
            return;
        }
        if (--run.polish_left > 0)
            return;
        log_polish(sim, run, current_t, "relaxed");
        run.Integration = sim.cfg.polish.then;
        sim.M3DG.BFGS_iter = 0;
        if (run.Integration == "BFGS-Normal")
            run.remesh_every = -1;
    }

    // Undo of the remeshing of one step (RemeshRollbackParams): armed before
    // the mesh is first changed, checked once the step's remeshing is done
    struct RemeshUndo
    {
        bool armed = false;
        MeshSnapshot snap;
        double E_before = 0.0;
        std::vector<double> terms_before;
    };

    void arm_undo(Simulation &sim, RemeshUndo &undo)
    {
        if (undo.armed || sim.cfg.remesh_rollback.max_rise == 0.0)
            return;
        undo.snap = take_snapshot(sim);
        undo.terms_before = energy_terms(sim, undo.E_before);
        undo.armed = true;
    }

    // Put the mesh from before the remeshing back and pause remeshing. `cause`
    // is the energy term that exploded (its values before/after) or
    // "failed_step" (the total energy before/after)
    void undo_remesh(Simulation &sim, RunState &run, RemeshUndo &undo, size_t current_t, const std::string &cause,
                     double before, double after, double E)
    {
        const RemeshRollbackParams &p = sim.cfg.remesh_rollback;
        const size_t nV_rejected = sim.mesh->nVertices();
        restore_snapshot(sim, undo.snap);
        undo.armed = false;
        sim.M3DG.BFGS_iter = 0;
        sim.Sim_handler.update_vertex_normals();
        run.n_rollback++;
        run.remesh_pause_until = current_t + p.wait + 1;
        std::cout << "Remesh at t = " << current_t << " undone (" << cause << ": " << before << " -> " << after
                  << "), no remeshing for " << p.wait << " steps\n";

        const std::string path = run.basic_name + "Rollback_log.txt";
        const bool fresh = !std::ifstream(path).good();
        std::ofstream log(path, std::ios_base::app);
        if (fresh)
            log << "# timestep cause value_before value_rejected E_before E_rejected nVertices_rejected nVertices_kept "
                   "(cause: the energy term that exploded, or failed_step with the total energies)\n";
        log << std::setprecision(12) << current_t << " " << cause << " " << before << " " << after << " "
            << undo.E_before << " " << E << " " << nV_rejected << " " << sim.mesh->nVertices() << "\n";
    }

    // Returns true if the remeshing exploded an energy term and was undone
    bool check_undo(Simulation &sim, RunState &run, RemeshUndo &undo, size_t current_t)
    {
        if (!undo.armed)
            return false;
        double E = 0.0;
        const std::vector<double> terms = energy_terms(sim, E);
        const int i = exploded_term(undo.terms_before, terms, sim.cfg.remesh_rollback.max_rise);
        if (i < 0)
            return false;
        undo_remesh(sim, run, undo, current_t, sim.Energies[i], undo.terms_before[i], terms[i], E);
        return true;
    }

    // The integrator step after an accepted remesh failed: undo the remesh
    // (and the step) instead of ending the run
    void undo_failed_step(Simulation &sim, RunState &run, RemeshUndo &undo, size_t current_t)
    {
        double E = 0.0;
        energy_terms(sim, E);
        undo_remesh(sim, run, undo, current_t, "failed_step", undo.E_before, E, E);
        sim.M3DG.remesh_flag = false;
        sim.M3DG.step_failed = false;
    }

    // A remesh that is kept leaves `undo` armed until the step after it is done
    void maybe_remesh(Simulation &sim, RunState &run, size_t current_t, RemeshUndo &undo)
    {
        if (!run.remesher || current_t < run.remesh_pause_until)
            return;

        // Collapse slivers every step
        bool flagSmallAngle = has_small_angle(sim, 0.2);
        bool polish = flagSmallAngle && start_polish(sim, run, current_t);
        if (flagSmallAngle)
        {
            arm_undo(sim, undo);
            fix_small_angles(sim);
        }

        // In the quality mode remesh_every is the longest interval (max_every)
        // and setting it to 1 still asks for a remesh on the next step.
        const bool quality = run.adapt_remesh && sim.cfg.quality_remesh;
        const RemeshQualityParams &qp = sim.cfg.remesh_quality;
        const size_t since = current_t - run.last_remesh;

        bool due = since > size_t(run.remesh_every) && run.remesh_every > 0;
        bool go = due || run.dt_sim == 0.0 || (flagSmallAngle && run.remesh_every < 0) || polish;

        MeshQuality before;
        bool measured = false;
        if (quality && !go && run.remesh_every > 0 && since >= size_t(qp.min_every))
        {
            before = measure_mesh_quality(*sim.mesh, *sim.geometry, sim.Options);
            measured = true;
            go = before.bad_fraction() - run.quality_baseline > qp.f_tol;
        }
        if (!go)
        {
            check_undo(sim, run, undo, current_t);
            return;
        }

        double E_before = 0.0;
        if (sim.cfg.remesh_log)
        {
            if (!measured)
                before = measure_mesh_quality(*sim.mesh, *sim.geometry, sim.Options);
            E_before = total_energy(sim);
        }

        arm_undo(sim, undo);
        if (has_small_angle(sim, sim.Options.angleThresh))
            fix_small_angles(sim);

        run.last_remesh = current_t;
        run.n_remesh++;
        run.remesh_op = remesh(*sim.mesh, *sim.geometry, sim.Options);
        sim.geometry->refreshQuantities();
        sim.M3DG.BFGS_iter = 0;
        sim.Sim_handler.update_vertex_normals();
        if (check_undo(sim, run, undo, current_t))
            return;

        MeshQuality after;
        if (quality || sim.cfg.remesh_log)
            after = measure_mesh_quality(*sim.mesh, *sim.geometry, sim.Options);

        double output = 0.0;
        if (quality)
        {
            // Defects the remesher leaves behind (vetoed collapses, flips from
            // the final smoothing) are not counted against the next interval
            run.quality_baseline = after.bad_fraction();
            run.remesh_every = qp.max_every;
        }
        else if (run.adapt_remesh)
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
        if (sim.cfg.remesh_log)
        {
            double E_after = total_energy(sim);
            std::ofstream log(run.basic_name + "Remesh_log.txt", std::ios_base::app);
            log << std::setprecision(12) << current_t << " " << since << " " << run.remesh_op << " "
                << sim.mesh->nVertices() << " " << before.bad_fraction() << " " << before.fraction(before.n_long) << " "
                << before.fraction(before.n_short) << " " << before.fraction(before.n_flip) << " "
                << after.bad_fraction() << " " << after.fraction(after.n_long) << " " << after.fraction(after.n_short)
                << " " << after.fraction(after.n_flip) << " " << E_before << " " << E_after << "\n";
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
        if (current_t == 0 || int(current_t) == Switch_t || run.resumed_first_step)
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
        if (current_t == 1 || int(current_t) == Switch_t || run.resumed_first_step)
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

    // Stopping criteria of the L-BFGS integrators (ConvergenceMonitor.h), called
    // after every successful step. Each closed window goes to Convergence_log.txt.
    // Returns true when the run is converged and should end.
    bool monitor_convergence(Simulation &sim, RunState &run, ConvergenceMonitor &monitor, size_t current_t)
    {
        const bool normal = run.Integration == "BFGS-Normal";
        if (run.Integration != "BFGS" && !normal)
        {
            run.monitored.clear();
            return false;
        }
        if (run.monitored != run.Integration)
        {
            monitor.reset();
            run.monitored = run.Integration;
        }

        Mem3DG &M3DG = sim.M3DG;
        if (!monitor.add_step(M3DG.E_step_start, M3DG.E_step_end))
            return false;

        const StoppingParams &p = monitor.params();
        const double E = M3DG.E_step_end;
        const double dE_rel = monitor.window_dE_rel();
        // BFGS moves along the full gradient and may not have normals yet
        const double f_rms = force_density_rms(*sim.mesh, *sim.geometry, sim.Sim_handler.Current_grad,
                                               normal ? &sim.Sim_handler.Vertex_normals : nullptr);
        const double KB = bending_modulus(sim.Sim_handler.Energies, sim.Sim_handler.Energy_constants);
        const double f_scale = force_density_scale(*sim.mesh, *sim.geometry, KB);
        const bool force_ok = p.normal_tol_g == 0.0 || f_rms < p.normal_tol_g * f_scale;
        const bool pass = normal ? dE_rel < p.normal_tol_E && force_ok : dE_rel < p.bfgs_tol_E;
        const bool converged = monitor.close_window(pass);

        std::ofstream log(run.basic_name + "Convergence_log.txt", std::ios_base::app);
        if (!run.convergence_log)
        {
            log << "# timestep integration E dE_rel f_rms f_scale passed_windows "
                   "(f_rms: normal force density for BFGS-Normal, full for BFGS; f_scale = KB/R^3)\n";
            if (KB == 0.0)
                std::cout << "No bending energy is on, the force density is measured in units of 1/R^3\n";
            run.convergence_log = true;
        }
        log << std::setprecision(10) << current_t << " " << run.Integration << " " << E << " " << dE_rel << " "
            << f_rms << " " << f_scale << " " << monitor.passed_windows() << "\n";

        if (converged && !normal && p.bfgs_switch)
        {
            // Same as the "BFGS-Normal" switch
            std::cout << "The BFGS energy has plateaued at t = " << current_t << ", switching to BFGS-Normal\n";
            run.Integration = "BFGS-Normal";
            M3DG.BFGS_iter = 0;
            run.remesh_every = -1;
        }
        if (converged && normal && p.normal_stop)
        {
            std::cout << "BFGS-Normal converged at t = " << current_t << ": relative energy change " << dE_rel
                      << " per window, normal force density " << f_rms / f_scale << " KB/R^3\n";
            return true;
        }
        return false;
    }

    // L-BFGS line searches that keep ending at the smallest step: first reset
    // the history, then give up on the integrator. Returns true to end the run.
    bool handle_stall(Simulation &sim, RunState &run, size_t current_t)
    {
        const bool normal = run.Integration == "BFGS-Normal";
        const int n = sim.cfg.stopping.stall_steps;
        if (n == 0 || (run.Integration != "BFGS" && !normal) || run.dt_sim >= 1e-10)
        {
            run.tiny_steps = 0;
            return false;
        }
        run.tiny_steps++;
        if (run.tiny_steps == n)
        {
            std::cout << run.Integration << " stalled at t = " << current_t << ", resetting the L-BFGS history\n";
            sim.M3DG.BFGS_iter = 0;
        }
        else if (run.tiny_steps >= 2 * n)
        {
            run.tiny_steps = 0;
            if (normal)
            {
                std::cout << "BFGS-Normal is still stalled at t = " << current_t << ", ending the run\n";
                return true;
            }
            std::cout << "BFGS is still stalled at t = " << current_t << ", switching to BFGS-Normal\n";
            run.Integration = "BFGS-Normal";
            sim.M3DG.BFGS_iter = 0;
            run.remesh_every = -1;
        }
        return false;
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

    // Last configuration of the run, labelled with the step the loop would have
    // run next. The mesh, a row of Output_data.txt and the bead rows are all
    // written for this same state, whatever the save interval, so a finished run
    // can be loaded by main_visualize and its final energies read.
    void save_last_step(Simulation &sim, const RunState &run, size_t step)
    {
        Mem3DG &M3DG = sim.M3DG;
        M3DG.discreteTs = step;
        sim.geometry->refreshQuantities();
        sim.Sim_handler.Calculate_gradient(); // bead forces and gradient norms of this state
        double tot_E = 0;
        sim.Sim_handler.Calculate_energies(&tot_E);

        Save_mesh(sim.mesh, sim.geometry, run.basic_name, step);
        std::ofstream Sim_data(run.output_file, std::ios_base::app);
        M3DG.write_output_row(Sim_data, run.time, sim.geometry->totalVolume(), sim.geometry->totalArea(), tot_E, run.dt_sim);
        M3DG.write_bead_rows(run.Bead_filenames);
        std::cout << "Final state saved as step " << step << ", total energy " << tot_E << "\n";
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

    // ---------------------------------------------------------------------
    // Resuming a run ("continue_sim": true plus "Subfolder")
    // ---------------------------------------------------------------------

    using json = nlohmann::json;

    // What a resume needs, read from the run folder before anything is touched
    struct ResumeInfo
    {
        std::string dir;
        size_t step = 0; // N: the last recorded step, the loop restarts here
        double time = 0.0;
        size_t steps = 0;       // steps to run from N
        size_t final_step = 0;  // N + steps
        std::string mesh_file;
        std::vector<Vector3> bead_pos;

        // Saved switches that are not dismissed, in the saved order
        std::vector<std::string> old_switches;
        std::vector<int> old_times;
        std::vector<std::string> dismissed;
        std::vector<std::string> pending;
        bool use_new_switches = false;
        std::vector<std::string> new_switches;
        std::vector<int> new_times; // absolute

        json effective; // goes to the new Input_file.json
    };

    bool file_exists(const std::string &path)
    {
        struct stat info;
        return stat(path.c_str(), &info) == 0;
    }

    // name.ext, or name_b.ext, name_c.ext... when it is taken
    std::string unused_name(const std::string &stem, const std::string &ext)
    {
        std::string candidate = stem + ext;
        for (char suffix = 'b'; file_exists(candidate) && suffix <= 'z'; suffix++)
            candidate = stem + "_" + suffix + ext;
        if (file_exists(candidate))
            throw std::runtime_error("Too many files named " + stem + "*" + ext);
        return candidate;
    }

    std::vector<std::string> split_words(const std::string &line)
    {
        std::istringstream in(line);
        std::vector<std::string> words;
        std::string w;
        while (in >> w)
            words.push_back(w);
        return words;
    }

    bool parse_number(const std::string &word, double &value)
    {
        char *end = nullptr;
        value = std::strtod(word.c_str(), &end);
        return end != word.c_str() && *end == '\0';
    }

    // Data rows of a text output file: the words of every line whose first
    // `key_col + 1` words are numbers (headers and # lines are skipped)
    std::vector<std::vector<std::string>> read_rows(const std::string &path, size_t key_col)
    {
        std::vector<std::vector<std::string>> rows;
        std::ifstream in(path);
        std::string line;
        while (std::getline(in, line))
        {
            std::vector<std::string> words = split_words(line);
            double v;
            if (words.size() > key_col && parse_number(words[0], v) && parse_number(words[key_col], v))
                rows.push_back(words);
        }
        return rows;
    }

    // Keep headers and the rows with step < first_dropped, the rest is redone
    void drop_rows_from(const std::string &path, size_t key_col, size_t first_dropped)
    {
        if (!file_exists(path))
            return;
        const std::string tmp = path + ".tmp";
        {
            std::ifstream in(path);
            std::ofstream out(tmp);
            std::string line;
            while (std::getline(in, line))
            {
                std::vector<std::string> words = split_words(line);
                double v = 0;
                bool is_row = words.size() > key_col && parse_number(words[0], v) && parse_number(words[key_col], v);
                if (!is_row || v < double(first_dropped))
                    out << line << "\n";
            }
        }
        std::rename(tmp.c_str(), path.c_str());
    }

    std::string join_names(const std::vector<std::string> &names)
    {
        std::string s;
        for (const std::string &n : names)
            s += (s.empty() ? "" : ", ") + n;
        return s.empty() ? "none" : s;
    }

    // Reads the run folder and the launch file; changes nothing on disk.
    ResumeInfo plan_resume(const SimConfig &cfg)
    {
        if (!cfg.loaded_from_subfolder)
            throw std::runtime_error("\"continue_sim\" is true but the file has no \"Subfolder\"");
        if (cfg.subfolder_missing)
            throw std::runtime_error("\"continue_sim\" is true but the run folder " + cfg.subfolder_dir + " does not exist");

        ResumeInfo r;
        r.dir = cfg.subfolder_dir;
        const json &launch = cfg.launch_raw;

        // N: the last step with a row in Output_data.txt and its mesh
        std::vector<std::vector<std::string>> rows = read_rows(r.dir + "Output_data.txt", 1);
        bool found = false;
        for (size_t i = rows.size(); i-- > 0 && !found;)
        {
            size_t step = size_t(std::stod(rows[i][1]));
            if (file_exists(r.dir + "membrane_" + std::to_string(step) + ".obj"))
            {
                r.step = step;
                r.time = std::stod(rows[i][0]);
                found = true;
            }
        }
        if (!found)
            throw std::runtime_error("No step of " + r.dir + "Output_data.txt has a membrane_<step>.obj");
        r.mesh_file = r.dir + "membrane_" + std::to_string(r.step) + ".obj";

        // Bead positions at N (the last row up to N when N has none)
        for (size_t b = 0; b < cfg.beads.size(); b++)
        {
            std::string file = r.dir + "Bead_" + std::to_string(b) + "_data.txt";
            std::vector<std::vector<std::string>> brows = read_rows(file, 3);
            const std::vector<std::string> *best = nullptr;
            for (const std::vector<std::string> &w : brows)
                if (size_t(std::stod(w[0])) <= r.step)
                    best = &w;
            if (!best)
                throw std::runtime_error("No row up to step " + std::to_string(r.step) + " in " + file);
            if (size_t(std::stod((*best)[0])) != r.step)
                std::cout << "Warning: " << file << " has no row for step " << r.step << ", using step " << (*best)[0] << "\n";
            r.bead_pos.push_back(Vector3({std::stod((*best)[1]), std::stod((*best)[2]), std::stod((*best)[3])}));
        }

        // How many steps to run
        if (launch.contains("extra_steps"))
        {
            r.steps = launch["extra_steps"].get<size_t>();
        }
        else if (launch.contains("timesteps"))
        {
            r.steps = launch["timesteps"].get<size_t>();
            std::cout << "\"extra_steps\" is not specified, using \"timesteps\" = " << r.steps
                      << " as the number of steps to run\n";
        }
        else
        {
            throw std::runtime_error("Resuming needs \"extra_steps\" or \"timesteps\" in the file");
        }
        r.final_step = r.step + r.steps;

        // Saved switches: fired (< N), pending (N .. timesteps) or dismissed
        // (after the planned end of the saved run, or never)
        for (size_t i = 0; i < cfg.switches.size(); i++)
        {
            const int t = cfg.switch_times[i];
            if (t < 0 || size_t(t) > cfg.timesteps)
            {
                r.dismissed.push_back(cfg.switches[i] + "@" + std::to_string(t));
                continue;
            }
            r.old_switches.push_back(cfg.switches[i]);
            r.old_times.push_back(t);
            if (size_t(t) >= r.step)
                r.pending.push_back(cfg.switches[i] + "@" + std::to_string(t));
        }

        std::vector<std::string> launch_names;
        std::vector<int> launch_times;
        if (launch.contains("Switches"))
        {
            launch_names = launch["Switches"].get<std::vector<std::string>>();
            if (!launch.contains("Switch_times"))
                throw std::runtime_error("The launch file has \"Switches\" but no \"Switch_times\"");
            launch_times = launch["Switch_times"].get<std::vector<int>>();
            if (launch_names.size() != launch_times.size())
                throw std::runtime_error("The launch file's \"Switches\" and \"Switch_times\" have different lengths");
            for (const std::string &name : launch_names)
                if (std::find(known_switches().begin(), known_switches().end(), name) == known_switches().end())
                    throw std::runtime_error("The launch file has an unknown switch " + name);
        }

        std::cout << "Resuming " << r.dir << " from step " << r.step << " (time " << r.time << "), running to step "
                  << r.final_step << "\n";
        if (!r.dismissed.empty())
            std::cout << "Dismissed saved switches (after the planned end " << cfg.timesteps << " or never): "
                      << join_names(r.dismissed) << "\n";
        if (!r.pending.empty())
        {
            std::cout << "The saved run did not reach all its switches, pending: " << join_names(r.pending) << "\n";
            if (!launch_names.empty())
                std::cout << "Warning: the Switches of the launch file are ignored (" << join_names(launch_names) << ")\n";
        }
        else
        {
            r.use_new_switches = true;
            std::cout << "All saved switches are done, using the launch file's: " << join_names(launch_names) << "\n";
            for (size_t i = 0; i < launch_names.size(); i++)
            {
                r.new_switches.push_back(launch_names[i]);
                r.new_times.push_back(int(r.step) + launch_times[i]);
            }
        }

        // The config that will be saved: switches accumulate so a later resume can replay them
        r.effective = cfg.raw;
        r.effective["timesteps"] = std::max(cfg.timesteps, r.final_step);
        r.effective["continue_sim"] = false;
        r.effective["resumed_from"] = r.step;
        std::vector<std::string> eff_names = r.old_switches;
        std::vector<int> eff_times = r.old_times;
        for (size_t i = 0; i < r.new_switches.size(); i++)
        {
            auto same = std::find(eff_names.begin(), eff_names.end(), r.new_switches[i]);
            if (same != eff_names.end())
            {
                std::cout << "Warning: " << r.new_switches[i] << " was already a switch of this run, its saved time is replaced\n";
                eff_times[same - eff_names.begin()] = r.new_times[i];
            }
            else
            {
                eff_names.push_back(r.new_switches[i]);
                eff_times.push_back(r.new_times[i]);
            }
        }
        if (eff_names.empty())
        {
            r.effective.erase("Switches");
            r.effective.erase("Switch_times");
        }
        else
        {
            r.effective["Switches"] = eff_names;
            r.effective["Switch_times"] = eff_times;
        }
        return r;
    }

    // Everything that changes the run folder, once the simulation has been built
    void commit_resume(const ResumeInfo &r, const std::string &launch_path)
    {
        // Step N is recomputed, so its rows (and later ones) go
        drop_rows_from(r.dir + "Output_data.txt", 1, r.step);
        for (size_t b = 0; b < r.bead_pos.size(); b++)
            drop_rows_from(r.dir + "Bead_" + std::to_string(b) + "_data.txt", 0, r.step);
        drop_rows_from(r.dir + "Remeshing_count.txt", 0, r.step);
        drop_rows_from(r.dir + "Remesh_log.txt", 0, r.step);

        // Input file history: Input_file.json is always the latest
        const std::string input = r.dir + "Input_file.json";
        if (file_exists(input))
            std::rename(input.c_str(), unused_name(r.dir + "Input_file_" + std::to_string(r.step), ".json").c_str());
        {
            std::ofstream out(input);
            out << r.effective.dump(4) << "\n";
        }
        copy_file(launch_path, unused_name(r.dir + "Resume_" + std::to_string(r.step), ".json"));

        std::ofstream log(r.dir + "Resume_log.txt", std::ios_base::app);
        log << "step " << r.step << " time " << r.time << " steps " << r.steps << " final " << r.final_step
            << " launch_file " << launch_path << " switches " << (r.use_new_switches ? "new" : "saved") << "\n";
    }

    // Silences std::cout while the saved switches are replayed
    struct CoutSilencer
    {
        std::ostringstream sink;
        std::streambuf *old;
        CoutSilencer() : old(std::cout.rdbuf(sink.rdbuf())) {}
        ~CoutSilencer() { std::cout.rdbuf(old); }
    };

    // Rebuilds what the saved switches did before step N (the integrator,
    // bead states, remesher settings, extra energies, the area target)
    void replay_state(Simulation &sim, RunState &run, size_t step)
    {
        CoutSilencer quiet;
        for (size_t t = 0; t < step; t++)
        {
            apply_switches(sim, run, t);
            update_area_target(sim);
        }
    }

} // namespace

int main(int argc, char **argv)
{
    std::string config_path;
    bool check_only;
    if (!parse_command_line(argc, argv, config_path, check_only))
        return EXIT_FAILURE;

    if (check_only)
    {
        // Only validate the input file
        try
        {
            SimConfig cfg = load_config(config_path);
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
    bool resuming = false;
    ResumeInfo resume;
    try
    {
        resuming = wants_continue(config_path);
        SimConfig loaded = load_config(config_path, resuming);
        if (resuming)
        {
            resume = plan_resume(loaded);
            for (size_t i = 0; i < loaded.beads.size(); i++)
                loaded.beads[i].pos = resume.bead_pos[i];
            // The saved schedule (without the dismissed switches) is replayed below
            loaded.switches = resume.old_switches;
            loaded.switch_times = resume.old_times;
        }
        build_simulation(loaded, sim, resuming ? resume.mesh_file : "");
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

    size_t start_t = 0;
    if (resuming)
    {
        try
        {
            commit_resume(resume, config_path);
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << "\n";
            return EXIT_FAILURE;
        }
        run.basic_name = resume.dir;
        std::cout << "The output directory is " << run.basic_name << " (appending)\n";
        run.output_file = run.basic_name + "Output_data.txt";
        run.Bead_filenames = open_output_files(sim, run.basic_name, true);
        start_t = resume.step;
    }
    else
    {
        try
        {
            run.basic_name = make_numbered_dir(cfg.first_dir);
            record_run(cfg.first_dir, run.basic_name, build_descriptive_name(sim), cfg.source_path);
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
    }
    M3DG.basic_name = run.basic_name;
    M3DG.BFGS_iter = 0;
    M3DG.Newton_iter = 0;
    if (resuming)
    {
        // Redo what the saved switches did before step N, then continue with the schedule
        replay_state(sim, run, start_t);
        if (resume.use_new_switches)
        {
            sim.cfg.switches = resume.new_switches;
            sim.cfg.switch_times = resume.new_times;
            run.Switch_times_map.clear();
            for (size_t i = 0; i < cfg.switches.size(); i++)
                run.Switch_times_map[cfg.switches[i]] = cfg.switch_times[i];
        }
        for (size_t i = 0; i < cfg.switches.size(); i++)
            std::cout << "The switch " << cfg.switches[i] << " happens at step " << cfg.switch_times[i] << "\n";
        std::cout << "The integrator is " << run.Integration << " at step " << start_t << "\n";
        run.time = resume.time;
        run.last_remesh = start_t;
        run.dt_sim = 1.0; // a zero dt_sim forces a remesh on the first step, the saved mesh is already remeshed
        run.resumed_first_step = true;
        M3DG.system_time = start_t;
        sim.Sim_handler.update_vertex_normals();
    }
    else
    {
        apply_initial_perturbations(sim, run);
    }

    if (!M3DG.boundary)
        sim.Sim_handler.Constraints.push_back("Volume_constraint");

    int coverage_index = index_of(sim.Energies, "Coverage");
    ConvergenceMonitor monitor(cfg.stopping);

    using clock = std::chrono::steady_clock;
    auto ms_since = [](clock::time_point t0)
    { return double(std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - t0).count()); };
    double remeshing_elapsed_time = 0, integrate_elapsed_time = 0, saving_mesh_time = 0;
    auto start = clock::now();
    const auto run_start = clock::now();

    std::ofstream Sim_data;
    const size_t Final_t = resuming ? resume.final_step : cfg.timesteps;
    size_t last_t = start_t;
    for (size_t current_t = start_t; current_t <= Final_t; current_t++)
    {
        last_t = current_t;
        M3DG.discreteTs = current_t;
        apply_switches(sim, run, current_t);
        update_area_target(sim);

        auto t0 = clock::now();
        RemeshUndo undo;
        maybe_remesh(sim, run, current_t, undo);
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
        run.resumed_first_step = false;

        for (Bead *b : M3DG.Beads)
            if (b->state == "manual")
                b->update_state();

        if (M3DG.small_TS && current_t - start_t > (Final_t - start_t) * 0.2 && cfg.finish_sim)
        {
            std::cout << "Ending sim due to small TS at t = " << current_t << "\n";
            break;
        }
        if (run.dt_sim < 0 && undo.armed)
        {
            undo_failed_step(sim, run, undo, current_t);
        }
        else if (run.dt_sim < 0 && M3DG.step_failed)
        {
            std::cout << "The line search hit a nan at timestep " << current_t << ", ending the run\n";
            break;
        }
        else if (run.dt_sim < 0)
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
            if (handle_stall(sim, run, current_t) || monitor_convergence(sim, run, monitor, current_t))
            {
                Sim_data.close();
                break;
            }
            advance_polish(sim, run, current_t);
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
    if (start_t <= Final_t)
        save_last_step(sim, run, last_t + 1);

    {
        // Cost side of the remeshing calibration (Scripts/calibrate_remesh_tol.py)
        std::ofstream timing(run.basic_name + (resuming ? "Timing_from_" + std::to_string(start_t) + ".txt" : "Timing.txt"));
        timing << "wall_ms " << ms_since(run_start) << "\nremesh_ms " << remeshing_elapsed_time << "\nintegrate_ms "
               << integrate_elapsed_time << "\nsave_ms " << saving_mesh_time << "\nremeshes " << run.n_remesh << "\nremesh_rollbacks " << run.n_rollback
               << "\nvertices " << sim.mesh->nVertices() << "\nlbfgs_pairs " << M3DG.lbfgs_pairs
               << "\nlbfgs_negative_sy " << M3DG.lbfgs_negative_sy << "\nlbfgs_restarts " << M3DG.lbfgs_restarts
               << "\nlbfgs_uphill " << M3DG.lbfgs_uphill << "\n";
    }
    if (M3DG.lbfgs_pairs > 0)
        std::cout << "L-BFGS: " << M3DG.lbfgs_pairs << " curvature pairs, " << M3DG.lbfgs_negative_sy
                  << " with s.y < 0, " << M3DG.lbfgs_restarts << " restarts, " << M3DG.lbfgs_uphill
                  << " uphill directions replaced\n";

    save_final_state(sim, run);
    return EXIT_SUCCESS;
}
