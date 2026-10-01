#include "SimConfig.h"

#include <sys/stat.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "geometrycentral/surface/meshio.h"

using json = nlohmann::json;

namespace
{

    template <class T, class... Args>
    std::unique_ptr<T> make_unique_ptr(Args &&...args)
    {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }

    const json &require(const json &j, const std::string &key, const std::string &where)
    {
        if (!j.contains(key) || j.at(key).is_null())
            throw std::runtime_error("Input file: missing required key \"" + key + "\" in " + where);
        return j.at(key);
    }

    Vector3 to_vector3(const json &j, const std::string &key, const std::string &where)
    {
        const json &v = require(j, key, where);
        if (!v.is_array() || v.size() < 3)
            throw std::runtime_error("Input file: \"" + key + "\" in " + where + " must be an array of 3 numbers");
        return Vector3({v[0].get<double>(), v[1].get<double>(), v[2].get<double>()});
    }

    BeadSpec parse_bead(const json &b, size_t index)
    {
        const std::string where = "Beads[" + std::to_string(index) + "]";
        BeadSpec spec;
        spec.gradient_order = b.value("gradient_order", std::string("Bead"));
        if (b.contains("Constraint"))
        {
            spec.constraint = b["Constraint"].get<std::string>();
            spec.constraint_constants = require(b, "Constraint_constants", where).get<std::vector<double>>();
        }

        spec.pos = to_vector3(b, "Pos", where);
        spec.radius = require(b, "radius", where).get<double>();
        spec.state = require(b, "state", where).get<std::string>();
        spec.inter_str = require(b, "inter_str", where).get<double>();
        spec.mem_inter = require(b, "mem_inter", where).get<std::string>();
        const auto &known = known_interactions();
        if (std::find(known.begin(), known.end(), spec.mem_inter) == known.end())
        {
            std::string list;
            for (const std::string &name : known)
                list += (list.empty() ? "" : ", ") + name;
            throw std::runtime_error("Input file: unknown mem_inter \"" + spec.mem_inter + "\" in " + where +
                                     ". Known: " + list);
        }

        if (b.contains("rc"))
        {
            spec.has_rc = true;
            spec.rc = b["rc"].get<double>();
        }
        spec.outside = b.value("outside", 1.0);
        if (b.contains("shift"))
        {
            spec.has_shift = true;
            spec.shift = b["shift"].get<double>();
        }
        if (spec.mem_inter == "Gravity" || spec.mem_inter == "Pinch")
        {
            spec.z_axis = require(b, "Z_Axis", where).get<std::vector<double>>();
            if (spec.z_axis.size() < 3)
                throw std::runtime_error("Input file: \"Z_Axis\" in " + where + " must have 3 entries");
        }

        spec.bonds = require(b, "bonds", where).get<std::vector<std::string>>();
        spec.bonds_constants = require(b, "bonds_constants", where).get<std::vector<std::vector<double>>>();
        spec.partners = require(b, "Beads", where).get<std::vector<int>>();
        if (spec.partners.size() != spec.bonds.size() || spec.bonds.size() != spec.bonds_constants.size())
            throw std::runtime_error("Input file: " + where + " needs one \"bonds\" type and one \"bonds_constants\" "
                                     "entry per bead listed in \"Beads\"");

        if (spec.state == "manual")
        {
            spec.has_velocity = true;
            spec.velocity = to_vector3(b, "Velocity", where);
            if (b.contains("FinalPos"))
            {
                spec.has_final_pos = true;
                spec.final_pos = to_vector3(b, "FinalPos", where);
            }
        }
        return spec;
    }

    json read_json(const std::string &path)
    {
        std::ifstream file(path);
        if (!file.is_open())
            throw std::runtime_error("Could not open input file " + path);
        return json::parse(file);
    }

} // namespace

const std::vector<std::string> &known_interactions()
{
    static const std::vector<std::string> names = {
        "Gravity", "Pinch", "Frenkel", "Frenkel_Normal_nopush", "Linear", "Linear_Normal",
        "LJ", "Shifted-LJ", "One_over_r_x", "One_over_r", "None"};
    return names;
}

const std::vector<std::string> &known_switches()
{
    static const std::vector<std::string> names = {
        "Newton", "Newton-Normal", "BFGS", "BFGS-Normal", "IpOpt", "IpOpt-Normal",
        "Freeze_beads", "Free_beads", "No_remesh", "Remesh_always", "Restore_remeshing",
        "Volume_constraint", "Save_all", "Finer_mesh", "Adapt_remesh", "Break_bonds"};
    return names;
}

SimConfig load_config(const std::string &path, bool resolve_subfolder)
{
    SimConfig cfg;
    cfg.source_path = path;
    json data = read_json(path);

    if (resolve_subfolder && data.contains("Subfolder"))
    {
        std::string first_dir = require(data, "first_dir", "the top level").get<std::string>();
        cfg.subfolder_dir = first_dir + data["Subfolder"].get<std::string>();
        cfg.loaded_from_subfolder = true;
        cfg.source_path = cfg.subfolder_dir + "/Input_file.json";
        std::cout << "Loading the run in " << cfg.subfolder_dir << "\n";
        data = read_json(cfg.source_path);
    }
    const std::string top = "the top level";

    cfg.init_file = require(data, "init_file", top).get<std::string>();
    cfg.first_dir = require(data, "first_dir", top).get<std::string>();
    cfg.timesteps = require(data, "timesteps", top).get<size_t>();
    cfg.save_interval = require(data, "save_interval", top).get<int>();
    cfg.recentering = require(data, "recentering", top).get<bool>();
    cfg.boundary = require(data, "boundary", top).get<bool>();

    if (data.contains("finish_sim"))
        cfg.finish_sim = data["finish_sim"].get<bool>();
    else
        std::cout << "The simulation will finish when the timestep decreases\n";
    cfg.saving_states = data.value("saving_states", false);

    if (data.contains("remeshing"))
        cfg.remeshing = data["remeshing"].get<bool>();
    else
        std::cout << "Warning: \"remeshing\" not given, remeshing is on\n";
    cfg.remesh_every = data.value("remesh_every", 1);
    cfg.adapt_remesh = data.value("adapt_remesh", true);
    cfg.count_remesh = data.value("Count_remesh", false);
    if (data.contains("remesher"))
    {
        const json &r = data["remesher"];
        cfg.remesher.present = true;
        cfg.remesher.size_max = require(r, "size_max", "remesher").get<double>();
        cfg.remesher.size_min = require(r, "size_min", "remesher").get<double>();
        cfg.remesher.refine_angle = require(r, "refine_angle", "remesher").get<double>();
        cfg.remesher.aspect_min = require(r, "aspect_min", "remesher").get<double>();
    }

    if (data.contains("Integration"))
        cfg.integration = data["Integration"].get<std::string>();
    else
        std::cout << "The integration method is not defined, using Gradient descent\n";
    cfg.bfgs_saved_states = data.value("BFGS_saved_states", 10);

    if (data.contains("Switches"))
    {
        cfg.switches = data["Switches"].get<std::vector<std::string>>();
        cfg.switch_times = require(data, "Switch_times", top).get<std::vector<int>>();
        if (cfg.switch_times.size() != cfg.switches.size())
            throw std::runtime_error("Input file: \"Switches\" and \"Switch_times\" have different lengths");
        for (const std::string &sw : cfg.switches)
        {
            const auto &known = known_switches();
            if (std::find(known.begin(), known.end(), sw) == known.end())
                std::cout << "Warning: unknown switch \"" << sw << "\" will be ignored\n";
        }
    }

    if (data.contains("Displacement"))
        cfg.displacement = to_vector3(data, "Displacement", top);
    cfg.rescale = data.value("rescale", 1.0);
    if (data.contains("Initial_noise"))
    {
        cfg.has_initial_noise = true;
        cfg.initial_noise = data["Initial_noise"].get<double>();
    }
    cfg.harmonic = data.contains("Harmonic");

    if (data.contains("Field"))
    {
        cfg.has_field = true;
        cfg.field = data["Field"].get<std::string>();
        cfg.field_vals = require(data, "Field_vals", top).get<std::vector<double>>();
    }
    // Timestep is only read together with backtrack (as before)
    if (data.contains("backtrack"))
    {
        cfg.has_backtrack = true;
        cfg.backtrack = data["backtrack"].get<bool>();
        cfg.timestep = data.value("Timestep", 1e-4);
    }
    if (data.contains("momentum"))
    {
        cfg.has_momentum = true;
        cfg.momentum = data["momentum"].get<bool>();
    }

    if (data.contains("Energies"))
    {
        size_t i = 0;
        for (const json &e : data["Energies"])
        {
            const std::string where = "Energies[" + std::to_string(i++) + "]";
            cfg.energies.push_back({require(e, "Name", where).get<std::string>(),
                                    require(e, "constants", where).get<std::vector<double>>()});
        }
    }
    if (data.contains("Beads"))
    {
        size_t i = 0;
        for (const json &b : data["Beads"])
            cfg.beads.push_back(parse_bead(b, i++));
        for (size_t k = 0; k < cfg.beads.size(); k++)
            for (int partner : cfg.beads[k].partners)
                if (partner < 0 || size_t(partner) >= cfg.beads.size() || size_t(partner) == k)
                    throw std::runtime_error("Input file: Beads[" + std::to_string(k) + "] is bonded to bead " +
                                             std::to_string(partner) + ", which does not exist");
    }

    cfg.raw = data;
    return cfg;
}

Simulation::~Simulation()
{
    // Everything holding MeshData must go before the mesh itself
    M3DG = Mem3DG();
    Sim_handler = E_Handler();
    Options = RemeshOptions();
    Beads.clear();
    Interaction_container.clear();
    delete geometry;
    delete mesh;
}

std::unique_ptr<Interaction> make_interaction(const BeadSpec &spec, ManifoldSurfaceMesh *mesh,
                                              VertexPositionGeometry *geometry,
                                              std::vector<double> &params)
{
    const std::string &type = spec.mem_inter;
    if (type == "Gravity" || type == "Pinch")
    {
        params.push_back(spec.pos.z);
        params.push_back(-1);
        params.push_back(spec.z_axis[0]);
        params.push_back(spec.z_axis[1]);
        params.push_back(spec.z_axis[2]);
        if (type == "Gravity")
            return make_unique_ptr<Gravity_Plane>(mesh, geometry, params);
        return make_unique_ptr<Pinch_Interaction>(mesh, geometry, params);
    }
    if (type == "Frenkel" || type == "Frenkel_Normal_nopush" || type == "Linear" ||
        type == "Linear_Normal" || type == "One_over_r_x")
    {
        params.push_back(spec.outside);
        if (type == "Frenkel")
            return make_unique_ptr<Frenkel>(mesh, geometry, params);
        if (type == "Frenkel_Normal_nopush")
            return make_unique_ptr<Frenkel_Normal>(mesh, geometry, params);
        if (type == "Linear")
            return make_unique_ptr<Linear>(mesh, geometry, params);
        if (type == "Linear_Normal")
            return make_unique_ptr<Linear_Normal>(mesh, geometry, params);
        return make_unique_ptr<One_over_r_Normal>(mesh, geometry, params);
    }
    if (type == "LJ")
    {
        params.push_back(spec.shift);
        return make_unique_ptr<LJ>(mesh, geometry, params);
    }
    if (type == "Shifted-LJ")
    {
        // LJ cut at rc and shifted so the energy is zero at the cutoff (WCA
        // when rc = 2^(1/6) sigma, the default). An explicit "shift" wins.
        double epsilon = params[0], sigma = params[1], rc = params[2];
        if (rc <= 0)
            throw std::runtime_error("Input file: Shifted-LJ needs a positive cutoff rc");
        double sr6 = pow(sigma / rc, 6);
        params.push_back(spec.has_shift ? spec.shift : -4 * epsilon * (sr6 * sr6 - sr6));
        return make_unique_ptr<LJ>(mesh, geometry, params);
    }
    if (type == "One_over_r")
        return make_unique_ptr<One_over_r>(mesh, geometry, params);
    if (type == "None")
        return make_unique_ptr<No_mem_Inter>();

    throw std::runtime_error("Input file: unknown bead interaction mem_inter = \"" + type + "\"");
}

namespace
{

    // Default cutoff when "rc" is not given
    double default_rc(const BeadSpec &spec)
    {
        if (spec.mem_inter == "LJ" || spec.mem_inter == "Shifted-LJ")
            return spec.radius * pow(2, 1.0 / 6.0);
        if (spec.mem_inter == "Frenkel" || spec.mem_inter == "Frenkel_Normal_nopush")
            return spec.radius * 2.0;
        return -1;
    }

    // Post-process the constants that depend on the initial mesh
    void setup_energies(Simulation &sim)
    {
        VertexPositionGeometry *geometry = sim.geometry;
        sim.V_bar = geometry->totalVolume();
        sim.A_bar = geometry->totalArea();

        for (const EnergySpec &energy : sim.cfg.energies)
        {
            std::vector<double> constants = energy.constants;
            std::cout << "The constants for " << energy.name << " are ";
            for (double c : constants)
                std::cout << c << " ";
            std::cout << " \n";

            if (energy.name == "Volume_constraint")
            {
                // KV V_bar, a negative V_bar means the current volume
                if (constants[1] < 0)
                    constants[1] = geometry->totalVolume();
                sim.V_bar = constants[1];
            }
            if (energy.name == "Area_constraint")
            {
                // KA A_bar nu dA. With nu > 0 the target comes from the reduced volume
                // and is approached in steps of dA
                if (constants[2] > 0)
                {
                    double nu = constants[2];
                    sim.A_bar = pow(36 * PI * sim.V_bar * sim.V_bar / (nu * nu), 1.0 / 3.0);
                    double area = geometry->totalArea();
                    sim.dA = constants[3];
                    if ((sim.A_bar - area) * sim.dA < 0)
                        sim.dA = -sim.dA;
                    constants[1] = area + (sim.dA / (fabs(sim.dA))) * std::min(fabs(sim.dA), fabs(sim.A_bar - area));
                    std::cout << "The target Area is " << sim.A_bar << " from nu " << nu << " (current " << area << ")\n";
                }
                else if (constants[1] > 0)
                {
                    sim.A_bar = constants[1];
                    sim.dA = 0;
                }
                else
                {
                    sim.dA = 0;
                    sim.A_bar = geometry->totalArea();
                    constants[1] = sim.A_bar;
                }
                std::cout << "The target area is " << sim.A_bar << "\n";
            }
            if (energy.name == "Membrane_tension" || energy.name == "Excess_tension")
            {
                double A0 = geometry->totalArea();
                constants[1] = A0 * constants[1];
            }
            sim.Energies.push_back(energy.name);
            sim.Energy_constants.push_back(constants);
        }
    }

    void setup_beads(Simulation &sim)
    {
        const std::vector<BeadSpec> &specs = sim.cfg.beads;
        // Reserve so the pointers handed out below stay valid
        sim.Beads.reserve(specs.size());
        sim.Interaction_container.reserve(specs.size());

        for (size_t i = 0; i < specs.size(); i++)
        {
            const BeadSpec &spec = specs[i];
            sim.Energies.push_back(spec.gradient_order);
            sim.Energy_constants.push_back({});

            std::vector<double> params = {spec.inter_str, spec.radius, spec.has_rc ? spec.rc : default_rc(spec)};
            sim.Interaction_container.push_back(make_interaction(spec, sim.mesh, sim.geometry, params));

            sim.Beads.push_back(Bead());
            Bead &bead = sim.Beads.back();
            bead.mesh = sim.mesh;
            bead.geometry = sim.geometry;
            bead.Pos = spec.pos;
            bead.strength = params[0];
            bead.sigma = params[1];
            bead.rc = params[2];
            bead.interaction = spec.mem_inter;
            bead.Bead_I = sim.Interaction_container.back().get();
            bead.Bead_I->Bead_1 = &bead;
            bead.Bead_id = i;

            bead.Bond_type = spec.bonds;
            bead.Interaction_constants_vector = spec.bonds_constants;
            bead.state = spec.state;
            bead.Constraint = spec.constraint;
            bead.Constraint_constants = spec.constraint_constants;
            if (spec.has_velocity)
                bead.Velocity = spec.velocity;
            if (spec.has_final_pos)
                bead.FinalPos = spec.final_pos;
            bead.CoverageForce = Vector3({0.0, 0.0, 0.0});

            std::cout << "Bead " << i << ": " << spec.mem_inter << " radius " << bead.sigma << " cutoff " << bead.rc
                      << " state " << bead.state << "\n";
        }

        // Bonds: bead i is connected to every bead listed in its "Beads" entry
        for (size_t i = 0; i < specs.size(); i++)
        {
            for (int partner : specs[i].partners)
                sim.Beads[i].Beads.push_back(&sim.Beads[partner]);
        }

        for (size_t i = 0; i < sim.Beads.size(); i++)
        {
            sim.Beads[i].Total_beads = sim.Beads.size();
            sim.Beads[i].Bead_id = i;
        }
    }

} // namespace

void build_simulation(const SimConfig &cfg, Simulation &sim)
{
    sim.cfg = cfg;

    std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
    std::unique_ptr<VertexPositionGeometry> geometry_uptr;
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(cfg.init_file);
    sim.mesh = mesh_uptr.release();
    sim.geometry = geometry_uptr.release();

    if (cfg.displacement.norm() > 0)
    {
        std::cout << "Displacing the membrane by " << cfg.displacement << " \n";
        sim.geometry->normalize(cfg.displacement);
    }
    sim.geometry->rescale(cfg.rescale);
    sim.geometry->refreshQuantities();

    setup_energies(sim);
    setup_beads(sim);

    Mem3DG &M3DG = sim.M3DG;
    E_Handler &Sim_handler = sim.Sim_handler;
    M3DG = Mem3DG(sim.mesh, sim.geometry);
    Sim_handler = E_Handler(sim.mesh, sim.geometry, sim.Energies, sim.Energy_constants);
    Sim_handler.Trgt_vol = sim.V_bar;
    Sim_handler.Trgt_area = sim.A_bar;
    M3DG.recentering = cfg.recentering;
    M3DG.boundary = cfg.boundary;
    Sim_handler.boundary = cfg.boundary;
    for (Bead &bead : sim.Beads)
    {
        M3DG.Add_bead(&bead);
        Sim_handler.Add_Bead(&bead);
    }
    M3DG.Sim_handler = &Sim_handler;

    M3DG.Field = cfg.has_field ? cfg.field : "None";
    M3DG.Field_vals = cfg.field_vals;
    M3DG.backtrack = cfg.backtrack;
    if (cfg.has_backtrack)
        M3DG.timestep = cfg.timestep;
    if (cfg.has_momentum)
        M3DG.momentum = cfg.momentum;
    M3DG.m = cfg.bfgs_saved_states;

    if (cfg.remesher.present)
    {
        sim.Options.max_absolute_length = cfg.remesher.size_max;
        sim.Options.min_absolute_length = cfg.remesher.size_min;
        sim.Options.refine_angle = cfg.remesher.refine_angle;
        sim.Options.aspect_min = cfg.remesher.aspect_min;
        sim.Options.maxIterations = 1;
    }

    std::cout << "The energy elements are ";
    for (const std::string &name : sim.Energies)
        std::cout << name << " ";
    std::cout << "\n";
}

std::string build_descriptive_name(const Simulation &sim, int Nsim)
{
    std::string Directory = "";
    std::stringstream stream;
    auto fixed4 = [&stream](double x)
    {
        stream.str(std::string());
        stream << std::fixed << std::setprecision(4) << x;
        return stream.str();
    };

    size_t bead_counter = 0;
    for (size_t z = 0; z < sim.Energies.size(); z++)
    {
        const std::string &name = sim.Energies[z];
        Directory += name + "_";
        if ((name == "Bead" || name == "H1_Bead" || name == "H2_Bead") && bead_counter < sim.Beads.size())
        {
            const Bead &bead = sim.Beads[bead_counter];
            Directory += "radius_" + fixed4(bead.sigma) + "_";
            Directory += bead.interaction + "_";
            Directory += "str_" + fixed4(bead.strength) + "_";
            if (bead.Constraint == "Radial")
                Directory += "theta_const_" + fixed4(bead.Constraint_constants[0]) + "_";
            bead_counter += 1;
        }
        for (double c : sim.Energy_constants[z])
            Directory += fixed4(c) + "_";
    }

    if (sim.cfg.has_field)
    {
        Directory += sim.M3DG.Field + "_";
        for (double v : sim.M3DG.Field_vals)
            Directory += fixed4(v) + "_";
    }

    bool bonds_exist = false;
    for (const Bead &bead : sim.Beads)
    {
        if (bead.Bond_type.size() > 0 && !bonds_exist)
        {
            Directory += "Bonds_";
            bonds_exist = true;
        }
        for (size_t j = 0; j < bead.Bond_type.size(); j++)
            Directory += bead.Bond_type[j] + "_" + fixed4(bead.Interaction_constants_vector[j][0]) + "_";
    }

    for (const std::string &sw : sim.cfg.switches)
        Directory += "Switch_" + sw + "_";

    Directory += "Nsim_" + std::to_string(Nsim) + "/";
    return Directory;
}

std::string make_numbered_dir(const std::string &first_dir)
{
    mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    for (int dir_counter = 1;; dir_counter++)
    {
        std::string candidate = first_dir + std::to_string(dir_counter) + "/";
        struct stat sb;
        if (stat(candidate.c_str(), &sb) == 0)
            continue;
        if (mkdir(candidate.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0)
            return candidate;
    }
}

std::vector<std::string> open_output_files(const Simulation &sim, const std::string &basic_name, bool append)
{
    std::ios_base::openmode mode = append ? std::ios_base::app : std::ios_base::out;

    std::ofstream Sim_data(basic_name + "Output_data.txt", mode);
    if (!append)
    {
        Sim_data << "time step Volume Area ";
        for (const std::string &name : sim.Energies)
            Sim_data << name << " ";
        Sim_data << " Total_E grad_norm backtrackstep\n";
    }
    Sim_data.close();

    std::vector<std::string> Bead_filenames;
    for (size_t i = 0; i < sim.Beads.size(); i++)
    {
        Bead_filenames.push_back(basic_name + "Bead_" + std::to_string(i) + "_data.txt");
        std::ofstream Bead_data(Bead_filenames[i], mode);
        if (!append)
            Bead_data << "####### This data is taken every " << sim.cfg.save_interval
                      << " steps just like the mesh radius is " << sim.Beads.back().sigma << " \n";
    }

    if (sim.cfg.count_remesh)
    {
        std::ofstream Remeshing_count(basic_name + "Remeshing_count.txt", std::ios_base::app);
        Remeshing_count << "#### timestep remeshing_operations nVertices nEdges nFaces \n";
    }
    return Bead_filenames;
}

void Save_mesh(ManifoldSurfaceMesh *mesh, VertexPositionGeometry *geometry, const std::string &basic_name,
               size_t current_t)
{
    std::ofstream o(basic_name + "membrane_" + std::to_string(current_t) + ".obj");
    o << "#This is a meshfile from a saved state\n";
    for (Vertex v : mesh->vertices())
    {
        Vector3 Pos = geometry->inputVertexPositions[v];
        o << "v " << Pos.x << " " << Pos.y << " " << Pos.z << "\n";
    }
    for (Face f : mesh->faces())
    {
        o << "f";
        for (Vertex v : f.adjacentVertices())
            o << " " << v.getIndex() + 1;
        o << "\n";
    }
}
