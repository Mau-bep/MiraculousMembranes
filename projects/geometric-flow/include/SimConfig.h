#pragma once

// Shared loading of the simulation input file (Input_file.json).
//
//   SimConfig cfg = load_config(path);        // parse + validate the json
//   Simulation sim;
//   build_simulation(cfg, sim);                // mesh, energies, beads, M3DG, E_Handler, remesher
//
// Both main_cluster and main_visualize go through these two calls, so a new
// option only needs to be added here.

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/remeshing.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <nlohmann/json.hpp>

#include "Beads.h"
#include "ConvergenceMonitor.h"
#include "Energy_Handler.h"
#include "Interaction.h"
#include "Mem-3dg.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

struct EnergySpec
{
    std::string name;
    std::vector<double> constants;
};

struct BeadSpec
{
    std::string gradient_order = "Bead"; // name used in the energy list
    std::string constraint = "None";
    std::vector<double> constraint_constants;

    Vector3 pos{0.0, 0.0, 0.0};
    double radius = 0.0;
    std::string state;
    double inter_str = 0.0;
    std::string mem_inter;

    bool has_rc = false;
    double rc = -1.0;
    double outside = 1.0; // Frenkel, Frenkel_Normal_nopush, Linear, Linear_Normal, One_over_r_x
    bool has_shift = false;
    double shift = 0.0;   // LJ; Shifted-LJ computes it from rc unless given
    std::vector<double> z_axis; // Gravity, Pinch

    std::vector<std::string> bonds;
    std::vector<std::vector<double>> bonds_constants;
    std::vector<int> partners; // "Beads": indices of the bonded beads

    bool has_velocity = false;
    Vector3 velocity{0.0, 0.0, 0.0};
    bool has_final_pos = false;
    Vector3 final_pos{0.0, 0.0, 0.0};
};

struct RemeshParams
{
    bool present = false;
    double size_max = 0.0;
    double size_min = 0.0;
    double refine_angle = 0.0;
    double aspect_min = 0.0;
};

// "adapt_remesh": "quality" (see RemeshMonitor.h): remesh when the fraction
// of edges the remesher would change has grown by more than f_tol since the
// last remesh, at most every min_every and at least every max_every steps.
struct RemeshQualityParams
{
    double f_tol = 0.01;
    int min_every = 1;
    int max_every = 100;
};

// "polish": what BFGS-Normal does when the sliver guard finds a corner below
// 0.2 rad. BFGS-Normal only moves vertices along the normals, so remeshing
// and continuing leaves the new vertices on the faceted surface. Instead the
// mesh is remeshed, relaxed with gd_steps of gradient descent (which also
// moves vertices tangentially) and handed to `then`. After max_cycles polishes
// (-1: no limit) the guard only remeshes, as before. gd_steps 0 turns it off.
struct PolishParams
{
    int gd_steps = 100;
    std::string then = "BFGS-Normal"; // or "BFGS"
    int max_cycles = -1;
};

// "remesh_rollback": a remesh (sliver collapses included) that raises any
// single energy term by more than max_rise * max(|term|, 1), or makes one
// NaN/inf, or after which the integrator step fails, is undone, and
// remeshing pauses for `wait` steps. max_rise 0 turns it off. Per term because the total hides it: at a nearly pinched bud neck
// remeshes folded an edge to ~179 deg and Bending_tan went 0.5 -> 1435 and
// 0.5 -> 128 (rises of 1435 and 128 in these units) while the surface tension
// term, 620, made the total rise only 250% and 22%. Measured 2026-10-06, the
// largest rise of a remesh that was fine: 0.63 (refinement, regression d),
// 0.37 (first remesh of a coarse initial mesh), 2.6 (polish at size_min 0.005).
// One broken remesh rose only 9.96 (0.7 -> 10.7): the step after it failed.
struct RemeshRollbackParams
{
    double max_rise = 10.0;
    int wait = 10;
};

struct SimConfig
{
    nlohmann::json raw;      // the parsed file, for anything not covered below
    std::string source_path; // file the config was read from

    // Visualize only: when "Subfolder" is given the real config is
    // first_dir/Subfolder/Input_file.json and outputs are appended there.
    // If that folder does not exist, subfolder_missing is set and this file
    // itself is the config (main_visualize offers to create the folder).
    bool loaded_from_subfolder = false;
    bool subfolder_missing = false;
    std::string subfolder_dir; // with trailing '/'
    // The file that was passed in when its Subfolder replaced it (main_cluster
    // resume reads extra_steps, timesteps and Switches from it); null otherwise.
    nlohmann::json launch_raw;

    std::string init_file;
    std::string first_dir;
    size_t timesteps = 0;
    int save_interval = 1;
    bool finish_sim = false;
    bool saving_states = false;

    bool remeshing = true;
    int remesh_every = 1;
    bool adapt_remesh = true;
    bool quality_remesh = false; // "adapt_remesh": "quality"
    RemeshQualityParams remesh_quality;
    bool count_remesh = false;
    bool remesh_log = false; // Remesh_log.txt: mesh quality and energy around every remesh
    RemeshParams remesher;

    std::string integration = "Gradient_descent";
    int bfgs_saved_states = 10;
    LbfgsOptions lbfgs;
    StoppingParams stopping;
    PolishParams polish;
    RemeshRollbackParams remesh_rollback;
    std::vector<std::string> switches;
    std::vector<int> switch_times;

    Vector3 displacement{0.0, 0.0, 0.0};
    double rescale = 1.0;
    bool has_initial_noise = false;
    double initial_noise = 0.0;
    bool harmonic = false;

    bool recentering = false;
    bool boundary = false;
    bool has_field = false;
    std::string field = "None";
    std::vector<double> field_vals;
    bool has_backtrack = false;
    bool backtrack = true;
    double timestep = 1e-4;
    bool has_momentum = false;
    bool momentum = false;

    std::vector<EnergySpec> energies;
    std::vector<BeadSpec> beads;
};

// Bead membrane interactions make_interaction can build (the "mem_inter" key).
const std::vector<std::string> &known_interactions();

// Names of the switches main_cluster knows how to apply.
const std::vector<std::string> &known_switches();

// Parse and validate. Throws std::runtime_error with a readable message on
// missing required keys. With resolve_subfolder=true a "Subfolder" key makes
// the loader read first_dir/Subfolder/Input_file.json instead (main_visualize),
// or set subfolder_missing if that folder does not exist.
SimConfig load_config(const std::string &path, bool resolve_subfolder = false);

// True when the file has "continue_sim": true (main_cluster then resumes the
// run in first_dir/Subfolder). Throws if the file cannot be read, or if it asks
// to continue a run whose folder it does not name or that does not exist.
bool wants_continue(const std::string &path);

// Everything a running simulation owns. Not copyable or movable: the beads,
// the energy handler and the integrator hold pointers into each other.
struct Simulation
{
    Simulation() = default;
    Simulation(const Simulation &) = delete;
    Simulation &operator=(const Simulation &) = delete;
    ~Simulation();

    SimConfig cfg;

    ManifoldSurfaceMesh *mesh = nullptr;
    VertexPositionGeometry *geometry = nullptr;

    std::vector<std::string> Energies;
    std::vector<std::vector<double>> Energy_constants;
    double V_bar = 0.0; // target volume
    double A_bar = 0.0; // final target area
    double dA = 0.0;    // per-step change of the area target (Area_constraint with nu > 0)

    std::vector<std::unique_ptr<Interaction>> Interaction_container;
    std::vector<Bead> Beads; // reserved up front, never grows after build

    Mem3DG M3DG;
    E_Handler Sim_handler;
    RemeshOptions Options;
};

// A copy of the membrane and of the per-element data that lives across steps
// (energy handler, integrator, remesher options), to undo a remesh. The mesh
// is declared first so the data on it is destroyed before it.
struct MeshSnapshot
{
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;
    VertexData<Vector3> previous_grad, current_grad;
    FaceData<double> face_reference;
    VertexData<double> H_vector_0, dH_vector;
    VertexData<int> integrator_no_remesh_v, options_no_remesh_v;
    EdgeData<int> options_no_remesh;
    std::vector<Vector3> bead_positions; // a failed step may have moved them
};

MeshSnapshot take_snapshot(const Simulation &sim);

// Put the snapshot in place of sim.mesh/sim.geometry everywhere they are held
// (integrator, energy handler, beads) and delete the current ones; the beads
// go back to their positions too. The vertex normals of BFGS-Normal are not
// kept: recompute them as after a remesh.
void restore_snapshot(Simulation &sim, MeshSnapshot &snap);

// The energy terms (sim.Energies order, beads included) and their sum in E;
// the handler's own per-term values are put back.
std::vector<double> energy_terms(Simulation &sim, double &E);

// Index of the first term that is not finite or rose by more than
// max_rise * max(|before|, 1) (RemeshRollbackParams), -1 if none.
int exploded_term(const std::vector<double> &before, const std::vector<double> &after, double max_rise);

// Load the mesh and set up energies, beads, integrator and remesher options.
// With resume_mesh set, the energy targets still come from init_file, but the
// membrane itself is read from resume_mesh (main_cluster resume).
void build_simulation(const SimConfig &cfg, Simulation &sim, const std::string &resume_mesh = "");

// Build the interaction for one bead (throws on an unknown mem_inter).
// params = {inter_str, radius, rc, ...type specific}
std::unique_ptr<Interaction> make_interaction(const BeadSpec &spec, ManifoldSurfaceMesh *mesh,
                                              VertexPositionGeometry *geometry,
                                              std::vector<double> &params);

// Command line of main_cluster/main_visualize: <input.json> [--check].
// A trailing number (the old Nsim argument) is accepted and ignored so that
// existing job scripts keep working. Returns false (after printing the usage)
// if the arguments are not understood.
bool parse_command_line(int argc, char **argv, std::string &config_path, bool &check_only);

// Descriptive run name: energies, constants, beads, bonds, switches.
// Recorded in first_dir/runs_index.txt by record_run().
std::string build_descriptive_name(const Simulation &sim);

// mkdir -p. Throws std::runtime_error if a directory cannot be created.
void make_dirs(const std::string &path);

// First free first_dir/<n>/ (created, together with first_dir if missing).
// Returns the directory with trailing '/'.
std::string make_numbered_dir(const std::string &first_dir);

// Append "<n>/ <descriptive name> <input file>" to first_dir/runs_index.txt,
// so the numbered run folders can be told apart without opening them.
void record_run(const std::string &first_dir, const std::string &run_dir, const std::string &descriptive_name,
                const std::string &config_path);

// Output_data.txt header and one Bead_<i>_data.txt per bead. With append=true
// (continuing a run) existing files are kept and no headers are written.
// Returns the bead file names.
std::vector<std::string> open_output_files(const Simulation &sim, const std::string &basic_name,
                                           bool append);

void Save_mesh(ManifoldSurfaceMesh *mesh, VertexPositionGeometry *geometry,
               const std::string &basic_name, size_t current_t);
