#pragma once

#include <Eigen/Core>
#include <omp.h>
#include "geometrycentral/surface/barycentric_coordinate_helpers.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/mutation_manager.h"
#include "geometrycentral/surface/remeshing.h"

#include "Beads.h"
#include <deque>
#include "Energy_Handler.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class Mem3DG
{

public:
  ManifoldSurfaceMesh *mesh = nullptr;
  VertexPositionGeometry *geometry = nullptr;
  E_Handler *Sim_handler = nullptr; // energies and gradients (not owned)
  std::vector<Bead *> Beads;        // external agents moved by the integrators (not owned)
  Bead Bead_1;                      // used by the test programs' constructor only
  VertexData<double> H_Vector_0;
  VertexData<double> dH_Vector;
  VertexData<int> No_remesh_list_v;
  std::string basic_name; // output folder, for the line search logs

  // Integration settings
  bool backtrack = true;
  double timestep = 1e-4; // fixed step when backtrack is false
  bool momentum = false;
  double learn_rate = 0.5;
  bool recentering = true;
  bool boundary = false;
  std::string Field = "None";
  std::vector<double> Field_vals;

  // State
  size_t discreteTs = 0;
  size_t system_time = 0;
  bool small_TS = false;
  bool remesh_flag = false;
  double grad_norm = 0.0;
  double Old_norm2 = 0;
  double Current_norm2 = 0;
  double A = 0.0;
  double V = 0.0;
  double E_Vol = 0.0;
  double E_Sur = 0.0;
  double E_Ben = 0.0;
  double E_Bead = 0.0;
  std::vector<double> Energy_vals;
  Vector3 Total_force{0.0, 0.0, 0.0};

  // Legacy flags still read by the test programs
  double Area_evol_steps = 100000;
  bool is_test = false;
  bool pulling = false;
  bool Save_SS = false;
  bool stop_increasing = false;
  double pulling_force = 0.0;

  // L-BFGS
  int BFGS_iter = 0;
  int m = 10;
  std::vector<Eigen::VectorXd> s_list;
  std::vector<Eigen::VectorXd> y_list;
  std::vector<double> rho_list;

  // Newton
  int Newton_iter = 0;
  double mu = 0.0;

  // Convergence estimators
  int N_data = 0;
  double mean_E = 0;
  double var_E = 0;
  double mean_Grad = 0;
  double var_Grad = 0;
  bool Turn_normal_iter = false;

  // constructors
  Mem3DG() {};
  Mem3DG(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo);
  Mem3DG(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, Bead input_Bead);
  Mem3DG(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, bool test);

  // Rule of five: declare destructor and disable copying to avoid
  // accidental shallow copies of raw pointers held by this class.
  ~Mem3DG();
  Mem3DG(const Mem3DG &) = delete;
  Mem3DG &operator=(const Mem3DG &) = delete;
  Mem3DG(Mem3DG &&) = default;
  Mem3DG &operator=(Mem3DG &&) = default;

  void Add_bead(Bead *bead);

  // Integrators: each takes one step and returns the time step used (negative on failure).
  double integrate(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data);
  double integrate_BFGS(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data);
  double integrate_BFGS_Normal(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data);
  double integrate_Newton(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data, std::vector<std::string> Constraints, std::vector<std::string> Data_filenames);
  double integrate_Newton_Normal(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data, std::vector<std::string> Constraints, std::vector<std::string> Data_filenames);
  // Newton step in the normal direction without moving the mesh (main_visualize)
  VertexData<Vector3> Newton_Normal_step(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data, std::vector<std::string> Constraints, std::vector<std::string> Data_filenames);

  // Line searches
  double Backtracking();
  double Backtracking_grad(Eigen::VectorXd pk, double Projection, double Current_grad_norm);
  double Backtracking_grad_Normal(Eigen::VectorXd pk, double Projection, double Current_grad_norm);
  double Backtracking_BFGS(VertexData<Vector3> Force, std::vector<Vector3> Bead_forces);

  // Remeshing sizing field (diagnostics and main_visualize)
  Eigen::Matrix2d Face_sizing(Face f);
  FaceData<double> Face_sizings();
  VertexData<double> Vert_sizing(FaceData<double> Face_sizings);

  // Finite difference checks used by the test programs
  virtual VertexData<Vector3> SurfaceGrad() const;
  VertexData<Vector3> Grad_Bead(std::ofstream &Gradient_file, bool Save, bool Projection);
  virtual VertexData<Vector3> Grad_tot_Area(std::ofstream &Gradient_file, bool Save) const;
  virtual bool Area_sanity_check();
};