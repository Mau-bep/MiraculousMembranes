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

// #include "mean-curvature-flow.h"

// #include <geometrycentral/utilities/eigen_interop_helpers.h>
// #include <geometrycentral/utilities/vector3.h>
using namespace geometrycentral;
using namespace geometrycentral::surface;

enum class RemeshBoundaryCondition2
{
  Fixed,
  Tangential,
  Free
};
enum class RemeshSmoothStyle2
{
  Circumcentric,
  Laplacian
};
struct RemeshOptions2
{
  double targetEdgeLength = -1;    // the target edge length in flat regions. If `targetEdgeLength` is negative, the target
                                   // edge length is set to relative the input mesh's mean edge length
  size_t maxIterations = 10;       // the maximum number of iterations to run for
  double curvatureAdaptation = 0;  // how much target length should vary due to curvature. Set curvatureAdaptation
                                   // to 0 if you want lengths to be approximately targetEdgeLength everywhere
  double minRelativeLength = 0.05; // the minimum possible edge length allowed in the output mesh. Defined relative to
                                   // targetEdgeLength
  double min_absolute_length =
      0.001; // the minimum possible edge length allowed in the output mesh, as an absolute number
  double max_absolute_length =
      0.2;                   // the maximum possible edge length allowed in the output mesh, as an absolute number
  double refine_angle = 0.7; // THe maximum dihedral angle allowed in the output mesh, in radians
  double aspect_min = 0.2;
  bool no_remesh_list = false;
  int numberOp = 0;
  float angleThresh = 0.15;

  std::vector<Edge> Remesh_list_e;   // list of edges to remesh. If empty, all edges are considered for remeshing
  std::vector<Vertex> Remesh_list_v; // list of vertices to remesh. If empty, all vertices are considered for remeshing
  std::vector<Face> Remesh_list_f;   // list of faces to remesh. If empty, all faces are considered for remeshing

  EdgeData<int> No_remesh_list;
  VertexData<int> No_remesh_list_v;

  RemeshSmoothStyle2 smoothStyle = RemeshSmoothStyle2::Circumcentric; // smoothing function to use
  RemeshBoundaryCondition2 boundaryCondition =
      RemeshBoundaryCondition2::Tangential; // allowed movement of boundary vertices
};

class Mem3DG
{

public:
  ManifoldSurfaceMesh *mesh;
  VertexPositionGeometry *geometry;
  // VertexData<Vector3> velocity;
  VertexData<double> H_Vector_0;
  VertexData<double> dH_Vector;
  Bead Bead_1;
  std::vector<Bead *> Beads;
  // EdgeData<int> No_remesh_list;
  VertexData<int> No_remesh_list_v;

  // VertexData<Vector3> PrevForce;
  // Vector<double> Delta_x;
  // Vector<double> Delta_y;
  // DenseMatrix<double> Quasi_H;
  // DenseMatrix<double> Inv_Quasi_H;
  std::string basic_name;

  double Area_evol_steps;
  double Old_norm2;
  double Current_norm2;
  size_t system_time;
  double grad_norm;
  bool is_test;
  bool pulling;
  bool Save_SS;
  bool stop_increasing;
  double pulling_force;
  bool small_TS;
  bool recentering;
  bool boundary;
  double A;
  double V;
  double E_Vol;
  double E_Sur;
  double E_Ben;
  double E_Bead;

  bool backtrack;
  size_t discreteTs;
  double timestep;
  bool remesh_flag;
  std::string Field;
  std::vector<double> Field_vals;

  E_Handler *Sim_handler;
  std::vector<double> Energy_vals;

  Vector3 Total_force;

  // BFGS PARAMETERS
  int BFGS_iter;
  int m;
  std::vector<Eigen::VectorXd> s_list;
  std::vector<Eigen::VectorXd> y_list;
  std::vector<double> rho_list;

  // Newton param
  int Newton_iter;
  double mu;
  // Convergence estimators
  int N_data;
  double mean_E;
  double var_E;
  double mean_Grad;
  double var_Grad;
  bool Turn_normal_iter;

  // Eigen::SimplicialLDLT<SparseMatrix<double>> solver_H;
  bool momentum;
  double learn_rate;

  // double pulling_offset;
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

  virtual VertexData<Vector3> buildFlowOperator(double h, double V_bar, double nu, double c0, double P0, double KA, double KB, double Kd);
  virtual VertexData<Vector3> buildFlowOperator(double V_bar, double P0, double KA, double KB, double h);
  virtual VertexData<Vector3> buildFlowOperator(double h, double V_bar, double P0, double KA);
  Vector3 computeHalfedgeMeanCurvatureVector(Halfedge he) const;
  Vector3 computeHalfedgeGaussianCurvatureVector(Halfedge he) const;
  Vector3 dihedralAngleGradient(Halfedge he, Vertex v) const;

  // virtual Vector3 cornerAngleGradient(Corner c, Vertex v) const;
  virtual VertexData<Vector3> OsmoticPressure() const;
  virtual VertexData<Vector3> SurfaceTension() const;
  virtual VertexData<Vector3> SurfaceGrad() const;
  SparseMatrix<double> H2_operator(bool CM, bool Vol_const, bool Area_const);
  SparseMatrix<double> H1_operator(bool CM, bool Vol_const, bool Area_const);

  virtual VertexData<Vector3> Bending(double H0) const;
  virtual double E_Volume_constraint(double P0, double V, double V_bar) const;
  virtual double E_Area_constraint(double KA, double A, double A_bar) const;
  virtual double E_Pressure(double P0, double V, double V_bar) const;
  virtual double E_Surface(double KA, double A, double A_bar) const;
  virtual double E_Bending(double H0, double KB) const;

  virtual VertexData<Vector3> Linear_force_field(double x0, double slope) const;

  virtual void Grad_Vol_dx(std::ofstream &Gradient_file, double P0, double V_bar, size_t index) const;
  void Grad_Bead_dx(std::ofstream &Gradient_file, bool Save);

  virtual VertexData<Vector3> Grad_Vol(std::ofstream &Gradient_file, double P0, double V_bar, bool Save) const;
  virtual VertexData<Vector3> Grad_Area(std::ofstream &Gradient_file, double A_bar, double KA, bool Save) const;
  virtual VertexData<Vector3> Grad_Bending(std::ofstream &Gradient_file, double H_bar, double KB, bool Save);
  virtual void Grad_Bending_2(std::ofstream &Gradient_file, double H_bar, double KB);

  virtual void Bending_test(std::ofstream &Analysis_file, double H_bar, double KB);
  VertexData<Vector3> Grad_Bead(std::ofstream &Gradient_file, bool Save, bool Projection);
  virtual VertexData<Vector3> Grad_tot_Area(std::ofstream &Gradient_file, bool Save) const;

  void Smooth_vertices();

  double Backtracking();
  double Backtracking_grad(Eigen::VectorXd pk, double Projection, double Current_grad_norm);
  double Backtracking_grad_Normal(Eigen::VectorXd pk, double Projection, double Current_grad_norm);
  double Backtracking_grad_Normal_2(Eigen::VectorXd pk, double Projection, double Current_grad_norm);

  double Backtracking_BFGS(VertexData<Vector3> Force, std::vector<Vector3> Bead_forces);
  double Backtracking(VertexData<Vector3> Force, double P0, double V_bar, double A_bar, double KA, double KB, double H_bar, bool bead, bool pulling);
  double Backtracking(VertexData<Vector3> Force, double D_P, double V_bar, double A_bar, double KA, double KB, double H_bar);
  double Backtracking(VertexData<Vector3> Force, double D_P, double V_bar, double KA);
  double Backtracking_field(VertexData<Vector3> Force, double D_P, double V_bar, double A_bar, double KA, double KB, double H_bar);
  virtual VertexData<Vector3> Project_force(VertexData<Vector3> Force) const;
  virtual bool Area_sanity_check();
  VertexData<double> Vert_sizing(FaceData<double> Face_sizings);
  FaceData<double> Face_sizings();
  Eigen::Matrix2d Face_sizing(Face f);
  EdgeData<double> Edge_sizing(VertexData<double> Vert_sizings);

  double integrate(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data);
  double integrate_BFGS(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data);
  double integrate_BFGS_Normal(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data);

  double integrate_Newton(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data, std::vector<std::string> Constraints, std::vector<std::string> Data_filenames);
  double integrate_Newton_Normal(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data, std::vector<std::string> Constraints, std::vector<std::string> Data_filenames);
  double integrate_Newton_Normal_Sherman(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data, std::vector<std::string> Constraints, std::vector<std::string> Data_filenames);

  VertexData<Vector3> Newton_Normal_step(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data, std::vector<std::string> Constraints, std::vector<std::string> Data_filenames);
  double integrate_implicit(std::vector<std::string> Energies, std::vector<std::vector<double>> Energy_constants, std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data);

  void Get_Energies(std::vector<std::string> Energies, std::vector<std::vector<double>> Energy_constants, double *NewE);

  double integrate(double h, double V_bar, double nu, double c0, double P0, double KA, double KB, double Kd, std::ofstream &Sim_data, double time, bool bead, std::vector<std::string> Bead_data_filenames, bool Save_output_data, bool pulling);
  double integrate(double h, double V_bar, double P0, double KA, std::ofstream &Sim_data, double time, bool Save);
  double integrate(double h, double V_bar, double nu, double c0, double P0, double KA, double KB, double Kd, std::ofstream &Sim_data, double time, bool Save);
  double integrate_field(double h, double V_bar, double nu, double P0, double KA, double KB, double slope, double x0, std::ofstream &Sim_data, double time, bool Save);

  double integrate_finite(double h, double V_bar, double nu, double c0, double P0, double KA, double KB, double Kd, std::ofstream &Sim_data, double time, std::ofstream &Gradient_file_vol, std::ofstream &Gradient_file_area, std::ofstream &Gradient_file_bending, bool Save);

  void Save_mesh(size_t current_t);

  int remesh(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, RemeshOptions options);
  int remesh(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, MutationManager &mm,
             RemeshOptions options);

  bool splitWorstEdges(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, MutationManager &mm,
                       RemeshOptions options);
  bool improveFaces(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, MutationManager &mm, RemeshOptions options);

  size_t fixDelaunay(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom);
  size_t fixDelaunay(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, MutationManager &mm);
  double smoothByLaplacian(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, MutationManager &mm,
                           double stepSize = 1, RemeshBoundaryCondition2 bc = RemeshBoundaryCondition2::Tangential);

  // Average positions of vertices based on surrounding triangle circumenters as in [Chen & Holst 2011]
  // Returns the average amount each vertex was moved by
  double smoothByCircumcenter(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, MutationManager &mm,
                              double stepSize = 1, RemeshBoundaryCondition2 bc = RemeshBoundaryCondition2::Tangential);
};