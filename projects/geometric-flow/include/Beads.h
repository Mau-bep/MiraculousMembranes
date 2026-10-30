#pragma once

#include <Eigen/Core>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
// Bead only stores an Interaction pointer; the .cpp files include Interaction.h
class Interaction;

using namespace geometrycentral;
using namespace geometrycentral::surface;

class Bead
{

public:
  int Bead_id = 0;
  int Total_beads = 0;
  Vector3 Pos{0.0, 0.0, 0.0};
  Vector3 Prev_Total_force{0.0, 0.0, 0.0};
  Vector3 Total_force{0.0, 0.0, 0.0};
  Vector3 CoverageForce{0.0, 0.0, 0.0};
  double sigma = 0.0;
  double strength = 0.0;
  // +1 the bead is outside the membrane (default), -1 inside the closed membrane ("outside": -1 in the input file):
  // the faces that the Coverage energy counts are the ones whose outward normal points to the bead (+1) or away from it (-1)
  double side = 1.0;

  double rc = -1.0;
  double prev_force = 0.0;
  double prev_E_stationary = 0.0;

  std::string state;
  Vector3 Velocity{0.0, 0.0, 0.0};
  Vector3 FinalPos{0.0, 0.0, 0.0};

  Vector3 Stopping_pos{0.0, 0.0, 0.0};

  ManifoldSurfaceMesh *mesh = nullptr;
  VertexPositionGeometry *geometry = nullptr;
  std::string interaction;
  std::vector<Bead *> Beads;
  std::vector<std::string> Bond_type;
  // std::vector<double> Interaction_constants;
  std::vector<std::vector<double>> Interaction_constants_vector;

  std::string Constraint = "None";
  std::vector<double> Constraint_constants;
  // Constraint "Normal" (planar membrane): in BFGS-Normal the bead only moves
  // along this unit vector, the normal of the plane. Set from Constraint_constants
  // [nx, ny, nz], or from the mean normal of the interior vertices when they are
  // empty. The constraint holds only while Axis_active, which integrate_BFGS_Normal
  // sets for its step: every other integrator moves the bead freely.
  Vector3 Constraint_axis{0.0, 0.0, 0.0};
  bool Axis_active = false;

  Interaction *Bead_I = nullptr;

  // constructors
  Bead() {};
  Bead(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, Vector3 Position, double sigma_bead, double strg);
  Bead(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, Vector3 Position, double input_sigma, double strg, double input_rc);
  Bead(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, Vector3 Position, std::vector<double> Energy_constants, Interaction *Interact, int id, int Number_beads);
  void Reasign_mesh(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo);

  VertexData<Vector3> Gradient();
  void Reset_bead(Vector3 Actual_pos);
  void Move_bead(double dt, Vector3 center);
  void Move_bead(double dt, Vector3 center, Vector3 Force);

  void Add_bead(Bead *bead, std::string Interaction, std::vector<double> Interaction_strength);

  // State "rigid": the bead keeps the length of each of its "Rigid" bonds
  // (bonds_constants [L]) and only moves perpendicular to them
  bool Has_rigid_bond() const;
  // A free bead with Constraint "Normal" inside BFGS-Normal: one coordinate along Constraint_axis, a variable of the integrator
  bool Moves_along_axis() const { return Axis_active && state == "default" && Constraint == "Normal"; }
  Vector3 Free_force(Vector3 Force) const;
  double Enforce_rigid_bonds();

  void update_state();
};

// Puts every rigid bond back at its length
void Enforce_rigid_bonds(const std::vector<Bead *> &beads);