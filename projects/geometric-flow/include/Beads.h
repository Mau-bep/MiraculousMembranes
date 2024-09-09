#pragma once

#include <Eigen/Core>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
// #include "Mem-3dg.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
 
class Bead {

  public:
    Vector3 Pos;
    Vector3 Total_force;
    double sigma;
    double strength;
    double pulling_speed;
    double rc;
    double prev_force;
    double prev_E_stationary;

    std::string state;
    Vector3 Velocity;

    Vector3 Stopping_pos;

    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;
    std::string interaction;
    std::vector<Bead *> Beads;
    std::vector<std::string>  Bond_type;
    // std::vector<double> Interaction_constants;
    std::vector<std::vector<double>> Interaction_constants_vector;

    // constructors
    Bead() {};
    Bead(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, Vector3 Position, double sigma_bead, double strg );
    Bead(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo,Vector3 Position,double input_sigma , double strg, double input_rc);

    void Reasign_mesh(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);

    VertexData<Vector3> Gradient();
    void Set_Force(Vector3 Force);
    double Energy();
    void Reset_bead(Vector3 Actual_pos);
    void Move_bead(double dt, Vector3 center);

    void Add_bead(Bead *bead, std::string Interaction, std::vector<double> Interaction_strength);

    void Bead_interactions();
    

};