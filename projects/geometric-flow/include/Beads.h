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
    
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;
    std::string interaction;

    // constructors
    Bead() {};
    Bead(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, Vector3 Position, double sigma_bead, double strg );


    VertexData<Vector3> Gradient();
    double Energy();
    void Reset_bead(Vector3 Actual_pos);
    void Move_bead(double dt, Vector3 center);

};