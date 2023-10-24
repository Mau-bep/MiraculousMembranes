#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include <Eigen/Core>
#include <Eigen/Cholesky>
using namespace geometrycentral;
using namespace geometrycentral::surface;

class MembraneFlow {

  public:
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;

    // constructors
    MembraneFlow() {};
    MembraneFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);

    virtual SparseMatrix<double> buildFlowOperator(const SparseMatrix<double>& M,double sigma, double kappa, double H0,double h) const;

    void integrate(double h, double sigma , double kappa, double H0, double P, double V0);
};