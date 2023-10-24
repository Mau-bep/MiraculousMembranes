#pragma once
#include <Eigen/Core>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class IMem3DG {

  public:
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;
    
    size_t index;
    // constructors
    IMem3DG() {};
    IMem3DG(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);

    virtual SparseMatrix<double> buildFlowOperator(const SparseMatrix<double>& M,double h,double nu,double V_bar,double P0, double KA,double KB) const;

    void integrate(double h,double nu,double V_bar,double P0, double KA,double KB);
};