#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class NormalFlow {

  public:
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;
    //constructors 
    NormalFlow() {}
    NormalFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);
    virtual SparseMatrix<double> buildFlowOperator(const SparseMatrix<double>& M, double h) const;  
    virtual SparseMatrix<double> buildFlowOperatorx(const SparseMatrix<double>& M, double h) const;
    virtual SparseMatrix<double> buildFlowOperatory(const SparseMatrix<double>& M, double h) const;
    virtual SparseMatrix<double> buildFlowOperatorz(const SparseMatrix<double>& M, double h) const;

    void integrate(double h);
};