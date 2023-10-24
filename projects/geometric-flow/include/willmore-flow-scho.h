#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class WillmoreFlowScho {

  public:
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;
    //constructors 
    WillmoreFlowScho() {}
    WillmoreFlowScho(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);
    virtual SparseMatrix<double> buildFlowOperator(const SparseMatrix<double>& M, double h) const;  
    
    void integrate(double h);
};