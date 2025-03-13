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
    double integrate_Sob(double h,double Kb, std::ofstream& Sim_data,bool Save_output_data);

    Vector3 computeHalfedgeMeanCurvatureVector(Halfedge he) const;
    Vector3 computeHalfedgeGaussianCurvatureVector(Halfedge he) const ;
    Vector3 dihedralAngleGradient(Halfedge he, Vertex v) const;
    virtual VertexData<Vector3> Bending(double H0) const;

    VertexData<Vector3> Sobolev_operator() const;
    VertexData<Vector3> SurfaceGrad() const;
    VertexData<Vector3> OsmoticPressure() const;

    virtual double E_Bending(double H0, double KB) const;

    double Backtracking(VertexData<Vector3> Force, double h );
};