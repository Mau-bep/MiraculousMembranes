#pragma once

#include <Eigen/Core>
#include <omp.h>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
// #include <geometrycentral/utilities/eigen_interop_helpers.h>
// #include <geometrycentral/utilities/vector3.h>
using namespace geometrycentral;
using namespace geometrycentral::surface;

class Mem3DG_p {

  public:
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;
    // VertexData<Vector3> velocity;
    // VertexData<double> H_Vector_0;
    // VertexData<double> dH_Vector;
    double Old_norm2;
    double Current_norm2;
    double system_time;
    double grad_norm;
    // constructors
    Mem3DG_p() {};
    Mem3DG_p(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);


    virtual VertexData<Vector3> buildFlowOperator(double h, double V_bar, double nu, double c0, double P0,double KA,double KB,double Kd) const;
    virtual Vector3 computeHalfedgeMeanCurvatureVector(Halfedge he) const;
    virtual Vector3 computeHalfedgeGaussianCurvatureVector(Halfedge he) const ;
    virtual Vector3 dihedralAngleGradient(Halfedge he, Vertex v) const;
    
    virtual VertexData<Vector3> OsmoticPressure(double D_P) const;
    virtual VertexData<Vector3> SurfaceTension(double lambda) const;
    virtual VertexData<Vector3> Bending(double H0,double KB) const;
    virtual double E_Pressure(double P0,double V ,double V_bar) const;
    virtual double E_Surface(double KA,double A, double A_bar) const;
    virtual double E_Bending(double H0, double KB) const;
    
    double Backtracking(VertexData<Vector3> Force,double D_P,double V_bar,double A_bar,double KA,double KB,double H_bar);
    double integrate(double h, double V_bar, double nu, double c0,double P0,double KA,double KB, double Kd ,std::ofstream& Sim_data, double time);
};