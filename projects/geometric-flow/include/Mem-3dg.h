#pragma once

#include <Eigen/Core>
#include <omp.h>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "Beads.h"
// #include "mean-curvature-flow.h"

// #include <geometrycentral/utilities/eigen_interop_helpers.h>
// #include <geometrycentral/utilities/vector3.h>
using namespace geometrycentral;
using namespace geometrycentral::surface;

class Mem3DG {

  public:
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;
    // VertexData<Vector3> velocity;
    VertexData<double> H_Vector_0;
    VertexData<double> dH_Vector;
    Bead Bead_1;

    // EdgeData<int> No_remesh_list;
    VertexData<int> No_remesh_list_v;

    // VertexData<Vector3> PrevForce;
    // Vector<double> Delta_x;
    // Vector<double> Delta_y;
    // DenseMatrix<double> Quasi_H;
    // DenseMatrix<double> Inv_Quasi_H;

    

    double Old_norm2;
    double Current_norm2;
    size_t system_time;
    double grad_norm;
    bool is_test;
    bool pulling;
    double pulling_force;
  
    // double pulling_offset;
    // constructors
    Mem3DG() {};
    Mem3DG(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);
    Mem3DG(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo,Bead input_Bead);
    Mem3DG(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, bool test);

    virtual VertexData<Vector3> buildFlowOperator(double h, double V_bar, double nu, double c0, double P0,double KA,double KB,double Kd);
    
    Vector3 computeHalfedgeMeanCurvatureVector(Halfedge he) const;
    Vector3 computeHalfedgeGaussianCurvatureVector(Halfedge he) const ;
    Vector3 dihedralAngleGradient(Halfedge he, Vertex v) const;
    
    // virtual Vector3 cornerAngleGradient(Corner c, Vertex v) const;
    virtual VertexData<Vector3> OsmoticPressure(double D_P) const;
    virtual VertexData<Vector3> SurfaceTension(double lambda) const;
    virtual VertexData<Vector3> SurfaceGrad() const;
    
    virtual VertexData<Vector3> Bending(double H0,double KB) const;
    virtual double E_Pressure(double P0,double V ,double V_bar) const;
    virtual double E_Surface(double KA,double A, double A_bar) const;
    virtual double E_Bending(double H0, double KB) const;
    
    virtual void Grad_Vol_dx(std::ofstream& Gradient_file,double P0, double V_bar,size_t index) const;
    void Grad_Bead_dx(std::ofstream& Gradient_file,bool Save);

    virtual VertexData<Vector3> Grad_Vol(std::ofstream& Gradient_file,double P0, double V_bar,bool Save) const;
    virtual VertexData<Vector3> Grad_Area(std::ofstream& Gradient_file, double A_bar,double KA, bool Save) const;
    virtual VertexData<Vector3> Grad_Bending(std::ofstream& Gradient_file, double H_bar,double KB,bool Save);
    virtual void Grad_Bending_2(std::ofstream& Gradient_file, double H_bar,double KB);

    virtual void Bending_test(std::ofstream& Analysis_file, double H_bar,double KB);
    VertexData<Vector3> Grad_Bead(std::ofstream& Gradient_file,bool Save,bool Projection);
    virtual VertexData<Vector3> Grad_tot_Area(std::ofstream& Gradient_file, bool Save) const;
    

    double Backtracking(VertexData<Vector3> Force,double D_P,double V_bar,double A_bar,double KA,double KB,double H_bar,bool bead) ;
    double Backtracking(VertexData<Vector3> Force,double D_P,double V_bar,double A_bar,double KA,double KB,double H_bar) ;
    virtual VertexData<Vector3> Project_force(VertexData<Vector3> Force) const; 
    virtual bool Area_sanity_check();  
    
    double integrate(double h, double V_bar, double nu, double c0,double P0,double KA,double KB, double Kd ,std::ofstream& Sim_data, double time, bool bead,std::ofstream& Bead_data,bool Save_output_data);
    
    double integrate(double h, double V_bar, double nu, double c0,double P0,double KA,double KB, double Kd ,std::ofstream& Sim_data, double time,bool Save);
    double integrate_finite(double h, double V_bar, double nu, double c0,double P0,double KA,double KB, double Kd,std::ofstream& Sim_data, double time,std::ofstream& Gradient_file_vol,std::ofstream& Gradient_file_area,std::ofstream& Gradient_file_bending,bool Save );
};