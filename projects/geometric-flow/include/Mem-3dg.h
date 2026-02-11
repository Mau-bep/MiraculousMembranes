#pragma once

#include <Eigen/Core>
#include <omp.h>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "Beads.h"
#include "Energy_Handler.h"

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
    std::vector<Bead *> Beads;
    // EdgeData<int> No_remesh_list;
    VertexData<int> No_remesh_list_v;

    // VertexData<Vector3> PrevForce;
    // Vector<double> Delta_x;
    // Vector<double> Delta_y;
    // DenseMatrix<double> Quasi_H;
    // DenseMatrix<double> Inv_Quasi_H;
    std::string basic_name;
    
    double Area_evol_steps;
    double Old_norm2;
    double Current_norm2;
    size_t system_time;
    double grad_norm;
    bool is_test;
    bool pulling;
    bool Save_SS;
    bool stop_increasing;
    double pulling_force;
    bool small_TS;
    bool recentering;
    bool boundary;
    double A;
    double V;
    double E_Vol;
    double E_Sur;
    double E_Ben;
    double E_Bead;

    bool backtrack;
    double timestep;

    std::string Field;
    std::vector<double> Field_vals;


    E_Handler* Sim_handler;
    std::vector<double> Energy_vals;

    Vector3 Total_force;

    // Eigen::SimplicialLDLT<SparseMatrix<double>> solver_H;
    bool momentum;
    double learn_rate;

    // double pulling_offset;
    // constructors
    Mem3DG() {};
    Mem3DG(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);
    Mem3DG(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo,Bead input_Bead);
    Mem3DG(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, bool test);

    void Add_bead(Bead *bead);

    virtual VertexData<Vector3> buildFlowOperator(double h, double V_bar, double nu, double c0, double P0,double KA,double KB,double Kd);
    virtual VertexData<Vector3> buildFlowOperator(double V_bar, double P0,double KA,double KB,double h);
    virtual VertexData<Vector3> buildFlowOperator(double h,double V_bar,double P0,double KA);
    Vector3 computeHalfedgeMeanCurvatureVector(Halfedge he) const;
    Vector3 computeHalfedgeGaussianCurvatureVector(Halfedge he) const ;
    Vector3 dihedralAngleGradient(Halfedge he, Vertex v) const;
    
    // virtual Vector3 cornerAngleGradient(Corner c, Vertex v) const;
    virtual VertexData<Vector3> OsmoticPressure() const;
    virtual VertexData<Vector3> SurfaceTension() const;
    virtual VertexData<Vector3> SurfaceGrad() const;
    SparseMatrix<double> H2_operator(bool CM, bool Vol_const, bool Area_const);
    SparseMatrix<double> H1_operator(bool CM, bool Vol_const, bool Area_const);
    
    virtual VertexData<Vector3> Bending(double H0) const;
    virtual double E_Volume_constraint(double P0,double V ,double V_bar) const;
    virtual double E_Area_constraint(double KA,double A, double A_bar) const;
    virtual double E_Pressure(double P0,double V ,double V_bar) const;
    virtual double E_Surface(double KA,double A, double A_bar) const;
    virtual double E_Bending(double H0, double KB) const;
    
    virtual VertexData<Vector3> Linear_force_field(double x0,double slope) const;


    virtual void Grad_Vol_dx(std::ofstream& Gradient_file,double P0, double V_bar,size_t index) const;
    void Grad_Bead_dx(std::ofstream& Gradient_file,bool Save);

    virtual VertexData<Vector3> Grad_Vol(std::ofstream& Gradient_file,double P0, double V_bar,bool Save) const;
    virtual VertexData<Vector3> Grad_Area(std::ofstream& Gradient_file, double A_bar,double KA, bool Save) const;
    virtual VertexData<Vector3> Grad_Bending(std::ofstream& Gradient_file, double H_bar,double KB,bool Save);
    virtual void Grad_Bending_2(std::ofstream& Gradient_file, double H_bar,double KB);

    virtual void Bending_test(std::ofstream& Analysis_file, double H_bar,double KB);
    VertexData<Vector3> Grad_Bead(std::ofstream& Gradient_file,bool Save,bool Projection);
    virtual VertexData<Vector3> Grad_tot_Area(std::ofstream& Gradient_file, bool Save) const;
    
    void Smooth_vertices();

    double Backtracking();
    double Backtracking_grad(Eigen::VectorXd pk, double Projection, double Current_grad_norm);
    double Backtracking_BFGS(VertexData<Vector3> Force);
    double Backtracking(VertexData<Vector3> Force,double P0,double V_bar,double A_bar,double KA,double KB,double H_bar,bool bead, bool pulling) ;
    double Backtracking(VertexData<Vector3> Force,double D_P,double V_bar,double A_bar,double KA,double KB,double H_bar) ;
    double Backtracking(VertexData<Vector3> Force,double D_P,double V_bar,double KA);
    double Backtracking_field(VertexData<Vector3> Force,double D_P,double V_bar,double A_bar,double KA,double KB, double H_bar);
    virtual VertexData<Vector3> Project_force(VertexData<Vector3> Force) const; 
    virtual bool Area_sanity_check();  
    VertexData<double> Vert_sizing(FaceData<double> Face_sizings);
    FaceData<double> Face_sizings();
    Eigen::Matrix2d Face_sizing(Face f);
    EdgeData<double> Edge_sizing(VertexData<double> Vert_sizings);

    double integrate(std::ofstream& Sim_data , double time, std::vector<std::string>Bead_data_filenames, bool Save_output_data);
    Eigen::MatrixXd integrate_BFGS(std::ofstream& Sim_data , double time, std::vector<std::string>Bead_data_filenames, bool Save_output_data, Eigen::MatrixXd Hessian);
    double integrate_Newton(std::ofstream& Sim_data , double time, std::vector<std::string>Bead_data_filenames, bool Save_output_data, std::vector<std::string> Constraints,std::vector<std::string> Data_filenames);
    
    double integrate_implicit(std::vector<std::string> Energies,  std::vector<std::vector<double>> Energy_constants,std::ofstream& Sim_data , double time, std::vector<std::string>Bead_data_filenames, bool Save_output_data);
    
    void Get_Energies(std::vector<std::string>Energies, std::vector<std::vector<double>> Energy_constants, double* NewE);

    double integrate(double h, double V_bar, double nu, double c0,double P0,double KA,double KB, double Kd ,std::ofstream& Sim_data, double time, bool bead,std::vector<std::string> Bead_data_filenames,bool Save_output_data,bool pulling);
    double integrate(double h, double V_bar,double P0,double KA,std::ofstream& Sim_data, double time,bool Save);
    double integrate(double h, double V_bar, double nu, double c0,double P0,double KA,double KB, double Kd ,std::ofstream& Sim_data, double time,bool Save);
    double integrate_field(double h, double V_bar, double nu, double P0,double KA,double KB, double slope,double x0 ,std::ofstream& Sim_data, double time,bool Save);
    
    double integrate_finite(double h, double V_bar, double nu, double c0,double P0,double KA,double KB, double Kd,std::ofstream& Sim_data, double time,std::ofstream& Gradient_file_vol,std::ofstream& Gradient_file_area,std::ofstream& Gradient_file_bending,bool Save );

    void Save_mesh(size_t current_t);
};