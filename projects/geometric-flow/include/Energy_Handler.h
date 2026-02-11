#pragma once

#include <Eigen/Core>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "Beads.h"
// #include "Mem-3dg.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
 

class E_Handler {
    public: 
            
        ManifoldSurfaceMesh* mesh;
        VertexPositionGeometry* geometry;

        std::vector<std::string> Energies;
        std::vector<std::vector<double>> Energy_constants;
        std::vector<std::string> Constraints;
        Eigen::VectorXd Lagrange_mult;
        int N_constraints;
        double Trgt_vol;
        double Trgt_area;
        bool boundary;
        std::vector<Bead *> Beads;
        
        std::vector<double> Energy_values;
        
        VertexData<Vector3> Previous_grad;
        VertexData<Vector3> Current_grad;
        FaceData<double> Face_reference;

        Eigen::MatrixXd Jacobian_constraints;

        std::vector<double> Gradient_norms;

        // Constructors
        E_Handler(){};
        E_Handler(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);
        E_Handler(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, std::vector<std::string> inputEnergyNames, std::vector<std::vector<double>> inputEnergyConstants);
        
        void Add_Bead(Bead* bead);

        void Add_Energy(std::string Energy_name, std::vector<double> Constants);

        
        SparseMatrix<double> H1_operator(bool CM, bool Vol_const, bool Area_const);
        SparseMatrix<double> H2_operator(bool CM, bool Vol_const, bool Area_const);
        
        void update_face_reference();

        virtual double E_Volume_constraint(std::vector<double> Constants) const;
        virtual double E_Area_constraint(std::vector<double> Constants) const;
        virtual double E_SurfaceTension(std::vector<double> Constants) const;
        virtual double E_Bending(std::vector<double> Constants) const;
        virtual double E_Bending_2(std::vector<double> Constants) const;
        virtual double E_Laplace(std::vector<double> Constants) const;
        virtual double E_Edge_reg(std::vector<double> Constants) const;
        virtual double E_Edge_reg_2(std::vector<double> Constants) const;
        virtual double E_Face_reg(std::vector<double> Constants) const;

        virtual VertexData<Vector3> F_Volume_constraint(std::vector<double> Constants) const;
        virtual VertexData<Vector3> F_SurfaceTension(std::vector<double> Constants) const;
        virtual VertexData<Vector3> F_Area_constraint(std::vector<double> Constants) const;
        virtual VertexData<Vector3> F_Bending(std::vector<double> Constants) const;
        virtual VertexData<Vector3> F_SurfaceTension_2(std::vector<double> Constants) const;
        virtual VertexData<Vector3> F_Bending_2(std::vector<double> Constants) const;
        virtual VertexData<Vector3> F_Volume_constraint_2(std::vector<double> Constants) const;
        virtual VertexData<Vector3> F_Volume(std::vector<double> Constants) const;
        virtual VertexData<Vector3> F_Laplace(std::vector<double> Constants) const;
        virtual VertexData<Vector3> F_Edge_reg(std::vector<double> Constants) const;
        virtual VertexData<Vector3> F_Edge_reg_2(std::vector<double> Constants) const;
        virtual VertexData<Vector3> F_Face_reg(std::vector<double> Constants) const;

        virtual SparseMatrix<double> H_SurfaceTension(std::vector<double> Constants);
        virtual SparseMatrix<double> H_Bending(std::vector<double> Constants);
        virtual SparseMatrix<double> H_Bending_2(std::vector<double> Constants);
        virtual SparseMatrix<double> H_Volume(std::vector<double> Constants);
        virtual SparseMatrix<double> H_Laplace(std::vector<double> Constants);
        virtual SparseMatrix<double> H_Edge_reg(std::vector<double> Constants);
        virtual SparseMatrix<double> H_Edge_reg_2(std::vector<double> Constants);
        virtual SparseMatrix<double> H_Face_reg(std::vector<double> Constants);
        

        virtual void Calculate_energies(double* E);
        virtual void Calculate_Lag_norm(double* Norm);
        // THis function saves the value of the gradient to current gradient but before saves current gradient to previous gradient
        virtual void Calculate_gradient();
        virtual void Calculate_Jacobian();
        
        virtual SparseMatrix<double> Calculate_Hessian(); 
        virtual SparseMatrix<double> Calculate_Hessian_E();
        virtual SparseMatrix<double> Calculate_Hessian_Constraints();


        virtual void Do_nothing();

    
};