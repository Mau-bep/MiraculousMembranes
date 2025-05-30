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
        bool boundary;
        std::vector<Bead *> Beads;
        
        std::vector<double> Energy_values;
        
        VertexData<Vector3> Previous_grad;
        VertexData<Vector3> Current_grad;

        std::vector<double> Gradient_norms;

        // Constructors
        E_Handler(){};
        E_Handler(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);
        E_Handler(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, std::vector<std::string> inputEnergyNames, std::vector<std::vector<double>> inputEnergyConstants);
        
        void Add_Bead(Bead* bead);

        void Add_Energy(std::string Energy_name, std::vector<double> Constants);

        
        SparseMatrix<double> H1_operator(bool CM, bool Vol_const, bool Area_const);
        SparseMatrix<double> H2_operator(bool CM, bool Vol_const, bool Area_const);
        
        virtual double E_Volume_constraint(std::vector<double> Constants) const;
        virtual double E_Area_constraint(std::vector<double> Constants) const;
        virtual double E_SurfaceTension(std::vector<double> Constants) const;
        virtual double E_Bending(std::vector<double> Constants) const;
        
        virtual VertexData<Vector3> F_Volume_constraint(std::vector<double> Constants) const;
        virtual VertexData<Vector3> F_SurfaceTension(std::vector<double> Constants) const;
        virtual VertexData<Vector3> F_Area_constraint(std::vector<double> Constants) const;
        virtual VertexData<Vector3> F_Bending(std::vector<double> Constants) const;
        
        virtual VertexData<Vector3> F_SurfaceTension_2(std::vector<double> Constants) const;
        virtual VertexData<Vector3> F_Bending_2(std::vector<double> Constants) const;
        virtual VertexData<Vector3> F_Volume_constraint_2(std::vector<double> Constants) const;

        virtual SparseMatrix<double> H_SurfaceTension(std::vector<double> Constants);
        virtual SparseMatrix<double> H_Bending(std::vector<double> Constants);


        virtual void Calculate_energies(double* E);
        // THis function saves the value of the gradient to current gradient but before saves current gradient to previous gradient
        virtual void Calculate_gradient();
        virtual void Do_nothing();

    
};