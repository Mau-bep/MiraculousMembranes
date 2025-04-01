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
        std::vector<Bead *> Beads;
        
        std::vector<double> Energy_values;
        
        VertexData<Vector3> Previous_grad;
        VertexData<Vector3> Current_grad;

        std::vector<double> Gradient_norms;

        // Constructors
        E_Handler(){};
        E_Handler(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);

        void Add_Bead(Bead *bead);

        void Add_Energy(std::string Energy_name, std::vector<double> Constants);

        
        SparseMatrix<double> H1_operator(bool CM, bool Vol_const, bool Area_const);
        SparseMatrix<double> H2_operator(bool CM, bool Vol_const, bool Area_const);
        
        virtual double E_Volume_constraint(double P0,double V ,double V_bar) const;
        virtual double E_Area_constraint(double KA,double A, double A_bar) const;
        virtual double E_SurfaceTension(double KA,double A, double A_bar) const;
        virtual double E_Bending(double H0, double KB) const;
        
        virtual VertexData<Vector3> F_Volume_constraint() const;
        virtual VertexData<Vector3> F_SurfaceTension() const;
        virtual VertexData<Vector3> F_Area_constraint() const;
        virtual VertexData<Vector3> F_Bending(double H0) const;


        void Calculate_energies(double *E);
        // THis function saves the value of the gradient to current gradient but before saves current gradient to previous gradient
        void Calculate_gradient();


    
};