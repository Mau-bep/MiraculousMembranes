#pragma once

#include <Eigen/Core>
#include <omp.h>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "Beads.h"

class Bead;

using namespace geometrycentral;
using namespace geometrycentral::surface;

class Interaction {

    public:
        ManifoldSurfaceMesh* mesh;
        VertexPositionGeometry* geometry;

        // The interaction needs to know which bead it is
        Bead* Bead_1;
        // What are the parameters for the interaction

        std::vector<double> Energy_constants;


        std::vector<std::string> E_Features;
        std::vector<int> E_Features_val;


        virtual double Bond_energy();
        virtual Vector3 Bond_force();
        virtual double Tot_Energy() = 0;
        virtual VertexData<Vector3>  Gradient() = 0;
        virtual SparseMatrix<double> Hessian() = 0;

        virtual std::vector<Eigen::Triplet<double>> Hessian_bonds_triplet();

        virtual double E_r(double r, std::vector<double> Energy_constants)  = 0;
        virtual double dE_r(double r, std::vector<double> Energy_constants) = 0;
        virtual double ddE_r(double r, std::vector<double> Energy_constants) = 0;
        // virtual Vector3 F_r(double r,Vector3 r_vec, std::vector<double> Energy_constants) = 0;
        




};

class No_mem_Inter: public Interaction{
    public:
    
            No_mem_Inter() {};
    
            double Tot_Energy() override {
                return Bond_energy();
            }
            virtual VertexData<Vector3> Gradient();
            virtual SparseMatrix<double> Hessian();

            virtual double E_r(double r, std::vector<double> Energy_constants) override {
                return 0.0;
            }
            virtual double dE_r(double r, std::vector<double> Energy_constants) override {
                return 0.0;
            }
            virtual double ddE_r(double r, std::vector<double> Energy_constants) override {
                return 0.0;
            }

};
class Integrated_Interaction: public Interaction {
    public:

        Integrated_Interaction() {};

        virtual double Tot_Energy();
        virtual VertexData<Vector3> Gradient();
        virtual SparseMatrix<double> Hessian();

        virtual double E_r(double r, std::vector<double> Energy_constants) = 0;
        virtual double dE_r(double r, std::vector<double> Energy_constants) = 0;
        virtual double ddE_r(double r, std::vector<double> Energy_constants) = 0;
        // virtual Vector3 F_r(double r, Vector3 r_vec, std::vector<double> Energy_constants) = 0;

};

class Normal_dot_Interaction: public Interaction {

    public:

        Normal_dot_Interaction() {};

        virtual double Tot_Energy();
        virtual VertexData<Vector3>  Gradient();
        virtual SparseMatrix<double> Hessian();
        virtual double E_r(double r, std::vector<double> Energy_constants) = 0;
        virtual double dE_r(double r, std::vector<double> Energy_constants) = 0;
        virtual double ddE_r(double r, std::vector<double> Energy_constants) = 0;
        // virtual Vector3 F_r(double r, Vector3 r_vec, std::vector<double> Energy_constants) = 0;


};


class Frenkel_Normal: public Normal_dot_Interaction{

    public:

    Frenkel_Normal(){}

    Frenkel_Normal(std::vector<double> params) {
        Energy_constants = params;
    }
    Frenkel_Normal(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, std::vector<double> params) {
        Energy_constants = params;
        mesh = inputMesh;
        geometry = inputGeo;
    }
    

    double E_r(double r, std::vector<double> Energy_constants) override {
        // r is the distance between the two beads
        // Energy_constants[0] is the strength of the interaction
        // Energy_constants[1] is the sigma of the interaction
        double sigma = Energy_constants[1];
        double epsilon = Energy_constants[0];
        double rc = Energy_constants[2]; // cutoff distance
        if( r >= rc) {
            return 0.0; // No interaction beyond cutoff
        }
        double rc2 = rc * rc;
        double r2 = r*r;
        double alpha = 2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3.0 );

        return epsilon*alpha*( (sigma*sigma/r2)-1  )*pow( (rc2/r2)-1 ,2.0);
        
    }
    double dE_r(double r, std::vector<double> Energy_constants) override {


        double epsilon = Energy_constants[0];
        double sigma = Energy_constants[1];
        double rc = Energy_constants[2]; // cutoff distance
        if (r >= rc) {
            return 0.0; // No interaction beyond cutoff
        }
        double rc2 = rc * rc;
        double r2 = r*r;
        double alpha = 2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3.0 );

        // Ok so now i need to calculate the derivative
        double Q1 = -2*alpha/(r2*r);
        double Q2 = rc2/r2 -1;
        double Q3 = sigma*sigma*( rc2/r2 -1) + 2*rc2 *( sigma*sigma/r2 -1);

        return epsilon*Q1*Q2*Q3;

    }
    double ddE_r(double r, std::vector<double> Energy_constants) override {

        double sigma = Energy_constants[1];
        double epsilon = Energy_constants[0];
        double rc = Energy_constants[2]; // cutoff distance
        if (r >= rc) {
            return 0.0; // No interaction beyond cutoff
        }
        double rc2 = rc * rc;
        double r2 = r*r;
        double alpha = 2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3.0 );

        // Ok so now i need to calculate the second derivative
        double Q1 = -2*alpha/(r2*r);
        double Q2 = rc2/r2 -1;
        double Q3 = sigma*sigma*( rc2/r2 -1) + 2*rc2 *( sigma*sigma/r2 -1);

        double dQ1 = 6*alpha/(r2*r2);
        double dQ2 = -2 *rc2/(r2*r);
        double dQ3 = sigma*sigma*(-2)*rc2/(r2*r) -  4*rc2*sigma*sigma/(r2*r);

        return dQ1*Q2*Q3 + Q1*dQ2*Q3 + Q1*Q2*dQ3;


    }


};

class Constant_Normal: public Normal_dot_Interaction{
    public:

    Constant_Normal(std::vector<double> params) {
        Energy_constants = params;
    }
    double E_r(double r, std::vector<double> Energy_constants) override {
     
    
     return 1.0;
    }
    double dE_r(double r, std::vector<double> Energy_constants) override {
       
        return 0.0; // Constant interaction, no gradient
    }
    double ddE_r(double r, std::vector<double> Energy_constants) override {
        

        return 0.0; // Constant interaction, no second derivative
    }

};
class One_over_r_Normal : public Normal_dot_Interaction{
    public:

        One_over_r_Normal(){}

        One_over_r_Normal(std::vector<double> params){
            Energy_constants = params;
        }

        One_over_r_Normal(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, std::vector<double> params) {
        Energy_constants = params;
        mesh = inputMesh;
        geometry = inputGeo;
        }

        double E_r(double r, std::vector<double> Energy_constants) override{
            //  1/r = 
            // Which are the energy constants
            double epsilon = Energy_constants[0];
            
            return epsilon/r;

        }

        double dE_r(double r, std::vector<double> Energy_constants) override{
            double epsilon = Energy_constants[0];

            return -epsilon/(r*r);

        }

        double ddE_r(double r, std::vector<double> Energy_constants) override{
            double epsilon = Energy_constants[0];
            
            return 2*epsilon/(r*r*r);
        }
        



};

class One_over_r : public Integrated_Interaction{
    public:

        One_over_r(){}

        One_over_r(std::vector<double> params){
            Energy_constants = params;
        }

        One_over_r(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, std::vector<double> params) {
        Energy_constants = params;
        mesh = inputMesh;
        geometry = inputGeo;
        }

        double E_r(double r, std::vector<double> Energy_constants) override{
            //  1/r = 
            // Which are the energy constants
            double epsilon = Energy_constants[0];
            
            return epsilon/r;

        }

        double dE_r(double r, std::vector<double> Energy_constants) override{
            double epsilon = Energy_constants[0];

            return -epsilon/(r*r);

        }

        double ddE_r(double r, std::vector<double> Energy_constants) override{
            double epsilon = Energy_constants[0];
            
            return 2*epsilon/(r*r*r);
        }
        



};

class LJ_Normal: public Normal_dot_Interaction {

    public:

        LJ_Normal(){}

        LJ_Normal(std::vector<double> params) {
            Energy_constants = params;
        }
        LJ_Normal(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, std::vector<double> params) {
            Energy_constants = params;
            mesh = inputMesh;
            geometry = inputGeo;
        }

        

        double E_r(double r, std::vector<double> Energy_constants) override {
            // r is the distance between the two beads
            // Energy_constants[0] is the strength of the interaction
            // Energy_constants[1] is the sigma of the interaction
            double sigma = Energy_constants[1];
            double epsilon = Energy_constants[0];
            double rc = Energy_constants[2]; // cutoff distance
            double shift = Energy_constants[3]; // shift value
            if( r >= rc) {
                return 0.0; // No interaction beyond cutoff
            }
            return 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6)) + shift;
        }
        double dE_r(double r, std::vector<double> Energy_constants) override {

            // r is the distance between the two beads
            // Energy_constants[0] is the strength of the interaction
            // Energy_constants[1] is the sigma of the interaction

            double sigma = Energy_constants[1];
            double epsilon = Energy_constants[0];
            double rc = Energy_constants[2]; // cutoff distance
            if (r >= rc) {
                return 0.0; // No interaction beyond cutoff
            }
            return -24 * epsilon * (2 * pow(sigma / r, 12) - pow(sigma / r, 6)) / r;
        }
        double ddE_r(double r, std::vector<double> Energy_constants) override {
            // r is the distance between the two beads
            // Energy_constants[0] is the strength of the interaction
            // Energy_constants[1] is the sigma of the interaction
            double sigma = Energy_constants[1];
            double epsilon = Energy_constants[0];
            double rc = Energy_constants[2]; // cutoff distance
            if (r >= rc) {
                return 0.0; // No interaction beyond cutoff
            }
            return 24 * epsilon * (2 * pow(sigma / r, 12) - pow(sigma / r, 6)) / (r * r) + 24 * epsilon * (24 * pow(sigma / r, 12) - 6 * pow(sigma / r, 6)) / ( r*r);
        }

        


};

class LJ : public Integrated_Interaction{
    public: 
        LJ(){}

        LJ(std::vector<double> params) {
            Energy_constants = params;
        }
        LJ(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, std::vector<double> params) {
            Energy_constants = params;
            mesh = inputMesh;
            geometry = inputGeo;
        }

        double E_r(double r, std::vector<double> Energy_constants) override {
            // r is the distance between the two beads
            double epsilon = Energy_constants[0];
            double sigma = Energy_constants[1];
            double rc = Energy_constants[2]; // cutoff distance
            double shift = Energy_constants[3]; // shift value
            if (r >= rc && rc > 0.0) {
                return 0.0; // No interaction beyond cutoff
            }
            return 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6))+shift;
        }
        double dE_r(double r, std::vector<double> Energy_constants) override {

            // r is the distance between the two beads
            // Energy_constants[0] is the strength of the interaction
            // Energy_constants[1] is the sigma of the interaction
            double sigma = Energy_constants[1];
            double epsilon = Energy_constants[0];
            double rc = Energy_constants[2]; // cutoff distance
            if (r >= rc && rc > 0) {
                return 0.0; // No interaction beyond cutoff
            }
            return -24 * epsilon * (2 * pow(sigma / r, 12) - pow(sigma / r, 6)) / r;
        }
        double ddE_r(double r, std::vector<double> Energy_constants) override {
            // r is the distance between the two beads
            // Energy_constants[0] is the strength of the interaction
            // Energy_constants[1] is the sigma of the interaction
            double sigma = Energy_constants[1];
            double epsilon = Energy_constants[0];
            double rc = Energy_constants[2]; // cutoff distance
            if (r >= rc && rc > 0) {
                return 0.0; // No interaction beyond cutoff
            }
            return 24 * epsilon * (2 * pow(sigma / r, 12) - pow(sigma / r, 6)) / (r * r) + 24 * epsilon * (24 * pow(sigma / r, 12) - 6 * pow(sigma / r, 6)) / ( r*r);
        }
        
};


