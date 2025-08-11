#pragma once 

#include <Eigen/Core>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include <fstream>
#include <omp.h>

using namespace geometrycentral;
using namespace geometrycentral::surface;

class Interaction {

    public:


        // What are the parameters for the interaction

        std::vector<double> Energy_constants;
        std::vector<double> Bead_params;

        std::vector<std::string> E_Features;
        std::vector<int> E_Features_val;


        virtual double Tot_Energy() = 0;

        virtual Vector3 Gradient() = 0;


        virtual double E_r(double r, std::vector<double> Energy_constants)  = 0;
        virtual Vector3 F_r(double r,Vector3 r_vec, std::vector<double> Energy_constants) = 0;
        



};


class Integrated_Interaction: public Interaction {
    public:
        Integrated_Interaction() {};

        virtual double Tot_Energy();
        virtual Vector3 Gradient();

        virtual double E_r(double r, std::vector<double> Energy_constants) = 0;
        virtual Vector3 F_r(double r, Vector3 r_vec, std::vector<double> Energy_constants) = 0;

};

class Normal_dot_Interaction: public Interaction {

    public:

        Normal_dot_Interaction() {};

        virtual double Tot_Energy();
        virtual Vector3 Gradient();

        virtual double E_r(double r, std::vector<double> Energy_constants) = 0;
        virtual Vector3 F_r(double r, Vector3 r_vec, std::vector<double> Energy_constants) = 0;


};



class LJ_Normal: public Normal_dot_Interaction {

    public:

        LJ_Normal(std::vector<double> params, std::vector<double> Bead_params) {
            Energy_constants = params;
            Bead_params = Bead_params;
        }

        double E_r(double r, std::vector<double> Energy_constants) override {
            // r is the distance between the two beads
            // Energy_constants[0] is the strength of the interaction
            // Energy_constants[1] is the sigma of the interaction
            double sigma = Energy_constants[1];
            double epsilon = Energy_constants[0];
            return 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6));
        }
        Vector3 F_r(double r, Vector3 r_vec, std::vector<double> Energy_constants) override {
            // r_vec is the vector between the two beads
            // Energy_constants[0] is the strength of the interaction
            // Energy_constants[1] is the sigma of the interaction
            double sigma = Energy_constants[1];
            double epsilon = Energy_constants[0];
            return -24 * epsilon * (2 * pow(sigma / r, 12) - pow(sigma / r, 6)) * (r_vec / r);
        }


};

class LJ : public Integrated_Interaction{
    public: 
        LJ(std::vector<double> params, std::vector<double> Bead_params) {
            Energy_constants = params;
            Bead_params = Bead_params;
        }
        double E_r(double r, std::vector<double> Energy_constants) override {
            // r is the distance between the two beads
            // Energy_constants[0] is the strength of the interaction
            // Energy_constants[1] is the sigma of the interaction
            double sigma = Energy_constants[1];
            double epsilon = Energy_constants[0];
            return 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6));
        }
        Vector3 F_r(double r, Vector3 r_vec, std::vector<double> Energy_constants) override {
            // r_vec is the vector between the two beads
            // Energy_constants[0] is the strength of the interaction
            // Energy_constants[1] is the sigma of the interaction
            double sigma = Energy_constants[1];
            double epsilon = Energy_constants[0];
            return -24 * epsilon * (2 * pow(sigma / r, 12) - pow(sigma / r, 6)) * (r_vec / r);
        }
};


