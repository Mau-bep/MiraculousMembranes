#pragma once 

#include <Eigen/Core>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include <fstream>
#include <omp.h>


class Interaction {

    public:


        // What are the parameters for the interaction

        std::vector<double> Energy_constants;
        std::vector<double> Bead_params;

        virtual double Tot_Energy();
        virtual Vector3 Tot_Force();


        virtual double E_r(double r, std::vector<double> Energy_constants);
        virtual Vector3 F_r(double r,Vector3 r_vec, std::vector<double> Energy_constants);
        



}