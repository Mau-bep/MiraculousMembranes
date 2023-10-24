

#include "Beads.h"
#include <fstream>
#include <omp.h>
#include <chrono>
#include <Eigen/Core>


using namespace std;
/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
Bead::Bead(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo,Vector3 Position,double input_sigma , double strg) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
    Pos = Position;
    sigma = input_sigma;
    strength = strg;
    interaction = "LJ";
    Total_force={0,0,0};
}


VertexData<Vector3> Bead::Gradient(){
    VertexData<Vector3> Force(*mesh,{0,0,0});

    Vector3 unit_r;
    Vector3 F;
    Total_force={0,0,0};
    double r;
    if(interaction=="LJ"){
        for(Vertex v : mesh->vertices()){
            unit_r= this->Pos- geometry->inputVertexPositions[v];
            r= unit_r.norm();
            
            unit_r=unit_r.unit();
            F=(4*strength*( -12*pow(sigma/r,12)+6*pow(sigma/r,6))/r)*unit_r;
            Force[v]=F;
            Total_force-=F;


        }

        return Force;
    }


    // std::cout<< "THe surface tension force in magnitude is: "<< -1*lambda*sqrt(Force.transpose()*Force) <<"\n";

    return Force;
}

double Bead::Energy() {
    double Total_E=0.0;
    double r;
    for(Vertex v : mesh->vertices()){
        r=(this->Pos-geometry->inputVertexPositions[v]).norm();

        Total_E+= 4*strength*(pow(sigma/r,12)-pow(sigma/r,6));

    }


    return Total_E;
}

void Bead::Reset_bead(Vector3 Actual_pos){
    this->Pos=Actual_pos;
    return;
}

void Bead::Move_bead(double dt,Vector3 center) {
    // We can add more things here later
    this->Pos=this->Pos+Total_force*dt -center;

    Total_force={0,0,0};

    return;
}
