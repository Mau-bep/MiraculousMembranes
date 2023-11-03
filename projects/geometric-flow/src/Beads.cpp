

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
    interaction = "Shifted-LJ";
    // interaction = "Shifted-LJ";
    Total_force={0,0,0};
}


VertexData<Vector3> Bead::Gradient(){
    VertexData<Vector3> Force(*mesh,{0,0,0});

    Vector3 unit_r;
    Vector3 F;
    Total_force={0,0,0};
    double r;
    

    if(interaction=="Spring"){
        double dual_area;
        
        for(Vertex v : mesh->vertices()){
            
            unit_r= this->Pos- geometry->inputVertexPositions[v];
            r= unit_r.norm();
            unit_r=unit_r.unit();
            // dual_area=geometry->barycentricDualArea(v);

            F=strength*( r-sigma)*unit_r;

            Force[v]=F;
            Total_force-=Force[v];

        }
        return Force;


    }


    // Vector3 Normal_v;
    double rc=2;
    double dual_area=0;
    if(interaction=="Shifted-LJ"){
        double rc2=rc*rc;
        double alpha;
        for(Vertex v : mesh->vertices()){
            unit_r= this->Pos- geometry->inputVertexPositions[v];
            r= unit_r.norm();
            if(r>rc){
                continue;
            }
            // std::cout<<"Are we here at some point?\n";

            // Normal_v=geometry->vertexNormalAreaWeighted(v);
            unit_r=unit_r.unit();
            alpha=2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3.0 );
            dual_area=geometry->barycentricDualArea(v);
            F=strength*alpha*(-2)*(  (sigma*sigma/(r*r*r))*pow(( (rc2/(r*r))-1),2.0 )+ 2*((sigma*sigma/(r*r))-1)*( (rc2/(r*r))-1 )*(rc2/(r*r*r)) )*unit_r;


            // F=(4*strength*( -12*pow(sigma/r,12)+6*pow(sigma/r,6))/r)*unit_r;
            // Force[v]=dot(F,Normal_v)*Normal_v;
            Force[v]=F;
            Total_force-=Force[v];

        }    


    return Force;
    }
    
    if(interaction=="LJ"){
        for(Vertex v : mesh->vertices()){
            unit_r= this->Pos- geometry->inputVertexPositions[v];
            r= unit_r.norm();
            // Normal_v=geometry->vertexNormalAreaWeighted(v);
            unit_r=unit_r.unit();

            F=(4*strength*( -12*pow(sigma/r,12)+6*pow(sigma/r,6))/r)*unit_r;
            // Force[v]=dot(F,Normal_v)*Normal_v;
            Force[v]=F;
            Total_force-=Force[v];

        }

        return Force;
    }


    // std::cout<< "THe surface tension force in magnitude is: "<< -1*lambda*sqrt(Force.transpose()*Force) <<"\n";

    return Force;
}


void Bead::Set_Force(Vector3 Force){
    
    this->Total_force=Force;
    return ;
}

double Bead::Energy() {
    double Total_E=0.0;
    double r;


    if(interaction=="Spring"){
        double dual_area;
        
        for(Vertex v : mesh->vertices()){
            
            r=(this->Pos-geometry->inputVertexPositions[v]).norm();

            // dual_area=geometry->barycentricDualArea(v);
            Total_E+=0.5*strength*(r-sigma)*(r-sigma);
        }
        return Total_E;


    }




    if(interaction=="Shifted-LJ"){
    double dual_area;
    double alpha;
    double rc2=4.0;
    double r2;
    for(Vertex v : mesh->vertices()){
        r2=(this->Pos-geometry->inputVertexPositions[v]).norm2();
        if(r2>rc2){
            continue;
        }
        alpha=2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3.0 );
        dual_area=geometry->barycentricDualArea(v);

        
        Total_E += strength*alpha*( (sigma*sigma/r2)-1  )*pow( (rc2/r2)-1 ,2.0);

        // Total_E+= 4*strength*(pow(sigma/r,12)-pow(sigma/r,6));

    }
    // std::cout<<"THe total energy is "<<Total_E <<" \n";
    return Total_E;
        
    }


    if(interaction=="LJ"){
    for(Vertex v : mesh->vertices()){
        r=(this->Pos-geometry->inputVertexPositions[v]).norm();
        Total_E+= 4*strength*(pow(sigma/r,12)-pow(sigma/r,6));

    }
    return Total_E;
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

    // Total_force={0,0,0};

    return;
}
