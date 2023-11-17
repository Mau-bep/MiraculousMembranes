

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
    // interaction = "LJ";
    interaction = "Shifted-LJ";
    Total_force={0,0,0};
}   


    VertexData<Vector3> Bead::Gradient(){
    
    
    
    VertexData<Vector3> Force(*mesh,{0,0,0});
    VertexData<double> E_v(*mesh,0.0);


    // The problem is that FaceData is not working lmao

    Vector3 unit_r={0.0,0.0,0.0};
    Vector3 F={0,0,0};
    Vector3 H_N={0,0,0};
    Vector<Vector3> Normals(mesh->nFaces());
    Total_force={0,0,0};
    double r=0;
    // double E_v=0;
    double dual_area=0;
    Halfedge he_grad;
    Vector3 F2={0,0,0};

    double rc=2;
    double rc2=rc*rc;
    double r2;
    double alpha;
    for(Face f : mesh->faces()){
        Normals[f.getIndex()]=geometry->faceNormal(f);
    }
    if(interaction=="LJ"){
    for(Vertex v : mesh->vertices()){
        unit_r= this->Pos- geometry->inputVertexPositions[v];
        r= unit_r.norm();
        E_v[v]=4*strength*(pow(sigma/r,12)-pow(sigma/r,6));
    }
    }
    if(interaction=="Shifted-LJ"){
        for(Vertex v : mesh->vertices()){
        r2=(this->Pos-geometry->inputVertexPositions[v]).norm2();
        if(r2>rc2){
            continue;
        }
        alpha=2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3.0 );
        dual_area=geometry->barycentricDualArea(v);
        // dual_area=geometry->circumcentricDualArea(v);
        
        E_v[v]= strength*alpha*( (sigma*sigma/r2)-1  )*pow( (rc2/r2)-1 ,2.0);
        }

    }



    if(interaction=="Spring"){
        
        
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
    
    // double dual_area=0;
    if(interaction=="Shifted-LJ"){
        Vector3 u;
        rc2=rc*rc;
        
        for(Vertex v : mesh->vertices()){
            F={0,0,0};
            unit_r= this->Pos- geometry->inputVertexPositions[v];
            r= unit_r.norm();
            if(r>rc){
                continue;
            }
            unit_r=unit_r.unit();
            
            alpha=2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3.0 );
            dual_area=geometry->barycentricDualArea(v);
            
            
            F=dual_area*strength*alpha*(-2)*(  (sigma*sigma/(r*r*r))*pow(( (rc2/(r*r))-1),2.0 )+ 2*((sigma*sigma/(r*r))-1)*( (rc2/(r*r))-1 )*(rc2/(r*r*r)) )*unit_r;
            // F=strength*alpha*(-2)*(  (sigma*sigma/(r*r*r))*pow(( (rc2/(r*r))-1),2.0 )+ 2*((sigma*sigma/(r*r))-1)*( (rc2/(r*r))-1 )*(rc2/(r*r*r)) )*unit_r;
            
            F2={0,0,0};
            
            for(Face f : v.adjacentFaces()){
                he_grad=f.halfedge();
                while(he_grad.vertex().getIndex()!=v.getIndex()){
                    he_grad=he_grad.next();
                }
                
                u=geometry->inputVertexPositions[he_grad.next().vertex()]-geometry->inputVertexPositions[he_grad.next().next().vertex()];
                Vector3 Grad_vec=(0.5)* cross(Normals[f.getIndex()],u);
                
                 F2+= (1.0/3.0)*(E_v[v]+E_v[he_grad.next().vertex()]+E_v[he_grad.next().next().vertex()]  )*Grad_vec;
            }

            // F=(4*strength*( -12*pow(sigma/r,12)+6*pow(sigma/r,6))/r)*unit_r;
            // Force[v]=dot(F,Normal_v)*Normal_v;
            Force[v]=F+F2;
            Total_force-=F;
           // 0.817069 0.667744 1.22363 

        }    

    // std::cout<<"MODIFIED LJ";

    return Force;
    }
    

    if(interaction=="LJ"){
        // Vector3 grad_A={0,0,0};
        
        Vector3 u;
        for(Vertex v : mesh->vertices()){
            F={0,0,0};
            
            unit_r= this->Pos- geometry->inputVertexPositions[v];
            r= unit_r.norm();
            // Normal_v=geometry->vertexNormalAreaWeighted(v);
            unit_r=unit_r.unit();
            // dual_area=geometry->circumcentricDualArea(v);
            dual_area=geometry->barycentricDualArea(v);

            // E_v=4*strength*(pow(sigma/r,12)-pow(sigma/r,6));
            // E_v=4;

            F=(4*dual_area*strength*( -12*pow(sigma/r,12)+6*pow(sigma/r,6))/r)*unit_r;//-(E_v)*H_N;

            // F=(4*strength*( -12*pow(sigma/r,12)+6*pow(sigma/r,6))/r)*unit_r;

            // I need to work this thing, its not trivial
            F2={0,0,0};
            for(Face f : v.adjacentFaces()){
                he_grad=f.halfedge();
                while(he_grad.vertex().getIndex()!=v.getIndex() ){
                    he_grad=he_grad.next();
                }
                
                u=geometry->inputVertexPositions[he_grad.next().vertex()]-geometry->inputVertexPositions[he_grad.next().next().vertex()];
                // This E_v is the one that is not so simple.
                // I need to 
                Vector3 Grad_vec=(0.5)* cross(Normals[f.getIndex()],u);
                
                F2+= (1.0/3.0)*(E_v[v]+E_v[he_grad.next().vertex()]+E_v[he_grad.next().next().vertex()]  )*Grad_vec;
            
            }


            
            // F2=E_v*(1.0/3) *-2*geometry->vertexNormalMeanCurvature(v);
            
            
            
            Force[v]=F+F2;
            Total_force-=F;

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

    double dual_area;
    if(interaction=="Spring"){
        
        
        for(Vertex v : mesh->vertices()){
            
            r=(this->Pos-geometry->inputVertexPositions[v]).norm();

            // dual_area=geometry->barycentricDualArea(v);
            
            Total_E+=0.5*strength*(r-sigma)*(r-sigma);
        }
        return Total_E;


    }


    if(interaction=="Shifted-LJ"){
    // double dual_area;
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
        // dual_area=geometry->circumcentricDualArea(v);

        
        Total_E += dual_area*strength*alpha*( (sigma*sigma/r2)-1  )*pow( (rc2/r2)-1 ,2.0);

        // Total_E += strength*alpha*( (sigma*sigma/r2)-1  )*pow( (rc2/r2)-1 ,2.0);

    }
    // std::cout<<"THe total energy is "<<Total_E <<" \n";
    return Total_E;
        
    }


    if(interaction=="LJ"){
    for(Vertex v : mesh->vertices()){
        r=(this->Pos-geometry->inputVertexPositions[v]).norm();
        // dual_area=geometry->circumcentricDualArea(v);
        dual_area=geometry->barycentricDualArea(v);

        // Total_E+= 4*strength*(pow(sigma/r,12)-pow(sigma/r,6));
        // Total_E+=2*dual_area;
    

        // if(v.getIndex()==3){
            Total_E+= 4*dual_area*strength*(pow(sigma/r,12)-pow(sigma/r,6));
        
        //     std::cout<<" The radial distance is "<<r <<"the dual area is "<< dual_area <<" and the total computed energy is "<<Total_E<<" \n";

        // }

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
