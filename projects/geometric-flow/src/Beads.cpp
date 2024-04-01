

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
    // interaction = "test_Full";
    interaction = "Shifted_LJ_Normal_var";
   
    // interaction = "test_angle_normal_r_normalized";
    // interaction = "test_angle_normal_r_normalized_LJ_Full";
    Total_force={0,0,0};
}   


    VertexData<Vector3> Bead::Gradient(){
    
    
    
    VertexData<Vector3> Force(*mesh,{0,0,0});
    VertexData<double> E_v(*mesh,0.0);


    // The problem is that FaceData is not working lmao

    Vector3 unit_r;
    Vector3 F;
    Vector3 H_N;
    Vector<Vector3> Normals(mesh->nFaces());
    Vector<Vector3> Angle_Normals(mesh->nVertices());
    Vector<double> Angle_Normals_norm(mesh->nVertices());
    // I need to add the unit_r here
    Vector<Vector3> Unit_rs(mesh->nVertices());
    VertexData<double> rs(*mesh,0.0);
    Vector<bool> Consider_vertex(mesh->nVertices());
    Vector<double> Dual_areas(mesh->nVertices());

    Total_force={0,0,0};
    double r;
    // double E_v=0;
    double dual_area;
    Halfedge he_grad;
    Vector3 F2;
    Vector3 F3;
    Vector3 F4;
    Vector3 F5;
    Vector3 F5T;
    Vector3 F2T;
    Vector3 F4T;
    Vector3 F3T;
    double rc=1.2;
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
    if(interaction=="Shifted-LJ" ||  interaction == "Shifted_LJ_Normal"||interaction == "Shifted_LJ_Normal_var" || interaction=="test_angle_normal_r_normalized"||interaction=="test_angle_normal_r_normalized_LJ"|| interaction=="test_angle_normal_r_normalized_LJ_Full"){
        // std::cout<<"Loading essential quantities\n";
        for(Vertex v : mesh->vertices()){

        Angle_Normals[v.getIndex()]=geometry->vertexNormalAngleWeighted(v);
        Angle_Normals_norm[v.getIndex()]=Angle_Normals[v.getIndex()].norm();
        unit_r=(this->Pos-geometry->inputVertexPositions[v]);
        rs[v]=unit_r.norm();
        Unit_rs[v.getIndex()]=unit_r/rs[v];
        
        r2=rs[v]*rs[v];
        Consider_vertex[v.getIndex()]=dot(Angle_Normals[v.getIndex()],Unit_rs[v.getIndex()])>0;
        
        Dual_areas[v.getIndex()]=geometry->barycentricDualArea(v);
        // Dual_areas[v.getIndex()]=geometry->circumcentricDualArea(v);

        
        if(rs[v]>rc){
            E_v[v]=0.0;
           
        }
        else{
        alpha=2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3.0 );
        dual_area=Dual_areas[v.getIndex()];
        // dual_area=geometry->circumcentricDualArea(v);
        
        E_v[v]= strength*alpha*( (sigma*sigma/r2)-1  )*pow( (rc2/r2)-1 ,2.0);
        // Angle_Normals[v.getIndex()]=geometry->vertexNormalAngleWeighted(v);
        }
        
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
            if(r<rc){
                
            
            unit_r=unit_r.unit();
            
            alpha=2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3 );
            dual_area=Dual_areas[v.getIndex()];
            
            
            F=dual_area*strength*alpha*(-2)*(  (sigma*sigma/(r*r*r))*pow(( (rc2/(r*r))-1),2 )+ 2*((sigma*sigma/(r*r))-1)*( (rc2/(r*r))-1 )*(rc2/(r*r*r)) )*unit_r;
            }

            F2={0,0,0};
            
            for(Face f : v.adjacentFaces()){
                // i want to check that any of the vertices in this face contain information
                
                he_grad=f.halfedge();
                if( rs[he_grad.vertex()]> rc && rs[he_grad.next().vertex()]>rc && rs[he_grad.next().next().vertex()]>rc){
                    continue;
                }

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




    if(interaction=="test_angle"){

        Vertex v_1 =mesh->vertex(1);
        Halfedge he = v_1.halfedge();
        
        // If the energy is equal to the angle in this corner, the gradient only depends on three vertices
        Vertex v_2 =v_1.halfedge().next().vertex();
        Vertex v_3 = v_1.halfedge().next().next().vertex();
        Vector3 Normal = geometry->faceNormal(he.face());
        Force[v_1]=geometry->Angle_grad(he.corner(),v_1,Normal);
        Force[v_2]=geometry->Angle_grad(he.corner(),v_2,Normal);
        Force[v_3]=geometry->Angle_grad(he.corner(),v_3,Normal);

        return Force;


    }
        if(interaction=="test_normal"){

        double Face_area;
        Vector3 Face_normal;
        Vector3 grad;
        Halfedge he;
        // for(Vertex v: mesh->vertices()){
        for(Face f: mesh->faces()){
            Face_area=geometry->faceArea(f);
            Face_normal=geometry->faceNormal(f);
            // So the energy is the multiplication of the Normal in the middle
            // Times a constant vector, therefore i have to do the gradient of the normal 

            he = f.halfedge();

            grad=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal),sqrt(3)*Vector3({0,1,0}))*Face_normal/(2*Face_area);

            Force[he.vertex()]+=grad;
            
            he = he.next();
            grad=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal),sqrt(3)*Vector3({0,1,0}))*Face_normal/(2*Face_area);
            Force[he.vertex()]+=grad;

            he = he.next();
            grad=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal),sqrt(3)*Vector3({0,1,0}))*Face_normal/(2*Face_area);
            Force[he.vertex()]+=grad;



        }
        // }
        return Force;


    }

        if(interaction=="test_angle_normal"){

        double Face_area;
        Vector3 Face_normal;
        Vector3 grad;
        // Halfedge he;
        Face f;
        Vector3 Vertex_Normal;
        Halfedge he2;
        
        
        // FOr the record, if i only consider one dot product the first term is correct in that vertex 
        // I now want to add the contributions of the neighbors to this point.
        for(Vertex v: mesh->vertices()){
            for(Halfedge he : v.outgoingHalfedges()){
            if(!he.isInterior()){
                continue;
            }
            f =he.face();
            Face_area=geometry->faceArea(f);
            Face_normal=geometry->faceNormal(f);
            
     
            grad=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal),sqrt(3)*Vector3({0,1,0}))*Face_normal/(2*Face_area);

            Force[he.vertex()]+=grad*geometry->angle(he.corner());
            Force[he.vertex()]+=geometry->Angle_grad(he.corner(),v,Face_normal)*dot(Face_normal,sqrt(3)*Vector3({0,1,0}));

            he2 = he.next();
            // grad=-1*dot(cross(geometry->inputVertexPositions[he2.next().next().vertex()]-geometry->inputVertexPositions[he2.next().vertex()],Face_normal),sqrt(3)*Vector3({0,1,0}))*Face_normal/(2*Face_area);
            Force[he.vertex()]+=grad*geometry->angle(he2.corner());
            Force[he.vertex()]+=geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,sqrt(3)*Vector3({0,1,0}));

            he2 = he2.next();
            // grad=-1*dot(cross(geometry->inputVertexPositions[he2.next().next().vertex()]-geometry->inputVertexPositions[he2.next().vertex()],Face_normal),sqrt(3)*Vector3({0,1,0}))*Face_normal/(2*Face_area);
            Force[he.vertex()]+=grad*geometry->angle(he2.corner());
            Force[he.vertex()]+=geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,sqrt(3)*Vector3({0,1,0}));
           


        }
        }
        // }
        // }
        // }
        return Force;


    }


        if(interaction=="test_angle_normal_r"){

        double Face_area;
        Vector3 Face_normal;
        Vector3 grad;
        // Halfedge he;
        Face f;
        Vector3 Vertex_Normal;
        Halfedge he2;
        Vector3 unit_r;
        double r;

        
        // FOr the record, if i only consider one dot product the first term is correct in that vertex 
        // I now want to add the contributions of the neighbors to this point.
        for(Vertex v: mesh->vertices()){
            // Vertex v=mesh->vertex(3);
            for(Halfedge he : v.outgoingHalfedges()){
            if(!he.isInterior()){
                continue;
            }
            f =he.face();
            Face_area=geometry->faceArea(f);
            Face_normal=geometry->faceNormal(f);
            
            unit_r=(Pos-geometry->inputVertexPositions[v]);
            r=unit_r.norm();
            unit_r=unit_r/r;
            grad=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal),unit_r)*Face_normal/(2*Face_area);

            Force[he.vertex()]+=grad*geometry->angle(he.corner());
            Force[he.vertex()]+=geometry->Angle_grad(he.corner(),v,Face_normal)*dot(Face_normal,unit_r);
            
            Force[he.vertex()]+=geometry->angle(he.corner())*(Face_normal-unit_r*dot(unit_r,Face_normal))/r;





            he2 = he.next();
            
            unit_r=(Pos-geometry->inputVertexPositions[he2.vertex()]);
            r=unit_r.norm();
            unit_r=unit_r/r;
            
            grad=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal),unit_r)*Face_normal/(2*Face_area);

            Force[he.vertex()]+=grad*geometry->angle(he2.corner());
            Force[he.vertex()]+=geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,unit_r);

            he2 = he2.next();
            
            unit_r=(Pos-geometry->inputVertexPositions[he2.vertex()]);
            r=unit_r.norm();
            unit_r=unit_r/r;

            grad=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal),unit_r)*Face_normal/(2*Face_area);
            Force[he.vertex()]+=grad*geometry->angle(he2.corner());
            Force[he.vertex()]+=geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,unit_r);
           


        }
        }
        // }
        // }
        // }
        return Force;


    }


        if(interaction=="test_angle_normal_r_normalized"){

        double Face_area;
        Vector3 Face_normal;
        Vector3 grad;
        Vector3 grad2;
        // Halfedge he;
        Face f;
        Vector3 Vertex_Normal;
        Vector3 Vertex_Angle_Normal;
        double Vertex_Angle_Normal_norm;
        Halfedge he2;
        Vector3 unit_r;
        double r;

        
        // FOr the record, if i only consider one dot product the first term is correct in that vertex 
        // I now want to add the contributions of the neighbors to this point.
        for(Vertex v: mesh->vertices()){
            // Vertex v=mesh->vertex(3);
            for(Halfedge he : v.outgoingHalfedges()){
            if(!he.isInterior()){
                continue;
            }
            f =he.face();
            Face_area=geometry->faceArea(f);
            Face_normal=geometry->faceNormal(f);
            
            // First vertex 
            Vertex_Angle_Normal=geometry->vertexNormalAngleWeighted(v);
            Vertex_Angle_Normal_norm=Vertex_Angle_Normal.norm();
            

            unit_r=(Pos-geometry->inputVertexPositions[v]);
            r=unit_r.norm();
            unit_r=unit_r/r;

            if(dot(Vertex_Angle_Normal,unit_r)>0 ){

            grad=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal),unit_r)*Face_normal/(2*Face_area);

            grad2=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal ),Vertex_Angle_Normal )*Face_normal/(2*Face_area);
            
            


            Force[he.vertex()]+=grad*geometry->angle(he.corner())/Vertex_Angle_Normal_norm;

            Force[he.vertex()]+=geometry->Angle_grad(he.corner(),v,Face_normal)*dot(Face_normal,unit_r)/Vertex_Angle_Normal_norm;
            
            Force[he.vertex()]+=geometry->angle(he.corner())*(Face_normal-unit_r*dot(unit_r,Face_normal))/(r*Vertex_Angle_Normal_norm);

            // Lets add the terms associated with the norm

            Force[he.vertex()]+=-1*grad2*geometry->angle(he.corner())*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);

            Force[he.vertex()]+=-1*geometry->Angle_grad(he.corner(),v,Face_normal)*dot(Face_normal,Vertex_Angle_Normal)*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);

            }

            // Second vertex
            
            he2 = he.next();
            Vertex_Angle_Normal=geometry->vertexNormalAngleWeighted(he2.vertex());
            Vertex_Angle_Normal_norm=Vertex_Angle_Normal.norm();


            unit_r=(Pos-geometry->inputVertexPositions[he2.vertex()]);
            r=unit_r.norm();
            unit_r=unit_r/r;
            if(dot(Vertex_Angle_Normal,unit_r)>0 ){

            grad=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal),unit_r)*Face_normal/(2*Face_area);

            grad2=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal ),Vertex_Angle_Normal )*Face_normal/(2*Face_area);


            Force[he.vertex()]+=grad*geometry->angle(he2.corner())/Vertex_Angle_Normal_norm;

            Force[he.vertex()]+=geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,unit_r)/Vertex_Angle_Normal_norm;

            // Now i add the terms associated with the norm

            Force[he.vertex()]+=-1*grad2*geometry->angle(he2.corner())*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);

            Force[he.vertex()]+=-1*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,Vertex_Angle_Normal)*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);
            }


            // Third vertex
            he2 = he2.next();
            
            Vertex_Angle_Normal=geometry->vertexNormalAngleWeighted(he2.vertex());
            Vertex_Angle_Normal_norm=Vertex_Angle_Normal.norm();

            unit_r=(Pos-geometry->inputVertexPositions[he2.vertex()]);
            r=unit_r.norm();
            unit_r=unit_r/r;

            if(dot(Vertex_Angle_Normal,unit_r)>0 ){
            grad=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal),unit_r)*Face_normal/(2*Face_area);

            grad2=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal ),Vertex_Angle_Normal )*Face_normal/(2*Face_area);

            
            Force[he.vertex()]+=grad*geometry->angle(he2.corner())/Vertex_Angle_Normal_norm;
            
            Force[he.vertex()]+=geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,unit_r)/Vertex_Angle_Normal_norm;

           // Now i add the terms associated with the norm

            Force[he.vertex()]+=-1*grad2*geometry->angle(he2.corner())*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);

            Force[he.vertex()]+=-1*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,Vertex_Angle_Normal)*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);
            
            }



        }
        }
  
        return Force;


    }
       if(interaction=="Shifted_LJ_Normal"){


        int v1_idx;
        int v2_idx;
        int v3_idx;
        double Face_area;
    
        double dual_area2;
        double dual_area3;
        Vector3 Face_normal;
        Vector3 grad;
        Vector3 grad2;
        // Halfedge he;
        Face f;
        Vector3 Vertex_Normal;
        Vector3 Vertex_Angle_Normal;
        Vector3 Vertex_Angle_Normal2;
        Vector3 Vertex_Angle_Normal3;
        double Vertex_Angle_Normal_norm;
        double Vertex_Angle_Normal_norm2;
        double Vertex_Angle_Normal_norm3;
        
        Halfedge he2;
        Vector3 unit_r;
        Vector3 unit_r1;
        Vector3 unit_r2;
        Vector3 unit_r3;
        double r;
        double r2;
        double r3;

        Vector3 u;
        rc=1.2;
        rc2=rc*rc;

        // FOr the record, if i only consider one dot product the first term is correct in that vertex 
        // I now want to add the contributions of the neighbors to this point.

        
        for(Vertex v : mesh->vertices()){
            v1_idx=v.getIndex();
        // Vertex v =mesh->vertex(0);
            Force[v1_idx]={0,0,0};
            
            F={0,0,0};
            unit_r= Unit_rs[v1_idx];
            r = rs[v1_idx];
            dual_area=Dual_areas[v1_idx];
            
            
            if(r<rc){
           
            // Interaction force

            alpha=2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3 );
            
            Vertex_Angle_Normal=Angle_Normals[v1_idx];
        
            if(Consider_vertex[v1_idx]){
            F+=dual_area*dot(Vertex_Angle_Normal,unit_r)*strength*alpha*(-2)*(  (sigma*sigma/(r*r*r))*pow(( (rc2/(r*r))-1),2 )+ 2*((sigma*sigma/(r*r))-1)*( (rc2/(r*r))-1 )*(rc2/(r*r*r)) )*unit_r/Angle_Normals_norm[v1_idx];
            Total_force-=F;

            }


            }
     

            // Vertex v=mesh->vertex(3);
            for(Halfedge he : v.outgoingHalfedges()){
            v2_idx=he.next().vertex().getIndex();
            v3_idx=he.next().next().vertex().getIndex();
        
            if(!he.isInterior()){
                continue;
            }
            
            f =he.face();
            Face_area=geometry->faceArea(f);
            Face_normal=geometry->faceNormal(f);
           

            if( rs[v1_idx]> rc && rs[v2_idx]>rc && rs[v3_idx]>rc){
                    // std::cout<<"are we skipping anything?\n";
                    continue;
            }

            unit_r = Unit_rs[v1_idx];
            r = rs[v1_idx];
            unit_r2 = Unit_rs[v2_idx];
            r2 = rs[v2_idx];
            unit_r3 = Unit_rs[v3_idx];
            r3 = rs[v3_idx];

            Vertex_Angle_Normal = Angle_Normals[v1_idx];
            Vertex_Angle_Normal_norm = Angle_Normals_norm[v1_idx];

            Vertex_Angle_Normal2 = Angle_Normals[v2_idx];
            Vertex_Angle_Normal_norm2 = Angle_Normals_norm[v2_idx];

            Vertex_Angle_Normal3 = Angle_Normals[v3_idx];
            Vertex_Angle_Normal_norm3 = Angle_Normals_norm[v3_idx];

            u=geometry->inputVertexPositions[v2_idx]-geometry->inputVertexPositions[v3_idx];
            Vector3 Grad_vec=(0.5)* cross(Normals[f.getIndex()],u);
            
         
            if(Consider_vertex[v1_idx]){
                F+=(1.0/3.0)*E_v[v]*dot(Vertex_Angle_Normal,unit_r)*Grad_vec/Vertex_Angle_Normal_norm;
                }
             
            if(Consider_vertex[v2_idx]){
               Force[v]+=(1.0/3.0)*E_v[v2_idx]*dot(Vertex_Angle_Normal2,unit_r2)*Grad_vec/Vertex_Angle_Normal_norm2;
            }

                if(Consider_vertex[v3_idx]){
               Force[v]+=(1.0/3.0)*E_v[v3_idx]*dot(Vertex_Angle_Normal3,unit_r3)*Grad_vec/Vertex_Angle_Normal_norm3;
                }

            
            dual_area=Dual_areas[v1_idx];
            
             // First vertex 
            if(Consider_vertex[v1_idx]){
                
            grad =-1*dot(cross(geometry->inputVertexPositions[v3_idx]-geometry->inputVertexPositions[v2_idx],Face_normal),unit_r)*Face_normal/(2*Face_area);

            grad2=-1*dot(cross(geometry->inputVertexPositions[v3_idx]-geometry->inputVertexPositions[v2_idx],Face_normal ),Vertex_Angle_Normal )*Face_normal/(2*Face_area);
            
            
            F+=dual_area*E_v[v1_idx]*grad*geometry->angle(he.corner())/Vertex_Angle_Normal_norm;
            
            F+=dual_area*E_v[v1_idx]*geometry->Angle_grad(he.corner(),v,Face_normal)*dot(Face_normal,unit_r)/Vertex_Angle_Normal_norm;
            
            F+=dual_area*E_v[v1_idx]*geometry->angle(he.corner())*(Face_normal-unit_r*dot(unit_r,Face_normal))/(r*Vertex_Angle_Normal_norm);
            
            Total_force-=dual_area*E_v[v1_idx]*geometry->angle(he.corner())*(Face_normal-unit_r*dot(unit_r,Face_normal))/(r*Vertex_Angle_Normal_norm);

            // // Lets add the terms associated with the norm
            
            
            F+=-1*dual_area*E_v[v1_idx]*grad2*geometry->angle(he.corner())*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);
            
            F+=-1*dual_area*E_v[v1_idx]*geometry->Angle_grad(he.corner(),v,Face_normal)*dot(Face_normal,Vertex_Angle_Normal)*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);
            
            }


            // // Second vertex

            
            he2 = he.next();
            if(Consider_vertex[v2_idx]){
            dual_area2=Dual_areas[v2_idx];
            

            
            grad=-1*dot(cross(geometry->inputVertexPositions[v3_idx]-geometry->inputVertexPositions[v2_idx],Face_normal),unit_r2)*Face_normal/(2*Face_area);

            grad2=-1*dot(cross(geometry->inputVertexPositions[v3_idx]-geometry->inputVertexPositions[v2_idx],Face_normal ),Vertex_Angle_Normal2 )*Face_normal/(2*Face_area);


            Force[v1_idx]+=dual_area2*E_v[v2_idx]*grad*geometry->angle(he2.corner())/Vertex_Angle_Normal_norm2;

            
            Force[v1_idx]+=dual_area2*E_v[v2_idx]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,unit_r2)/Vertex_Angle_Normal_norm2;

       
            // Now i add the terms associated with the norm
           
            Force[v1_idx]+=-1*dual_area2*E_v[v2_idx]*grad2*geometry->angle(he2.corner())*dot(Vertex_Angle_Normal2,unit_r2)/pow(Vertex_Angle_Normal_norm2,3);
            // F5+=-1*dual_area2*E_v[he2.vertex()]*grad2*geometry->angle(he2.corner())*dot(Vertex_Angle_Normal2,unit_r)/pow(Vertex_Angle_Normal_norm2,3);

            
            Force[v1_idx]+=-1*dual_area2*E_v[v2_idx]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,Vertex_Angle_Normal2)*dot(Vertex_Angle_Normal2,unit_r2)/pow(Vertex_Angle_Normal_norm2,3);

            // F5+=-1*dual_area2*E_v[he2.vertex()]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,Vertex_Angle_Normal2)*dot(Vertex_Angle_Normal2,unit_r)/pow(Vertex_Angle_Normal_norm2,3);
            }

            // // Third vertex
             he2 = he2.next();
            if(Consider_vertex[v3_idx]){
     
            dual_area3=Dual_areas[v3_idx];

            // if(dot(Vertex_Angle_Normal3,unit_r)>0 && rs[he2.vertex().getIndex()]<1.2){

                
            grad=-1*dot(cross(geometry->inputVertexPositions[v3_idx]-geometry->inputVertexPositions[v2_idx],Face_normal),unit_r3)*Face_normal/(2*Face_area);

            grad2=-1*dot(cross(geometry->inputVertexPositions[v3_idx]-geometry->inputVertexPositions[v2_idx],Face_normal ),Vertex_Angle_Normal3 )*Face_normal/(2*Face_area);

            
            Force[v1_idx]+=dual_area3*E_v[v3_idx]*grad*geometry->angle(he2.corner())/Vertex_Angle_Normal_norm3;
            
            
          
            Force[v1_idx]+=dual_area3*E_v[v3_idx]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,unit_r3)/Vertex_Angle_Normal_norm3;

    
        //    // Now i add the terms associated with the norm

       
            Force[v1_idx]+=-1*dual_area3*E_v[v3_idx]*grad2*geometry->angle(he2.corner())*dot(Vertex_Angle_Normal3,unit_r3)/pow(Vertex_Angle_Normal_norm3,3);

       
     
            Force[v1_idx]+=-1*dual_area3*E_v[v3_idx]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,Vertex_Angle_Normal3)*dot(Vertex_Angle_Normal3,unit_r3)/pow(Vertex_Angle_Normal_norm3,3);

            
            }

        }
        
      
        Force[v]+=F;
        F={0,0,0};
        F2={0,0,0};
        F3={0,0,0};
        F4={0,0,0};
        dual_area=0;
        // THis is between the iteration of the halfedges and the vertex
        


        


        // Force[v]+=F;
        }
        // std::cout<<"\n\n";
        // 
        // std::cout<<"Vertex number "<< v.getIndex()<< "Dot product equal to"<< dot(geometry->vertexNormalAngleWeighted(v),unit_r)<<"distance r ="<< rs[v.getIndex()]<<"SHould i consider?  \t "<< ((dot(geometry->vertexNormalAngleWeighted(v),(Pos-geometry->inputVertexPositions[v]).unit())>0) && (rs[v.getIndex()]<1.2)) <<"\n";
        // }
        // std::cout<<"\n\n\n";
        return Force;


    }


    if(interaction=="Shifted_LJ_Normal_var"){


        int v1_idx;
        int v2_idx;
        int v3_idx;
        double Face_area;
    
        double dual_area2;
        double dual_area3;
        Vector3 Face_normal;
        Vector3 grad;
        Vector3 grad2;
        // Halfedge he;
        Face f;
        Vector3 Vertex_Normal;
        Vector3 Vertex_Angle_Normal;
        Vector3 Vertex_Angle_Normal2;
        Vector3 Vertex_Angle_Normal3;
        double Vertex_Angle_Normal_norm;
        double Vertex_Angle_Normal_norm2;
        double Vertex_Angle_Normal_norm3;

        
        

        Halfedge he2;
        Vector3 unit_r;
        Vector3 unit_r1;
        Vector3 unit_r2;
        Vector3 unit_r3;
        double r;
        double r2;
        double r3;

        Vector3 u;
        rc=1.2;
        rc2=rc*rc;

        // FOr the record, if i only consider one dot product the first term is correct in that vertex 
        // I now want to add the contributions of the neighbors to this point.

        for(Face f : mesh->faces()){

            Face_normal=Normals[f.getIndex()];
            Face_area = geometry->faceArea(f);
            
            for(Halfedge he : f.adjacentHalfedges()){
              
                v1_idx = he.vertex().getIndex();
                v2_idx = he.next().vertex().getIndex();
                v3_idx = he.next().next().vertex().getIndex();
                unit_r = Unit_rs[v1_idx];
                unit_r2 = Unit_rs[v2_idx];
                unit_r3 = Unit_rs[v3_idx];
                if((dot(unit_r,Face_normal)<0 || rs[v1_idx]>1.2 ) && (dot(unit_r2,Face_normal)<0|| rs[v2_idx]>1.2 )&& (dot(unit_r3,Face_normal)<0 || rs[v3_idx]>1.2) ){
                    break;
                }
              
                u=geometry->inputVertexPositions[v2_idx]-geometry->inputVertexPositions[v3_idx];
                Vector3 Grad_vec=(0.5)* cross(Normals[f.getIndex()],u);

                // I need to add the consider restriction here

                if(dot(unit_r,Face_normal)>0){
                    Force[v1_idx]+=(1.0/3.0)*E_v[v1_idx]*dot(Face_normal,unit_r)*Grad_vec;
                }
                if(dot(unit_r2,Face_normal)>0){
                    Force[v1_idx]+=(1.0/3.0)*E_v[v2_idx]*dot(Face_normal,unit_r2)*Grad_vec;
                }
                if(dot(unit_r3,Face_normal)>0){
                    Force[v1_idx]+=(1.0/3.0)*E_v[v3_idx]*dot(Face_normal,unit_r3)*Grad_vec;
                }

                // We have the area gradient correctly now we will do the energy of interaction
                
                r=rs[v1_idx];
                if(rs[v1_idx]<1.2 && dot(unit_r,Face_normal)>0){
                // 
                
                alpha=2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3 );

                F=dot(Face_normal,unit_r)*(-2)*(Face_area/3)*strength*alpha*(  (sigma*sigma/(r*r*r))*pow(( (rc2/(r*r))-1),2 )+ 2*((sigma*sigma/(r*r))-1)*( (rc2/(r*r))-1 )*(rc2/(r*r*r)) )*unit_r;
                // F=(Face_area/3)*dot(Face_normal,unit_r)*strength*alpha*(-2)*(  (sigma*sigma/(r*r*r))*pow(( (rc2/(r*r))-1),2 )+ 2*((sigma*sigma/(r*r))-1)*( (rc2/(r*r))-1 )*(rc2/(r*r*r)) )*unit_r;
                Total_force-=F;
                Force[v1_idx]+=F;
                }


                // Vertex 1
                if(dot(unit_r,Face_normal)>0){
                    grad =-1*dot(cross(geometry->inputVertexPositions[v3_idx]-geometry->inputVertexPositions[v2_idx],Face_normal),unit_r)*Face_normal/(2*Face_area);            
            
                    Force[v1_idx]+=(Face_area/3.0)*E_v[v1_idx]*grad;    
                    Force[v1_idx]+=(Face_area/3.0)*E_v[v1_idx]*(Face_normal-unit_r*dot(unit_r,Face_normal))/(r);
                
                    Total_force-=(Face_area/3.0)*E_v[v1_idx]*(Face_normal-unit_r*dot(unit_r,Face_normal))/(r);
                }
                // Vertex 2
                
                if(dot(unit_r2,Face_normal)>0){
                    grad=-1*dot(cross(geometry->inputVertexPositions[v3_idx]-geometry->inputVertexPositions[v2_idx],Face_normal),unit_r2)*Face_normal/(2*Face_area);

                    Force[v1_idx]+=(Face_area/3.0)*E_v[v2_idx]*grad;
                }
                // Force[v1_idx]+=grad;

                // Vertex 3
                if(dot(unit_r3,Face_normal)>0){
                    grad=-1*dot(cross(geometry->inputVertexPositions[v3_idx]-geometry->inputVertexPositions[v2_idx],Face_normal),unit_r3)*Face_normal/(2*Face_area);
            
                    Force[v1_idx]+=(Face_area/3.0)*E_v[v3_idx]*grad;
                }
                // Force[v1_idx]+=grad;
            
            
                
            }   
            // The thing i want to do now is to get all the terms 
            
            // The energy term is  (area/3)*E[v]*(N . r)
             

        }
        return Force;


    }






       if(interaction=="test_angle_normal_r_normalized_LJ_Full"){

        double Face_area;
    
        double dual_area2;
        double dual_area3;
        Vector3 Face_normal;
        Vector3 grad;
        Vector3 grad2;
        // Halfedge he;
        Face f;
        Vector3 Vertex_Normal;
        Vector3 Vertex_Angle_Normal;
        Vector3 Vertex_Angle_Normal2;
        Vector3 Vertex_Angle_Normal3;
        double Vertex_Angle_Normal_norm;
        double Vertex_Angle_Normal_norm2;
        double Vertex_Angle_Normal_norm3;
        
        Halfedge he2;
        Vector3 unit_r;
        Vector3 unit_r1;
        Vector3 unit_r2;
        Vector3 unit_r3;
        double r;

        Vector3 u;
        rc=1.2;
        rc2=rc*rc;

        // FOr the record, if i only consider one dot product the first term is correct in that vertex 
        // I now want to add the contributions of the neighbors to this point.
        // for(Vertex v: mesh->vertices()){
        
        for(Vertex v : mesh->vertices()){
        
        // Vertex v =mesh->vertex(0);
            Force[v]={0,0,0};
            
            F={0,0,0};
            unit_r= this->Pos- geometry->inputVertexPositions[v];
            r = unit_r.norm();
            unit_r=unit_r.unit();
            dual_area=geometry->barycentricDualArea(v);
            if(v.getIndex()==3){
                std::cout<<"THe distance is "<< r<<"\n";
            }
            if(r<rc){
           
            // unit_r=unit_r.unit();
            
            alpha=2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3 );
            
            
            // First vertex 
            Vertex_Angle_Normal=geometry->vertexNormalAngleWeighted(v).unit();
           
            // std::cout<<"THe dot product is "<<dot(unit_r,Vertex_Angle_Normal)<<"\n";
            // if(dot(Vertex_Angle_Normal,unit_r)>0){
            if(v.getIndex()==3){
                std::cout<<"The force 1 was "<< F <<"\n";
            }
            F+=dual_area*dot(Vertex_Angle_Normal.unit(),unit_r.unit())*strength*alpha*(-2)*(  (sigma*sigma/(r*r*r))*pow(( (rc2/(r*r))-1),2 )+ 2*((sigma*sigma/(r*r))-1)*( (rc2/(r*r))-1 )*(rc2/(r*r*r)) )*unit_r;
            
            // if(dot(Vertex_Angle_Normal,unit_r)>0){
            // Force[v]+=F;

            if(v.getIndex()==3){
                std::cout<<"The force 1 is "<< F <<"\n";
            }
            // }
            
            // }

            }
            Total_force-=F;

            // std::cout<<" Vertex " << v.getIndex()<<" Has this dot product "<< dot(unit_r,Vertex_Angle_Normal)<<" and this radius " << r <<"\n" ;

            F2T=Vector3({0,0,0});
            F3T=Vector3({0,0,0});
            F4T=Vector3({0,0,0});
            F5T=Vector3({0,0,0});

            // Vertex v=mesh->vertex(3);
            for(Halfedge he : v.outgoingHalfedges()){
                
        
            if(!he.isInterior()){
                if(v.getIndex()==3){
                std::cout<<"THis never should happen\n";
                }
                continue;
            }
            
            
            f =he.face();
            Face_area=geometry->faceArea(f);
            Face_normal=geometry->faceNormal(f);
           

            // if( rs[he.vertex()]> rc && rs[he.next().vertex()]>rc && rs[he.next().next().vertex()]>rc){
            //         std::cout<<"are we skipping anything?\n";
            //         continue;
            // }
            unit_r= this->Pos- geometry->inputVertexPositions[v];
            r = unit_r.norm();
            unit_r=unit_r.unit();
            unit_r2= (Pos-geometry->inputVertexPositions[he.next().vertex()]).unit();
            unit_r3=(Pos-geometry->inputVertexPositions[he.next().next().vertex()]).unit();

            Vertex_Angle_Normal=geometry->vertexNormalAngleWeighted(v);
            Vertex_Angle_Normal_norm=Vertex_Angle_Normal.norm();

            Vertex_Angle_Normal2=geometry->vertexNormalAngleWeighted(he.next().vertex());
            Vertex_Angle_Normal_norm2=Vertex_Angle_Normal2.norm();

            Vertex_Angle_Normal3=geometry->vertexNormalAngleWeighted(he.next().next().vertex());
            Vertex_Angle_Normal_norm3=Vertex_Angle_Normal3.norm();

            u=geometry->inputVertexPositions[he.next().vertex()]-geometry->inputVertexPositions[he.next().next().vertex()];
            Vector3 Grad_vec=(0.5)* cross(Normals[f.getIndex()],u);
            
            F2={0,0,0};
            // if(dot(Vertex_Angle_Normal,unit_r)>0 && rs[v.getIndex()]<1.2 ){
                // F2+=(1.0/3.0)*E_v[v]*dot(Vertex_Angle_Normal,unit_r)*Grad_vec/Vertex_Angle_Normal_norm;
                F+=(1.0/3.0)*E_v[v]*dot(Vertex_Angle_Normal,unit_r)*Grad_vec/Vertex_Angle_Normal_norm;
                
                if(rs[v]>1.2){
                    std::cout<<"The radius is "<< rs[v] << "and the energy is different than 0 :O "<< E_v[v]<<"\n";
                }
                // Force[v]+=(1.0/3.0)*E_v[v]*dot(Vertex_Angle_Normal,unit_r)*Grad_vec/Vertex_Angle_Normal_norm;    
            // }
            F4={0,0,0};
          
            // if(dot(Vertex_Angle_Normal2,unit_r2)>0 && rs[he.next().vertex().getIndex()]<1.2){
              
                if(rs[he.next().vertex()]>1.2){
                    std::cout<<"The radius is "<< rs[he.next().vertex()] << "and the energy is different than 0 :v "<< E_v[he.next().vertex()]<<"\n";
                }
               Force[v]+=(1.0/3.0)*E_v[he.next().vertex()]*dot(Vertex_Angle_Normal2,unit_r2)*Grad_vec/Vertex_Angle_Normal_norm2;
           
            // }

            // if(dot(Vertex_Angle_Normal3,unit_r3)>0 && rs[he.next().next().vertex().getIndex()]<1.2){
              
                if(rs[he.next().next().vertex()]>1.2){
                    std::cout<<"The radius is "<< rs[he.next().next().vertex()] << "and the energy is different than 0 :x"<< E_v[he.next().next().vertex()]<<"\n";
                }
               Force[v]+=(1.0/3.0)*E_v[he.next().next().vertex()]*dot(Vertex_Angle_Normal3,unit_r3)*Grad_vec/Vertex_Angle_Normal_norm3;
           
            // }


            unit_r=(Pos-geometry->inputVertexPositions[v]);
            r=unit_r.norm();
            unit_r=unit_r/r;
            // std::cout<<"THe dot product is "<<dot(unit_r,Vertex_Angle_Normal)<< " and the distance is"<< r <<" \n";
            

            
            dual_area=geometry->barycentricDualArea(v);
            
             // First vertex 
            F3={0,0,0};
            
            // if(dot(Vertex_Angle_Normal,unit_r)>0 ){
                
            grad =-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal),unit_r)*Face_normal/(2*Face_area);

            grad2=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal ),Vertex_Angle_Normal )*Face_normal/(2*Face_area);
            

            
            F+=dual_area*E_v[he.vertex()]*grad*geometry->angle(he.corner())/Vertex_Angle_Normal_norm;
            // Force[he.vertex()]+=dual_area*E_v[he.vertex()]*grad*geometry->angle(he.corner())/Vertex_Angle_Normal_norm;
            
            
            F+=dual_area*E_v[he.vertex()]*geometry->Angle_grad(he.corner(),v,Face_normal)*dot(Face_normal,unit_r)/Vertex_Angle_Normal_norm;
            // Force[he.vertex()]+=dual_area*E_v[he.vertex()]*geometry->Angle_grad(he.corner(),v,Face_normal)*dot(Face_normal,unit_r)/Vertex_Angle_Normal_norm;
            
            
            F+=dual_area*E_v[he.vertex()]*geometry->angle(he.corner())*(Face_normal-unit_r*dot(unit_r,Face_normal))/(r*Vertex_Angle_Normal_norm);
            // Force[he.vertex()]+=dual_area*E_v[he.vertex()]*geometry->angle(he.corner())*(Face_normal-unit_r*dot(unit_r,Face_normal))/(r*Vertex_Angle_Normal_norm);

            Total_force-=dual_area*E_v[he.vertex()]*geometry->angle(he.corner())*(Face_normal-unit_r*dot(unit_r,Face_normal))/(r*Vertex_Angle_Normal_norm);

            // // Lets add the terms associated with the norm
            
            
            F+=-1*dual_area*E_v[he.vertex()]*grad2*geometry->angle(he.corner())*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);
            // Force[he.vertex()]+=-1*dual_area*E_v[he.vertex()]*grad2*geometry->angle(he.corner())*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);

            
            F+=-1*dual_area*E_v[he.vertex()]*geometry->Angle_grad(he.corner(),v,Face_normal)*dot(Face_normal,Vertex_Angle_Normal)*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);
            // Force[he.vertex()]+=-1*dual_area*E_v[he.vertex()]*geometry->Angle_grad(he.corner(),v,Face_normal)*dot(Face_normal,Vertex_Angle_Normal)*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);
            
            // // }


            // // Second vertex

            
            he2 = he.next();
            
            dual_area2=geometry->barycentricDualArea(he2.vertex());
         
        //     Vertex_Angle_Normal=geometry->vertexNormalAngleWeighted(he2.vertex());
        //     Vertex_Angle_Normal_norm=Vertex_Angle_Normal.norm();


            unit_r=(Pos-geometry->inputVertexPositions[he2.vertex()]);
            r=unit_r.norm();
            unit_r=unit_r/r;
            

           //if(dot(Vertex_Angle_Normal2,unit_r)>0 && rs[he2.vertex().getIndex()]<1.2){
            
            grad=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal),unit_r)*Face_normal/(2*Face_area);

            grad2=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal ),Vertex_Angle_Normal2 )*Face_normal/(2*Face_area);


            Force[he.vertex()]+=dual_area2*E_v[he2.vertex()]*grad*geometry->angle(he2.corner())/Vertex_Angle_Normal_norm2;
            // F5+=dual_area2*E_v[he2.vertex()]*grad*geometry->angle(he2.corner())/Vertex_Angle_Normal_norm2;

            
            Force[he.vertex()]+=dual_area2*E_v[he2.vertex()]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,unit_r)/Vertex_Angle_Normal_norm2;

            // F5+=dual_area2*E_v[he2.vertex()]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,unit_r)/Vertex_Angle_Normal_norm2;

            // Now i add the terms associated with the norm
           
            Force[he.vertex()]+=-1*dual_area2*E_v[he2.vertex()]*grad2*geometry->angle(he2.corner())*dot(Vertex_Angle_Normal2,unit_r)/pow(Vertex_Angle_Normal_norm2,3);
            // F5+=-1*dual_area2*E_v[he2.vertex()]*grad2*geometry->angle(he2.corner())*dot(Vertex_Angle_Normal2,unit_r)/pow(Vertex_Angle_Normal_norm2,3);

            
            Force[he.vertex()]+=-1*dual_area2*E_v[he2.vertex()]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,Vertex_Angle_Normal2)*dot(Vertex_Angle_Normal2,unit_r)/pow(Vertex_Angle_Normal_norm2,3);

            // F5+=-1*dual_area2*E_v[he2.vertex()]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,Vertex_Angle_Normal2)*dot(Vertex_Angle_Normal2,unit_r)/pow(Vertex_Angle_Normal_norm2,3);
            // }

        //     // // Third vertex
             he2 = he2.next();
            
        //     Vertex_Angle_Normal=geometry->vertexNormalAngleWeighted(he2.vertex());
        //     Vertex_Angle_Normal_norm=Vertex_Angle_Normal.norm();

            unit_r=(Pos-geometry->inputVertexPositions[he2.vertex()]);
            r=unit_r.norm();
            unit_r=unit_r.unit();
            dual_area3=geometry->barycentricDualArea(he2.vertex());

            // if(dot(Vertex_Angle_Normal3,unit_r)>0 && rs[he2.vertex().getIndex()]<1.2){

                
            grad=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal),unit_r)*Face_normal/(2*Face_area);

            grad2=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal ),Vertex_Angle_Normal3 )*Face_normal/(2*Face_area);

            
            Force[he.vertex()]+=dual_area3*E_v[he2.vertex()]*grad*geometry->angle(he2.corner())/Vertex_Angle_Normal_norm3;
            
            // F5+=dual_area3*E_v[he2.vertex()]*grad*geometry->angle(he2.corner())/Vertex_Angle_Normal_norm3;
            
          
            Force[he.vertex()]+=dual_area3*E_v[he2.vertex()]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,unit_r)/Vertex_Angle_Normal_norm3;

            // F5+=dual_area3*E_v[he2.vertex()]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,unit_r)/Vertex_Angle_Normal_norm3;

        //    // Now i add the terms associated with the norm

       
            Force[he.vertex()]+=-1*dual_area3*E_v[he2.vertex()]*grad2*geometry->angle(he2.corner())*dot(Vertex_Angle_Normal3,unit_r)/pow(Vertex_Angle_Normal_norm3,3);

            // F5+=-1*dual_area3*E_v[he2.vertex()]*grad2*geometry->angle(he2.corner())*dot(Vertex_Angle_Normal3,unit_r)/pow(Vertex_Angle_Normal_norm3,3);

     
            Force[he.vertex()]+=-1*dual_area3*E_v[he2.vertex()]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,Vertex_Angle_Normal3)*dot(Vertex_Angle_Normal3,unit_r)/pow(Vertex_Angle_Normal_norm3,3);

            // F5+=-1*dual_area3*E_v[he2.vertex()]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,Vertex_Angle_Normal3)*dot(Vertex_Angle_Normal3,unit_r)/pow(Vertex_Angle_Normal_norm3,3);
            
            // }




        // Force[v]+=F2+F3;
        // F2T+=F2;
        // F3T+=F3;
        // F4T+=F4;
        // F5T+=F5;
        }
        
        if(v.getIndex()==3){
        std::cout<<"The force at this point should be <-0.710611 4.98089 -6.90241> and it is : " << Force[v]<<"\n";
        std::cout<<"I am gonna ad a contribution that should be < 227.246 36.5989 -212.896  > and it is : "<< F <<"\n";
        

        

        }
        Force[v]+=F;
        F={0,0,0};
        F2={0,0,0};
        F3={0,0,0};
        F4={0,0,0};
        dual_area=0;
        // THis is between the iteration of the halfedges and the vertex
        


        


        // Force[v]+=F;
        }
        // std::cout<<"\n\n";
        // 
        // std::cout<<"Vertex number "<< v.getIndex()<< "Dot product equal to"<< dot(geometry->vertexNormalAngleWeighted(v),unit_r)<<"distance r ="<< rs[v.getIndex()]<<"SHould i consider?  \t "<< ((dot(geometry->vertexNormalAngleWeighted(v),(Pos-geometry->inputVertexPositions[v]).unit())>0) && (rs[v.getIndex()]<1.2)) <<"\n";
        // }
        // std::cout<<"\n\n\n";
        return Force;


    }



       if(interaction=="test_angle_normal_r_normalized_LJ"){

        double Face_area;
        double dual_area2;
        double dual_area3;
        Vector3 Face_normal;
        Vector3 grad;
        Vector3 grad2;
        // Halfedge he;
        Face f;
        Vector3 Vertex_Normal;
        Vector3 Vertex_Angle_Normal;
        Vector3 Vertex_Angle_Normal2;
        Vector3 Vertex_Angle_Normal3;
        double Vertex_Angle_Normal_norm;
        Halfedge he2;
        Vector3 unit_r;
        Vector3 unit_r1;
        Vector3 unit_r2;
        Vector3 unit_r3;
        double r;

        Vector3 u;
        rc2=rc*rc;

        // FOr the record, if i only consider one dot product the first term is correct in that vertex 
        // I now want to add the contributions of the neighbors to this point.
        for(Vertex v: mesh->vertices()){
        
        
        // Vertex v =mesh->vertex(6);
            // std::cout<<"Flag1\n";
            F={0,0,0};
            unit_r= this->Pos- geometry->inputVertexPositions[v];
            r= unit_r.norm();
            // std::cout<<"The value of r is "<<r<<"\n";
            if(r<rc){
                
            
            unit_r=unit_r.unit();
            
            alpha=2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3 );
            dual_area=Dual_areas[v.getIndex()];
            

            // First vertex 
            Vertex_Angle_Normal=geometry->vertexNormalAngleWeighted(v).unit();
           
            // std::cout<<"THe dot product is "<<dot(unit_r,Vertex_Angle_Normal)<<"\n";

            F=dot(Vertex_Angle_Normal,unit_r)*strength*alpha*(-2)*(  (sigma*sigma/(r*r*r))*pow(( (rc2/(r*r))-1),2 )+ 2*((sigma*sigma/(r*r))-1)*( (rc2/(r*r))-1 )*(rc2/(r*r*r)) )*unit_r;
            
            }
            // 
            if(dot(Vertex_Angle_Normal,unit_r)>0){
            Force[v]+=F;
            }
            

            // Vertex v=mesh->vertex(3);
            for(Halfedge he : v.outgoingHalfedges()){
                F2={0,0,0};
            if(!he.isInterior()){
                continue;
            }
            f =he.face();
            Face_area=geometry->faceArea(f);
            Face_normal=geometry->faceNormal(f);
           

            // if( rs[he.vertex()]> rc && rs[he.next().vertex()]>rc && rs[he.next().next().vertex()]>rc){
            //         std::cout<<"are we skipping anything?\n";
            //         continue;
            // }

            unit_r2= (Pos-geometry->inputVertexPositions[he.next().vertex()]).unit();
            unit_r3=(Pos-geometry->inputVertexPositions[he.next().next().vertex()]).unit();

            Vertex_Angle_Normal=geometry->vertexNormalAngleWeighted(v);
            Vertex_Angle_Normal_norm=Vertex_Angle_Normal.norm();

            Vertex_Angle_Normal2=geometry->vertexNormalAngleWeighted(he.next().vertex()).unit();
            Vertex_Angle_Normal3=geometry->vertexNormalAngleWeighted(he.next().next().vertex()).unit();


            u=geometry->inputVertexPositions[he.next().vertex()]-geometry->inputVertexPositions[he.next().next().vertex()];
            Vector3 Grad_vec=(0.5)* cross(Normals[f.getIndex()],u);
                
            // if(dot(Vertex_Angle_Normal,unit_r)>0 && rs[v.getIndex()]<rc){
            //     F2+=(1.0/3.0)*E_v[v]*dot(Vertex_Angle_Normal,unit_r)*Grad_vec;    
            // }
          
            // if(dot(Vertex_Angle_Normal2,unit_r2)>0){
            //     F2+=(1.0/3.0)*E_v[he.next().vertex()]*dot(Vertex_Angle_Normal2,unit_r2)*Grad_vec;
            // }

            // if(dot(Vertex_Angle_Normal3,unit_r3)>0){
            //     F2+=(1.0/3.0)*E_v[he.next().next().vertex()]*dot(Vertex_Angle_Normal3,unit_r3)*Grad_vec;
            // }

            // Force[v]+=F2;

            unit_r=(Pos-geometry->inputVertexPositions[v]);
            r=unit_r.norm();
            unit_r=unit_r/r;


            
            
            
             // First vertex 
            if(dot(Vertex_Angle_Normal,unit_r)>0 ){
                
            grad=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal),unit_r)*Face_normal/(2*Face_area);

            grad2=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal ),Vertex_Angle_Normal )*Face_normal/(2*Face_area);
            
            

            Force[he.vertex()]+=E_v[he.vertex()]*grad*geometry->angle(he.corner())/Vertex_Angle_Normal_norm;

            Force[he.vertex()]+=E_v[he.vertex()]*geometry->Angle_grad(he.corner(),v,Face_normal)*dot(Face_normal,unit_r)/Vertex_Angle_Normal_norm;
            
            Force[he.vertex()]+=E_v[he.vertex()]*geometry->angle(he.corner())*(Face_normal-unit_r*dot(unit_r,Face_normal))/(r*Vertex_Angle_Normal_norm);

            // Lets add the terms associated with the norm

            Force[he.vertex()]+=-1*E_v[he.vertex()]*grad2*geometry->angle(he.corner())*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);

            Force[he.vertex()]+=-1*E_v[he.vertex()]*geometry->Angle_grad(he.corner(),v,Face_normal)*dot(Face_normal,Vertex_Angle_Normal)*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);

            }
            // // // Second vertex

            // dual_area2=Dual_areas[he.next().vertex().getIndex()];
         
            he2 = he.next();
            Vertex_Angle_Normal=geometry->vertexNormalAngleWeighted(he2.vertex());
            Vertex_Angle_Normal_norm=Vertex_Angle_Normal.norm();


            unit_r=(Pos-geometry->inputVertexPositions[he2.vertex()]);
            r=unit_r.norm();
            unit_r=unit_r/r;
            if(dot(Vertex_Angle_Normal,unit_r)>0 ){

            grad=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal),unit_r)*Face_normal/(2*Face_area);

            grad2=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal ),Vertex_Angle_Normal )*Face_normal/(2*Face_area);


            Force[he.vertex()]+=E_v[he2.vertex()]*grad*geometry->angle(he2.corner())/Vertex_Angle_Normal_norm;

            Force[he.vertex()]+=E_v[he2.vertex()]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,unit_r)/Vertex_Angle_Normal_norm;

            // Now i add the terms associated with the norm

            Force[he.vertex()]+=-1*E_v[he2.vertex()]*grad2*geometry->angle(he2.corner())*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);

            Force[he.vertex()]+=-1*E_v[he2.vertex()]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,Vertex_Angle_Normal)*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);
            }

            // // Third vertex
             he2 = he2.next();
            
            Vertex_Angle_Normal=geometry->vertexNormalAngleWeighted(he2.vertex());
            Vertex_Angle_Normal_norm=Vertex_Angle_Normal.norm();

            unit_r=(Pos-geometry->inputVertexPositions[he2.vertex()]);
            r=unit_r.norm();
            unit_r=unit_r/r;

            if(dot(Vertex_Angle_Normal,unit_r)>0 ){
            grad=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal),unit_r)*Face_normal/(2*Face_area);

            grad2=-1*dot(cross(geometry->inputVertexPositions[he.next().next().vertex()]-geometry->inputVertexPositions[he.next().vertex()],Face_normal ),Vertex_Angle_Normal )*Face_normal/(2*Face_area);

            
            Force[he.vertex()]+=E_v[he2.vertex()]*grad*geometry->angle(he2.corner())/Vertex_Angle_Normal_norm;
            
            Force[he.vertex()]+=E_v[he2.vertex()]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,unit_r)/Vertex_Angle_Normal_norm;

           // Now i add the terms associated with the norm

            Force[he.vertex()]+=-1*E_v[he2.vertex()]*grad2*geometry->angle(he2.corner())*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);

            Force[he.vertex()]+=-1*E_v[he2.vertex()]*geometry->Angle_grad(he2.corner(),v,Face_normal)*dot(Face_normal,Vertex_Angle_Normal)*dot(Vertex_Angle_Normal,unit_r)/pow(Vertex_Angle_Normal_norm,3);
            
            }





        }
        }
  
        return Force;


    }

    if(interaction == "test_unit_r"){
        // For this i need? 
        double r;
        Vector3 unit_r;
        Vector3 Example=sqrt(3)*Vector3({1,1,1});
        Vector3 grad;


        for(Vertex v : mesh->vertices()){
            unit_r=Pos-geometry->inputVertexPositions[v];
            r=unit_r.norm();
            unit_r=unit_r.unit();
            grad=(Example-unit_r*dot(unit_r,Example))/r;
            Force[v]=grad;
        }

        return Force;


    }




    




    
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
    double rc2=1.2*1.2;
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


    if(interaction=="test_angle"){

        Vertex v_1 =mesh->vertex(1);
        Corner c = v_1.halfedge().corner();
        
        Total_E+=geometry->angle(c);

    }
    
    if(interaction == "test_normal"){
        
        // for(Vertex v : mesh->vertices()){
        for(Face f : mesh->faces()){
            Total_E+=dot(geometry->faceNormal(f),sqrt(3)*Vector3({0,1,0}));
        }
        // }
    }
    if(interaction == "test_angle_normal"){
        
        // Total_E=dot(geometry->vertexNormalAngleWeighted(mesh->vertex(3)),sqrt(3)*Vector3({0,1,0}));
        // for(Vertex v : mesh->vertices()){
        for(Vertex v : mesh->vertices()){
            Total_E+=dot(geometry->vertexNormalAngleWeighted(v),sqrt(3)*Vector3({0,1,0}));
        }
        return Total_E;
        // }
    }
    if(interaction == "test_angle_normal_r"){
        
        Vector3 unit_r;
        for(Vertex v : mesh->vertices()){
            unit_r=(Pos-geometry->inputVertexPositions[v]).unit();
            Total_E+=dot(geometry->vertexNormalAngleWeighted(v),unit_r);
        }
        return Total_E;
        // }
    }

    if(interaction == "test_unit_r"){
        Vector3 unit_r;
        
        for(Vertex v : mesh->vertices()){

            unit_r=(Pos-geometry->inputVertexPositions[v]).unit();
            Total_E+=dot(unit_r,sqrt(3)*Vector3({1,1,1}));

        }
        return Total_E;
    }
    Total_E=0;


        if(interaction == "test_angle_normal_r_normalized")
        {
        double val;
        Vector3 unit_r;
        Vector3 Angle_normal;
        

     
        for(Vertex v : mesh->vertices()){
        // Vertex v = mesh->vertex(3);
            unit_r=(Pos-geometry->inputVertexPositions[v]).unit();
            
            Angle_normal = geometry->vertexNormalAngleWeighted(v);
            // std::cout<<"Dot product is "<< dot(unit_r,Angle_normal)<<"\n";
            if(dot(unit_r,Angle_normal)>0){
                
                Total_E+=dot(unit_r,Angle_normal)/Angle_normal.norm();
                
            }
        }

        return Total_E;
        }
        Total_E=0;
        
          if(interaction == "test_angle_normal_r_normalized_LJ")
        {
        double val;
        Vector3 unit_r;
        Vector3 Angle_normal;
        double r2;
        double rc2=1.2*1.2;
        double alpha;
        for(Vertex v : mesh->vertices()){
        // Vertex v = mesh->vertex(3);
            
            
            unit_r=(Pos-geometry->inputVertexPositions[v]).unit();
            
            Angle_normal = geometry->vertexNormalAngleWeighted(v);
            // std::cout<<"Dot product is "<< dot(unit_r,Angle_normal)<<"\n";
            if(dot(unit_r,Angle_normal)>0){
            
            unit_r=Pos-geometry->inputVertexPositions[v];
            r2=unit_r.norm2();
            unit_r=unit_r.unit();

            if(r2<rc2){
            
            alpha=2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3.0 );
            dual_area=geometry->barycentricDualArea(v);
        
        
            Total_E += dual_area*strength*alpha*( (sigma*sigma/r2)-1  )*pow( (rc2/r2)-1 ,2.0)*dot(unit_r,Angle_normal)/Angle_normal.norm();

            }
        }
        }
        

        return Total_E;
        }
        
        Total_E=0;

        if(interaction == "Shifted_LJ_Normal")
        {
        double val;
        Vector3 unit_r;
        Vector3 Angle_normal;
        double r2;
        double r;
        double rc2=1.2*1.2;
        double alpha;
        for(Vertex v : mesh->vertices()){
        
            Angle_normal = geometry->vertexNormalAngleWeighted(v);
           
            unit_r=Pos-geometry->inputVertexPositions[v];
            r2=unit_r.norm2();
            r=unit_r.norm();
            unit_r=unit_r/r;
            
            if(r<1.2 && dot(unit_r,Angle_normal)>0){
        
            
            alpha=2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3 );
            dual_area=geometry->barycentricDualArea(v);
            // std::cout<<"We are here at some point right?\n";
            Total_E += dual_area*strength*alpha*( (sigma*sigma/r2)-1  )*pow( (rc2/r2)-1 ,2)*dot(unit_r,Angle_normal)/Angle_normal.norm();
            }

        }
        return Total_E;
        
        }
        
                
        if(interaction == "Shifted_LJ_Normal_var")
        {
        double val;
        Vector3 unit_r;
        Vector3 unit_r2;
        Vector3 unit_r3;
        
        Vector3 Angle_normal;
        Vector3 Face_normal;
        double face_area;
        double r2;
        double r;
        double rc2=1.2*1.2;
        double alpha;
        for (Face f : mesh->faces()){
            Face_normal= geometry->faceNormal(f);
            face_area = geometry->faceArea(f);
            for(Vertex v : f.adjacentVertices()){
                unit_r = Pos-geometry->inputVertexPositions[v];
                r=unit_r.norm();
                unit_r=unit_r/r;
                r2=r*r;
                if(r<1.2 && dot(Face_normal,unit_r)>0){
                // if(r<1.2){
                    
                    alpha=2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3 );
                    // Total_E+=(face_area/3.0);
                    // Total_E+=(face_area/3.0)*strength*alpha*( (sigma*sigma/r2)-1  )*pow( (rc2/r2)-1 ,2);
                    Total_E+=(face_area/3.0)*strength*alpha*( (sigma*sigma/r2)-1  )*pow( (rc2/r2)-1 ,2)*dot(Face_normal,unit_r);
                    // Total_E+=dot(Face_normal,unit_r);
                    
                }
            }
        }

        return Total_E;
        
        }
        
        
    
    
    
        Total_E=0;
    
        if(interaction == "test_angle_normal_r_normalized_LJ_Full")
        {
        double val;
        Vector3 unit_r;
        Vector3 Angle_normal;
        double r2;
        double r;
        double rc2=1.2*1.2;
        double alpha;
        for(Vertex v : mesh->vertices()){
        
            Angle_normal = geometry->vertexNormalAngleWeighted(v);
            // std::cout<<"Dot product is "<< dot(unit_r,Angle_normal)<<"\n";
            unit_r=Pos-geometry->inputVertexPositions[v];
            r2=unit_r.norm2();
            r=unit_r.norm();
            unit_r=unit_r/r;
            
            if(r<1.2){
        
            // if(r2<rc2){
            
            alpha=2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3 );
            dual_area=geometry->barycentricDualArea(v);
            
            Total_E += dual_area*strength*alpha*( (sigma*sigma/r2)-1  )*pow( (rc2/r2)-1 ,2)*dot(unit_r,Angle_normal)/Angle_normal.norm();
            
            // }
        }
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
