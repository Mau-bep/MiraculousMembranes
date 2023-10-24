// Implement member functions for MeanCurvatureFlow class.
#include "Mem-3dg.h"
#include <fstream>
// #include <Eigen/Core>
// #include <geometrycentral/utilities/eigen_interop_helpers.h>
// #include <geometrycentral/utilities/vector3.h>
/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
Mem3DG::Mem3DG(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
    velocity = VertexData<Vector3>(*mesh,{0,0,0});
    Old_norm2=0;
    Current_norm2=0;
    H_Vector_0 = VertexData<double>(*mesh,0.0);
    dH_Vector = VertexData<double> (*mesh,0.0);
    system_time=0.0;
    

}

/* 
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
VertexData<Vector3> Mem3DG::buildFlowOperator(double h, double V_bar, double nu, double c0,double P0,double KA,double KB, double Kd) {

    // Lets get our target area and curvature
    double V=geometry->totalVolume();
    double A_bar=4*PI*pow(3*V_bar/(4*PI*nu),2.0/3.0);
    double H_bar=sqrt(4*PI/A_bar)*c0/2.0; //Coment this with another comment
    size_t index;
    for(Vertex v : mesh->vertices()){
      index=v.getIndex();
      H_Vector_0[index]=geometry->scalarMeanCurvature(v)/geometry->barycentricDualArea(v);
      dH_Vector[index]=(H_bar-H_Vector_0[index])/50.0;

    }




    double A=geometry->totalArea();
    double D_P=-P0*(V-V_bar)/V_bar/V_bar;
    double lambda=KA*(A-A_bar )/A_bar;
    
    
    // This lines are for the bunny i need to delete them later
    // lambda=KA;
    // KB=0;
    // return (SurfaceTension(lambda)+OsmoticPressure(D_P));


    return (Bending(H_bar,KB,Kd)+OsmoticPressure(D_P)+SurfaceTension(lambda));
    // 
    // +SurfaceTension(lambda)
}


VertexData<Vector3> Mem3DG::OsmoticPressure(double D_P) const {

    // You have the face normals
    
    size_t index;
    Vector3 Normal;
    size_t N_vert=mesh->nVertices();

    VertexData<Vector3> Force(*mesh);
    for(Vertex v : mesh->vertices()) {
     // do science here
        index=v.getIndex();
        Normal={0,0,0};
        for(Face f : v.adjacentFaces()) {
            Normal+=geometry->faceArea(f)*geometry->faceNormal(f);

        }
        // Force[v.getIndex()]=D_P*Normal/3.0;
        Force[v.getIndex()]=D_P*Normal/3.0;
    }

    // std::cout<< "THe osmotic pressure force in magnitude is: "<< D_P*sqrt(Force.transpose()*Force) <<"\n";
    return Force;
}


VertexData<Vector3> Mem3DG::SurfaceTension(double lambda) const {

    size_t index;
    Vector3 Normal;
    size_t N_vert=mesh->nVertices();
    VertexData<Vector3> Force(*mesh);
    

        
    for(Vertex v : mesh->vertices()) {

            index=v.getIndex();
            Normal=geometry->vertexNormalMeanCurvature(v);
            Force[index]=Normal;
        }

    // std::cout<< "THe surface tension force in magnitude is: "<< -1*lambda*sqrt(Force.transpose()*Force) <<"\n";
    return -1*lambda*Force;
}
Vector3 Mem3DG::computeHalfedgeMeanCurvatureVector(Halfedge he) const {
   size_t fID = he.face().getIndex();
   size_t fID_he_twin = he.twin().face().getIndex();
   Vector3 areaGrad{0,0,0};

   Vector3 EdgeVector = geometry->inputVertexPositions[he.next().next().vertex()] - geometry->inputVertexPositions[he.next().vertex()];
   Vector3 EdgeVector2 = geometry->inputVertexPositions[he.twin().vertex()] - geometry->inputVertexPositions[he.twin().next().next().vertex()];
   
    areaGrad +=
        0.25 * cross(geometry->faceNormals[fID], EdgeVector );
  
    areaGrad += 0.25 * cross(geometry->faceNormals[fID_he_twin],
                                 EdgeVector2);
  
  return areaGrad / 2;
}


Vector3 Mem3DG::computeHalfedgeGaussianCurvatureVector(Halfedge he) const {
  Vector3 gaussVec{0, 0, 0};
  if (!he.edge().isBoundary()) {
    // gc::Vector3 eji{} = -vecFromHalfedge(he, *vpg);
    gaussVec = 0.5 * geometry->dihedralAngle(he)*( -1* geometry->inputVertexPositions[he.next().vertex()]-geometry->inputVertexPositions[he.vertex()]).unit();
  }
  else{
    std::cout<<"This mean gaussian curvature doesnt work";
  }
  return gaussVec;
}




Vector3 Mem3DG::dihedralAngleGradient(Halfedge he, Vertex v) const {
    // std::cout<< he.edge().isBoundary();


    double l = geometry->edgeLengths[he.edge()];
    // double l = 1.0;

    if (he.edge().isBoundary()) {
    return Vector3{0, 0, 0};
  } else if (he.vertex() == v) {
    return (geometry->cotan(he.next().next()) *
                geometry->faceNormals[he.face()] +
            geometry->cotan(he.twin().next()) *
                geometry->faceNormals[he.twin().face()]) /
           l;
  } else if (he.next().vertex() == v) {
    return (geometry->cotan(he.twin().next().next()) *
                geometry->faceNormals[he.twin().face()] +
            geometry->cotan(he.next()) *
                geometry->faceNormals[he.face()]) /
           l;
  } else if (he.next().next().vertex() == v) {
    return (-(geometry->cotan(he.next().next()) +
              geometry->cotan(he.next())) *
            geometry->faceNormals[he.face()]) /
           l;
  } else {
    // mem3dg_runtime_error("Unexpected combination of halfedge and vertex!");
    std::cout<< "THe dihedral angle gradient is not working\n";
    return Vector3{0, 0, 0};
  }
std::cout<< "THis is impossible to print\n";
return Vector3{0, 0, 0};
}

Vector3 Mem3DG::cornerAngleGradient(Corner c, Vertex v) const {
  Halfedge he = c.halfedge();
  Vector3 n = geometry->faceNormals[c.face()];
  Vector3 ej = geometry->inputVertexPositions[he.next().vertex()] - geometry->inputVertexPositions[he.vertex()];

  Vector3 ei = geometry->inputVertexPositions[he.next().next().vertex()] - geometry->inputVertexPositions[he.next().vertex()];

  Vector3 ek = geometry->inputVertexPositions[he.vertex()] - geometry->inputVertexPositions[he.next().next().vertex()];


  if (c.vertex() == v) { // vi
    Vector3 grad_anglek = -cross(n, ej).normalize() / norm(ej);
    Vector3 grad_anglej = -cross(n, ek).normalize() / norm(ek);
    return -(grad_anglek + grad_anglej);
  } else if (he.next().vertex() == v) { // vk
    return -cross(n, ej).normalize() / norm(ej);
  } else if (he.next().next().vertex() == v) { // vj
    return -cross(n, ek).normalize() / norm(ek);
  } else {
    // mem3dg_runtime_error("Unexpected combination of corner and vertex!");
    std::cout<< "THe corner angle gradient is not working\n";
    return Vector3{0, 0, 0};
  }
}



VertexData<Vector3> Mem3DG::Bending(double H0,double KB, double Kd) const {

    size_t index;
    size_t neigh_index;
    size_t N_vert=mesh->nVertices();
    Vector3 e_ji;
    Vector3 Hij;
    Vector3 Kij;
    Vector3 Sij_1;
    Vector3 Sij_2;

    Vector3 F1={0,0,0};
    Vector3 F2={0,0,0};
    Vector3 F3={0,0,0};
    Vector3 F4={0,0,0};
    // Vector3 F5={0,0,0};
    // Vector3 F6={0,0,0};

    Vector3 Position_1;
    Vector3 Position_2;


    VertexData<Vector3> Force(*mesh);

    Vector<double> Scalar_MC(N_vert);
    // Vector<double> Vertex_area(N_vert);
    double factor;

    for(Vertex v : mesh->vertices()) {
        index=v.getIndex();
        Scalar_MC.coeffRef(index)=geometry->scalarMeanCurvature(v)/geometry->barycentricDualArea(v);
        // Vertex_area.coeffRef(index)=geometry->barycentricDualArea(v);
        }   
    
    // std::cout<< Scalar_MC<< "\n";
    // std::cout<<"The spontaneous is"<< H0 <<"\n";
    for(Vertex v: mesh->vertices()){
        F1={0,0,0};
        F2={0,0,0};
        F3={0,0,0};
        F4={0,0,0};
        // F5={0,0,0};
        // F6={0,0,0};

        index=v.getIndex();
        Position_1=geometry->inputVertexPositions[v];
        for(Halfedge he: v.outgoingHalfedges()){

            
            neigh_index=he.tipVertex().getIndex();
            Position_2=geometry->inputVertexPositions[neigh_index];
            e_ji=Position_1-Position_2;  
          
            Kij=computeHalfedgeGaussianCurvatureVector(he);
            factor=-1*(Scalar_MC[index]-(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0))-1*(Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
            F1=F1+factor*Kij;


            // Hij=0.5*(geometry->cotan(he.twin())+geometry->cotan(he))*e_ji;
            // Lets try the code one 
            Hij=2*computeHalfedgeMeanCurvatureVector(he);

            factor=(1/3.0)*(Scalar_MC[index]*Scalar_MC[index] -(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0)*(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0))+(2.0/3.0)*(Scalar_MC[neigh_index]*Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0)*(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));

            F2=F2+factor*Hij;

            // Sij_1=  0.5*(geometry->cotan(he.twin().next().next())* geometry->faceNormal(he.twin().face()) + geometry->cotan(he.next()) * geometry->faceNormal(he.face()));
            
            // THis says the implementation, but uses outgoing edges
            // Sij_1=  0.5*(geometry->cotan(he.next().next())* geometry->faceNormal(he.face()) + geometry->cotan(he.twin()) * geometry->faceNormal(he.twin().face()));

            Sij_1 = geometry->edgeLengths[he.edge()] * dihedralAngleGradient(he,he.vertex());



            factor=-1*(Scalar_MC[index]-(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0));
            F3=F3+factor*Sij_1;

            Sij_2= (geometry->edgeLengths[he.twin().edge()]*dihedralAngleGradient(he.twin(),he.vertex())+ geometry->edgeLengths[he.next().edge()]*dihedralAngleGradient(he.next(),he.vertex()) + geometry->edgeLengths[he.twin().next().next().edge()]*dihedralAngleGradient(he.twin().next().next(), he.vertex()));


            // Sij_2=-0.5*(   geometry->cotan(he.twin())*geometry->faceNormal(he.twin().face()) + geometry->cotan(he)*geometry->faceNormal(he.face()));
            
            factor= -1*(Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
            F4=F4+factor*Sij_2;



            // Deviatoric_mean=-1*(Kd *Scalar_MC[index]+Kd * Scalar_MC[neigh_index]      )*Kij-1*(Kd *(-1*Scalar_MC[index]*Scalar_MC[index])/3 +Kd*(-1*Scalar_MC[neigh_index]*Scalar_MC[neigh_index])*2.0/3.0 )*Hij-1*(Kd*Scalar_MC[index] * Sij_1 +Kd * Scalar_MC[neigh_index]*Sij_2);

            // Deviatoric_gauss=-1* Kd* cornerAngleGradient(he.corner(),he.vertex())-1*  Kd*cornerAngleGradient(he.next().corner(),he.vertex()) -1* Kd * cornerAngleGradient(he.twin().corner(),he.vertex());

            // F5=Deviatoric_gauss+Deviatoric_mean;

    }
    

    // Force1.coeffRef(index)=F1.x;
    // Force1.coeffRef(index+N_vert)=F1.y;
    // Force1.coeffRef(index+2*N_vert)=F1.z;
    
    Force[index]=F1+F2+F3+F4;

    

    }

    // std::cout<< "THe bending force in magnitude is: "<< KB*sqrt(Force1.transpose()*Force1) <<"\n";

    return KB*Force;
}


// 

double Mem3DG::E_Pressure(double P0,double V, double V_bar) const {


    // double V = geometry->totalVolume();
    
    return -1*0.5*P0*(V-V_bar);
}

double Mem3DG::E_Surface(double KA,double A, double A_bar) const {


    // return 0.5*KA*A*A;
  
    return 0.5*KA*(A-A_bar)*(A-A_bar)/A_bar;
}


double Mem3DG::E_Bending(double H0,double KB) const{
    size_t index;
    double Eb=0;
    double H;
    
    for(Vertex v : mesh->vertices()) {
        index=v.getIndex();
        // Scalar_MC.coeffRef(index)
        H=abs(geometry->scalarMeanCurvature(v)/geometry->barycentricDualArea(v)-H0);
        Eb+=KB*H*H*geometry->barycentricDualArea(v);
        
        }   
    
    return Eb;
}



/*
Performs backtracking
Input:

Returns:


*/
double Mem3DG::Backtracking(VertexData<Vector3> Force,VertexData<Vector3> velocity,double init_time,double D_P,double V_bar,double A_bar,double KA,double KB,double H_bar) const {
double c1=1e-4;
double rho=0.5;
double alpha=1e-3;
double positionProjection = 0;
double A=geometry->totalArea();
double V=geometry->totalVolume();
double E_Vol=E_Pressure(D_P,V,V_bar);
double E_Sur=E_Surface(KA,A,A_bar);
double E_Ben=E_Bending(H_bar,KB);
double previousE=E_Vol+E_Sur+E_Ben;
double NewE;
VertexData<Vector3> initial_pos(*mesh);
initial_pos= geometry->inputVertexPositions;
double time=init_time+alpha;
// Zeroth iteration
double Projection=0;


geometry->inputVertexPositions+=alpha * velocity;
// std::cout<< geometry->inputVertexPositions[0]<<"and the other "<< initial_pos[0]<<"\n";
// for(Vertex v : mesh->vertices()) {
//   Projection+=dot(Force[v.getIndex()],velocity[v.getIndex()]);

// }
// if(Projection<=0){
//   std::cout<<"Is it necessary to update the velocity?\n";
//   Projection=0;
  for(Vertex v : mesh->vertices()) {
  velocity=Force;
  Projection+=Force[v.getIndex()].norm2();
  } 
// }
geometry->refreshQuantities();

A=geometry->totalArea();
V=geometry->totalVolume();
E_Vol=E_Pressure(D_P,V,V_bar);
E_Sur=E_Surface(KA,A,A_bar);
E_Ben=E_Bending(H_bar,KB);
NewE=E_Vol+E_Sur+E_Ben;




// std::cout<<"THe Old E is "<<previousE << "THe new energy is" << NewE<<"and the projection is "<<Projection<<" \n";
// geometry->refreshQuantities();
size_t counter=0;
while(true){
  // if(true){
  if( NewE<= previousE - c1*alpha*Projection ) {
    break;


  }


  // if(alpha<1e-20){
  //   alpha=-1.0;
  //   break;
  // }


  alpha*=rho;

  for(Vertex vi : mesh->vertices()){
    geometry->inputVertexPositions[vi.getIndex()]= initial_pos[vi.getIndex()]+alpha*velocity[vi.getIndex()];
  }

  geometry->refreshQuantities();
  time=init_time+alpha;
  
  A=geometry->totalArea();
  V=geometry->totalVolume();
  E_Vol=abs(E_Pressure(D_P,V,V_bar));
  E_Sur=abs(E_Surface(KA,A,A_bar));
  E_Ben=E_Bending(H_bar,KB);
  NewE=E_Vol+E_Sur+E_Ben;
  // std::cout<<"Alpha is "<< alpha<<"and the new energy is"<< NewE << "\n";
  
  // std::cout<<"THe energy changed to"<<NewE<<"\n";
  



}



return alpha;

}





/*
 * Performs mean curvature flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
double Mem3DG::integrate(double h, double V_bar, double nu, double c0,double P0,double KA,double KB, double Kd,std::ofstream& Sim_data, double time ) {



    // Vector<double> Total_force=buildFlowOperator(h,V_bar,nu,c0,P0,KA,KB,Kd);
    VertexData<Vector3> Force(*mesh);
    Force=buildFlowOperator(h,V_bar,nu,c0,P0,KA,KB,Kd);
    // std::cout<<"WHat is the first force"<<Force[0]<<"\n";
    // std::cout<<"is this true?"<<"\n";
    if(time==0.0){
      // std::cout<<"YES"<<"\n";
      velocity=Force;
      // std::cout<<"WHat is the first velocity"<<velocity[0]<<"\n";
      for(Vertex v : mesh->vertices()){
          Old_norm2+=Force[v.getIndex()].norm2();

      }
    }
    else {
      // std::cout<<"NOPE\n";
      for(Vertex v : mesh->vertices()){
          Current_norm2+=Force[v.getIndex()].norm2();

      }
      velocity*=Current_norm2/Old_norm2;
      velocity+=Force;
    }




    Vector3 Update;
    size_t vindex;
    size_t Nvert=mesh->nVertices();


    // I want to print the Volume, the area, VOl_E, Area_E, Bending_E 
    double V = geometry->totalVolume();
    double A = geometry->totalArea();
    double A_bar=4*PI*pow(3*V_bar/(4*PI*nu),2.0/3.0);
    double H_bar=sqrt(4*PI/A_bar)*c0/2.0; //Coment this with another comment
    double D_P=-1*P0*(V-V_bar);
    double lambda=KA*(A-A_bar )/A_bar;    

    // lambda=KA;
    // KB=0;


    double E_Vol=E_Pressure(D_P,V,V_bar);
    double E_Sur=E_Surface(KA,A,A_bar);
    double E_Ben=E_Bending(H_bar,KB);
    double backtrackstep;

    Sim_data << V_bar<<" "<< A_bar<<" "<< time <<" "<< V<<" " << A<<" " << E_Vol << " " << E_Sur << " " << E_Ben << " "<<"\n";

    
    // if(time<1e-5){
    // backtrackstep=1e-6;
    // }
    // else{
    // std::cout<< "Is backtracking costly?\n";
    backtrackstep=Backtracking(Force,velocity,10,D_P,V_bar,A_bar,KA,KB,H_bar);
    // }
    // std::cout<<"Our job is to find out\n";
    // time+=backtrackstep;
    system_time+=backtrackstep;
    // std::cout<< "EL tstep is"<< backtrackstep<<"\n";

    for (Vertex v : mesh->vertices()) {
        vindex=v.getIndex();

        // Update= { Total_force[vindex],Total_force[vindex+Nvert],Total_force[vindex+2*Nvert]  };
        // Update=h*Force[vindex];
        Update=backtrackstep*Force[vindex];
        geometry->inputVertexPositions[v] =geometry->inputVertexPositions[v]+ Update ; // placeholder
    }
    
  
  return backtrackstep;
}