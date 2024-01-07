// Implement member functions for MeanCurvatureFlow class.
#include "Mem-3dg.h"
#include "Beads.h"
#include <fstream>
#include <omp.h>

#include <chrono>
#include <Eigen/Core>
// #include <geometrycentral/utilities/eigen_interop_helpers.h>
// #include <geometrycentral/utilities/vector3.h>
/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
using namespace std;
Mem3DG::Mem3DG(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
    // velocity = VertexData<Vector3>(*mesh,{0,0,0});
    Old_norm2=0;
    Current_norm2=0;
    H_Vector_0 = VertexData<double>(*mesh,0.0);
    dH_Vector = VertexData<double> (*mesh,0.0);
    // H_target = VertexData<double> (*mesh,0.0);
    system_time=0;
    grad_norm=0.0;
    is_test=false;
    pulling=false;
    

}
Mem3DG::Mem3DG(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, Bead input_Bead) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
    // velocity = VertexData<Vector3>(*mesh,{0,0,0});
    Old_norm2=0;
    Current_norm2=0;
    H_Vector_0 = VertexData<double>(*mesh,0.0);
    dH_Vector = VertexData<double> (*mesh,0.0);
    // H_target = VertexData<double> (*mesh,0.0);
    system_time=0;
    grad_norm=0.0;
    is_test=false;
    Bead_1=input_Bead;
    pulling=false;
}



Mem3DG::Mem3DG(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo,bool test) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
    // velocity = VertexData<Vector3>(*mesh,{0,0,0});
    Old_norm2=0;
    Current_norm2=0;
    H_Vector_0 = VertexData<double>(*mesh,0.0);
    dH_Vector = VertexData<double> (*mesh,0.0);
    // H_target = VertexData<double> (*mesh,0.0);
    system_time=0;
    grad_norm=0.0;
    is_test=test;
    pulling=false;
}
/* 
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
VertexData<Vector3> Mem3DG::buildFlowOperator(double h, double V_bar, double nu, double c0,double P0,double KA,double KB, double Kd)  {

    // Lets get our target area and curvature
    double V=geometry->totalVolume();
    double A_bar=4*PI*pow(3*V_bar/(4*PI*nu),2.0/3.0);
    double H_bar=sqrt(4*PI/A_bar)*c0/2.0; //Coment this with another comment
    
    // size_t index;
    // for(Vertex v : mesh->vertices()){
    //   index=v.getIndex();
    //   H_Vector_0[index]=geometry->scalarMeanCurvature(v)/geometry->barycentricDualArea(v);
    //   dH_Vector[index]=(H_bar-H_Vector_0[index])/50.0;

    // }

    double A=geometry->totalArea();
    double D_P=-P0*(V-V_bar)/V_bar/V_bar;
    double lambda=KA*(A-A_bar )/A_bar;
    
    
    // This lines are for the bunny i need to delete them later
    // lambda=KA;
    // KB=0;
    // return (SurfaceTension(lambda)+OsmoticPressure(D_P));


    return (Bending(H_bar,KB)+OsmoticPressure(D_P)+SurfaceTension(lambda));
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
            Normal={0,0,0};
            // for(Halfedge he: v.outgoingHalfedges()){
            //   Normal+=2*computeHalfedgeMeanCurvatureVector(he);
            // }
            // index=v.getIndex();

            Normal=2*geometry->vertexNormalMeanCurvature(v);
            Force[v]=Normal;
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
        0.25 * cross(geometry->faceNormal(he.face()), EdgeVector );
  
    areaGrad += 0.25 * cross(geometry->faceNormal(he.twin().face()),
                                 EdgeVector2);
  
  return areaGrad/2 ;
}




VertexData<Vector3> Mem3DG::SurfaceGrad() const {

    size_t index;
    Vector3 Normal;
    size_t N_vert=mesh->nVertices();
    VertexData<Vector3> Force(*mesh);
    Vector3 u;
    Halfedge he_grad;
        
    for(Vertex v : mesh->vertices()) {
            Normal={0,0,0};
            Force[v]={0,0,0};
            Force[v]=-1*2*geometry->vertexNormalMeanCurvature(v);
            
        }

    // std::cout<< "THe surface tension force in magnitude is: "<< -1*lambda*sqrt(Force.transpose()*Force) <<"\n";
    return Force;
}








Vector3 Mem3DG::computeHalfedgeGaussianCurvatureVector(Halfedge he) const {
  Vector3 gaussVec{0, 0, 0};
  if (!he.edge().isBoundary()) {
    // gc::Vector3 eji{} = -vecFromHalfedge(he, *vpg);
    gaussVec = 0.5 * geometry->dihedralAngle(he)*( -1* geometry->inputVertexPositions[he.next().vertex()]+geometry->inputVertexPositions[he.vertex()] ).unit();
  }
  else{
    std::cout<<"This mean gaussian curvature doesnt work";
  }
  return gaussVec;
}




Vector3 Mem3DG::dihedralAngleGradient(Halfedge he, Vertex v) const {
    // std::cout<< he.edge().isBoundary();


    double l = geometry->edgeLength(he.edge());



    if (he.edge().isBoundary()) {
    return Vector3{0, 0, 0};
  } else if (he.vertex() == v) { //This is only used for the SIJ_1
    return (geometry->cotan(he.next().next()) *
                geometry->faceNormal(he.face()) +
            geometry->cotan(he.twin().next()) *
                geometry->faceNormal(he.twin().face())) /l;
  } else if (he.next().vertex() == v) { //This is for the firt s term
    return (geometry->cotan(he.twin().next().next()) *
                geometry->faceNormal(he.twin().face()) +
            geometry->cotan(he.next()) *
                geometry->faceNormal(he.face())) /
           l;
  } else if (he.next().next().vertex() == v) { //Este ocurre para el segundo termino
    return (-(geometry->cotan(he.next().next()) +
              geometry->cotan(he.next())) *
            geometry->faceNormal(he.face())) /
           l;
  } else {
    // mem3dg_runtime_error("Unexpected combination of halfedge and vertex!");
    std::cout<< "THe dihedral angle gradient is not working\n";
    return Vector3{0, 0, 0};
  }
std::cout<< "THis is impossible to print\n";
return Vector3{0, 0, 0};
}





VertexData<Vector3> Mem3DG::Bending(double H0,double KB) const {

    
    size_t neigh_index;
    size_t N_vert=mesh->nVertices();

    Vector3 Hij;
    Vector3 Kij;
    Vector3 Sij_1;
    Vector3 Sij_2;

    Vector3 F1={0,0,0};
    Vector3 F2={0,0,0};
    Vector3 F3={0,0,0};
    Vector3 F4={0,0,0};


    Vector3 Position_1;
    Vector3 Position_2;


    VertexData<Vector3> Force(*mesh);




    VertexData<double> Scalar_MC(*mesh,0.0);
    double factor;
    size_t index1;
    for(Vertex v1 : mesh->vertices()) {
        index1=v1.getIndex();
        Scalar_MC[index1]=geometry->scalarMeanCurvature(v1)/geometry->barycentricDualArea(v1);
        }   
    

    
    auto start = chrono::steady_clock::now();
    double H0i;
    double H0j;
    size_t index;
    for(Vertex v: mesh->vertices()){
        F1={0,0,0};
        F2={0,0,0};
        F3={0,0,0};
        F4={0,0,0};
        index=v.getIndex();
        // H0i= (system_time<50? (H_Vector_0[index]+H0)/2.0: H0);
        H0i=H0;
        Position_1=geometry->inputVertexPositions[v];
        for(Halfedge he: v.outgoingHalfedges()){

            
            neigh_index=he.tipVertex().getIndex();
            // H0j= (system_time<50? (H_Vector_0[neigh_index]+H0)/2: H0);
            H0j=H0;

            Position_2=geometry->inputVertexPositions[neigh_index];
      
            


            Kij=computeHalfedgeGaussianCurvatureVector(he);
            // factor=-1*(Scalar_MC[index]-(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0))-1*(Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
            factor=-1*(Scalar_MC[index]-H0i)-1*(Scalar_MC[neigh_index]-H0j);
            
            F1=F1+factor*Kij;


            Hij=2*computeHalfedgeMeanCurvatureVector(he);

            // factor=(1/3.0)*(Scalar_MC[index]*Scalar_MC[index] -(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0)*(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0))+(2.0/3.0)*(Scalar_MC[neigh_index]*Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0)*(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
            factor=(1/3.0)*(Scalar_MC[index]*Scalar_MC[index] -H0i*H0i)+(2.0/3.0)*(Scalar_MC[neigh_index]*Scalar_MC[neigh_index]-H0j*H0j);

            F2=F2+factor*Hij;


            Sij_1 =  geometry->edgeLength(he.edge()) * dihedralAngleGradient(he,he.vertex());

            // factor=-1*(Scalar_MC[index]-(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0));
            factor=-1*(Scalar_MC[index]-H0i);


            F3=F3+factor*Sij_1;

            Sij_2=(geometry->edgeLength(he.twin().edge())*dihedralAngleGradient(he.twin(),he.vertex())+ geometry->edgeLength(he.next().edge())*dihedralAngleGradient(he.next(),he.vertex()) + geometry->edgeLength(he.twin().next().next().edge())*dihedralAngleGradient(he.twin().next().next(), he.vertex()));
            
            // Sij_2=-1*( geometry->cotan(he.next().next())*geometry->faceNormal(he.face()) + geometry->cotan(he.twin())*geometry->faceNormal(he.twin().face()));

            // factor= -1*(Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
            factor= -1*(Scalar_MC[neigh_index]-H0j);
            
            F4=F4+factor*Sij_2;

    }
    



    Force[index]=F1+F2+F3+F4;

    }




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
double Mem3DG::Backtracking(VertexData<Vector3> Force,double D_P,double V_bar,double A_bar,double KA,double KB,double H_bar,bool bead) {
double c1=5e-4;
double rho=0.7;
double alpha=1e-2;
double positionProjection = 0;
double A=geometry->totalArea();
double V=geometry->totalVolume();
double E_Vol = E_Pressure(D_P,V,V_bar);
double E_Sur = E_Surface(KA,A,A_bar);
double E_Ben = E_Bending(H_bar,KB);
double E_Bead = Bead_1.Energy();
double previousE=E_Vol+E_Sur+E_Ben+E_Bead;
double NewE;
VertexData<Vector3> initial_pos(*mesh);
initial_pos= geometry->inputVertexPositions;
Vector3 Bead_init = this->Bead_1.Pos;
// std::cout<<"THe current energy is "<<previousE <<"Is this awful?\n";
// Zeroth iteration
double Projection=0;
Vector3 center; 
    



geometry->inputVertexPositions+=alpha * Force;
geometry->normalize(Vector3({0.0,0.0,0.0}),false);
geometry->refreshQuantities();
center = geometry->centerOfMass();


this->Bead_1.Move_bead(alpha,center);




// std::cout<< geometry->inputVertexPositions[0]<<"and the other "<< initial_pos[0]<<"\n";

for(Vertex v : mesh->vertices()) {
  Projection+=Force[v.getIndex()].norm2();
  } 


grad_norm=Projection;


A=geometry->totalArea();
V=geometry->totalVolume();
E_Vol=E_Pressure(D_P,V,V_bar);
E_Sur=E_Surface(KA,A,A_bar);
E_Ben=E_Bending(H_bar,KB);
E_Bead=Bead_1.Energy();
NewE=E_Vol+E_Sur+E_Ben+E_Bead;

if(std::isnan(E_Vol)){
  std::cout<<"E vol is nan\n";
}
if(std::isnan(E_Sur)){
  std::cout<<"E sur is nan\n";
}
if(std::isnan(E_Ben)){
  std::cout<<"E ben is nan\n";
}


size_t counter=0;
while(true){
  // if(true){
  if( NewE<= previousE - c1*alpha*Projection && Bead_1.Total_force.norm()*alpha<0.1  ) {
    break;

    }
  if(std::isnan(E_Vol)){
      std::cout<<"E vol is nan\n";
    }
  if(std::isnan(E_Sur)){
      std::cout<<"E sur is nan\n";
    }
  if(std::isnan(E_Ben)){
      std::cout<<"E ben is nan\n";
    }
  if(std::isnan(E_Bead)){
      std::cout<<"E bead is nan\n";
    }

  if(std::isnan(NewE)){
      alpha=-1.0;
      break;
    }


  alpha*=rho;
  if(alpha<1e-8){
    std::cout<<"THe timestep got small so the simulation will end \n";
    alpha=-1.0;
    // continue;
    break;
  }
  // for(Vertex vi : mesh->vertices()){
  //   geometry->inputVertexPositions[vi.getIndex()]= initial_pos[vi.getIndex()]+alpha*Force[vi.getIndex()];
  // }
  geometry->inputVertexPositions = initial_pos+alpha*Force;
  geometry->normalize(Vector3({0.0,0.0,0.0}),false);
  geometry->refreshQuantities();
  center = geometry->centerOfMass();

  this->Bead_1.Reset_bead(Bead_init);
  this->Bead_1.Move_bead(alpha,center);

  geometry->refreshQuantities();  
  // std::cout<<"THe old energy is "<< previousE <<"\n";
  // std::cout<<"Alpha is "<< alpha<<"and the new energy is"<< NewE << "\n";
  // std::cout<<"The projection is :"<<Projection<<"\n";
  // // std::cout<<"THe energy changed to"<<NewE<<"\n";
  // std::cout<< "Volume E"<<E_Vol <<"Surface E" << E_Sur <<"\n";
  
  A=geometry->totalArea();
  V=geometry->totalVolume();
  E_Vol=E_Pressure(D_P,V,V_bar);
  E_Sur=E_Surface(KA,A,A_bar);
  E_Ben=E_Bending(H_bar,KB);
  E_Bead=Bead_1.Energy();
  NewE=E_Vol+E_Sur+E_Ben+E_Bead;
  // std::cout<<"THe old energy is "<< previousE <<"\n";
  // std::cout<<"Alpha is "<< alpha<<"and the new energy is"<< NewE << "\n";
  // std::cout<<"The projection is :"<<Projection<<"\n";
  // // std::cout<<"THe energy changed to"<<NewE<<"\n";
  // std::cout<< "Volume E"<<E_Vol <<"Surface E" << E_Sur <<"\n";



}



return alpha;

}


double Mem3DG::Backtracking(VertexData<Vector3> Force,double D_P,double V_bar,double A_bar,double KA,double KB,double H_bar) {
double c1=1e-4;
double rho=0.7;
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
// Zeroth iteration
double Projection=0;

// std::cout<<"THe initial E is "<<previousE<<"\n";
geometry->inputVertexPositions+=alpha * Force;
// std::cout<< geometry->inputVertexPositions[0]<<"and the other "<< initial_pos[0]<<"\n";

for(Vertex v : mesh->vertices()) {
  Projection+=Force[v.getIndex()].norm2();
  } 


grad_norm=Projection;

geometry->refreshQuantities();

A=geometry->totalArea();
V=geometry->totalVolume();
E_Vol=E_Pressure(D_P,V,V_bar);
E_Sur=E_Surface(KA,A,A_bar);
E_Ben=E_Bending(H_bar,KB);

NewE=E_Vol+E_Sur+E_Ben;

if(std::isnan(E_Vol)){
  std::cout<<"E vol is nan\n";
}
if(std::isnan(E_Sur)){
  std::cout<<"E sur is nan\n";
}
if(std::isnan(E_Ben)){
  std::cout<<"E ben is nan\n";
}


size_t counter=0;
while(true){
  // if(true){
  // std::cout<<"THe new energy is "<<NewE <<"\n";
  if( NewE<= previousE - c1*alpha*Projection ) {
    break;

    }
  if(std::isnan(E_Vol)){
  std::cout<<"E vol is nan\n";
    }
if(std::isnan(E_Sur)){
  std::cout<<"E sur is nan\n";
    }
if(std::isnan(E_Ben)){
  std::cout<<"E ben is nan\n";
    }

  if(std::isnan(NewE)){
    std::cout<<"The energy got Nan\n";
    
    alpha=-1.0;
    break;
  }


  alpha*=rho;
  if(alpha<1e-8){
    // std::cout<<"THe timestep got small\n";
    if(system_time<100000){
      // std::cout<<"But the area evolution is not complete yet\n";
      break;
    }
    else{
    std::cout<<"THe simulation will stop because the timestep got smaller than 1e-8 \n";
    alpha=-1.0;
    // continue;
    break;
    }
  }
  // for(Vertex vi : mesh->vertices()){
  //   geometry->inputVertexPositions[vi.getIndex()]= initial_pos[vi.getIndex()]+alpha*Force[vi.getIndex()];
  // }
  geometry->inputVertexPositions = initial_pos+alpha*Force;
  geometry->refreshQuantities();  
  // std::cout<<"THe old energy is "<< previousE <<"\n";
  // std::cout<<"Alpha is "<< alpha<<"and the new energy is"<< NewE << "\n";
  // std::cout<<"The projection is :"<<Projection<<"\n";
  // // std::cout<<"THe energy changed to"<<NewE<<"\n";
  // std::cout<< "Volume E"<<E_Vol <<"Surface E" << E_Sur <<"\n";

  

 
  
  A=geometry->totalArea();
  V=geometry->totalVolume();
  E_Vol=E_Pressure(D_P,V,V_bar);
  E_Sur=E_Surface(KA,A,A_bar);
  E_Ben=E_Bending(H_bar,KB);
  NewE=E_Vol+E_Sur+E_Ben;
  // std::cout<<"THe old energy is "<< previousE <<"\n";
  // std::cout<<"Alpha is "<< alpha<<"and the new energy is"<< NewE << "\n";
  // std::cout<<"The projection is :"<<Projection<<"\n";
  // // std::cout<<"THe energy changed to"<<NewE<<"\n";
  // std::cout<< "Volume E"<<E_Vol <<"Surface E" << E_Sur <<"\n";



}

// for (Edge e : mesh->edges()){
//   std::cout<< e.remesh;

// }
// std::cout<<"\n";
  if(pulling){
    // std::cout<<"THere is pulling right?\n";
  VertexData<Vector3> Horizontal_pull(*mesh,Vector3({1.0,0.0,0.0}));
  
    // std::cout<<"\n";
  geometry->inputVertexPositions+=alpha*pulling_force*No_remesh_list_v*Horizontal_pull;
  geometry->refreshQuantities();
  }
return alpha;

}


VertexData<Vector3> Mem3DG::Project_force(VertexData<Vector3> Force) const{

  VertexData<Vector3> Projected_Force(*mesh,Vector3({0,0,0}));
  Vector3 Normal_v;
  double projection=0;
  for (Vertex v : mesh->vertices()){
    Normal_v=geometry->vertexNormalAreaWeighted(v);
    projection=dot(Force[v],Normal_v);
    // std::cout<<"The projection is negative?"<< (projection<0) <<" \n";
    Projected_Force[v]= (projection>0)*projection*Normal_v;

  }
  return Force;
  // return Projected_Force;
}

/*
 * Performs integration with the bead
 *
 * Input: The timestep <h>.
 * Returns:
 */
double Mem3DG::integrate(double h, double V_bar, double nu, double c0,double P0,double KA,double KB, double Kd,std::ofstream& Sim_data, double time, bool bead,std::ofstream& Bead_data,bool Save_output_data) {
//, Beads Bead_1 

    if(bead){
      Bead_data<<Bead_1.Pos.x <<" "<< Bead_1.Pos.y << " "<< Bead_1.Pos.z <<" \n";

    }
    // Vector<double> Total_force=buildFlowOperator(h,V_bar,nu,c0,P0,KA,KB,Kd);
    VertexData<Vector3> Force(*mesh);
    Force=buildFlowOperator(h,V_bar,nu,c0,P0,KA,KB,Kd);//+Bead_1.Gradient();
    VertexData<Vector3> Bead_force = Bead_1.Gradient();
    // VertexData<Vector3> Bead_force=Project_force(Bead_1.Gradient());
    Force+=Bead_force;




    // std::cout<<"The bead has initial position of  "<<Bead_1.Pos.x  << " "<< Bead_1.Pos.y << " "<< Bead_1.Pos.z <<" \n";
    // std::cout<<"The sigma for the interaction is "<< Bead_1.sigma <<" and the strength of the potential is "<< Bead_1.strength <<" \n";  

    // 
    
    // std::cout<<Bead_1.Pos.x <<" \n";

    // This force is the gradient basically
    double alpha=h;

    // I have the forces for almost everything i just need the bead.

    size_t vindex;
    size_t Nvert=mesh->nVertices();


    // I want to print the Volume, the area, VOl_E, Area_E, Bending_E 
    double V = geometry->totalVolume();
    double A = geometry->totalArea();
    double A_bar=4*PI*pow(3*V_bar/(4*PI*nu),2.0/3.0);
    double H_bar=sqrt(4*PI/A_bar)*c0/2.0; //Coment this with another comment
    double D_P=-1*P0*(V-V_bar)/V_bar/V_bar;
    double lambda=KA*(A-A_bar )/A_bar;    



    double E_Vol = E_Pressure(D_P,V,V_bar);
    double E_Sur = E_Surface(KA,A,A_bar);
    double E_Ben = E_Bending(H_bar,KB);
    double E_bead = Bead_1.Energy();
    
    double backtrackstep;

    // std::cout<<" The position of vertex 120 is "<<geometry->inputVertexPositions[120].x<<" "<<geometry->inputVertexPositions[120].y<< " "<<geometry->inputVertexPositions[120].z<<" \n"; 

    // if(system_time< 50){
    //   backtrackstep=h;
    // }
    // else{
    backtrackstep=Backtracking(Force,D_P,V_bar,A_bar,KA,KB,H_bar,bead);
    // }
    // std::cout<<" The position of vertex 120 is "<<geometry->inputVertexPositions[120].x<<" "<<geometry->inputVertexPositions[120].y<< " "<<geometry->inputVertexPositions[120].z<<" \n";
    
        
    // geometry->inputVertexPositions=geometry->inputVertexPositions+backtrackstep*Force;

    // Bead_1.Move_bead(backtrackstep,center);

    if(Save_output_data || backtrackstep<0){
    Sim_data << V_bar<<" "<< A_bar<<" "<< time <<" "<< V<<" " << A<<" " << E_Vol << " " << E_Sur << " " << E_Ben <<" " << E_bead << " "<< grad_norm<<" " << backtrackstep<< " "<< Bead_1.Pos.x << " "<< Bead_1.Pos.y << " "<< Bead_1.Pos.z <<" \n";
    }
    system_time+=1;
    
    
  
  return backtrackstep;
}

double Mem3DG::integrate(double h, double V_bar, double nu, double c0,double P0,double KA,double KB, double Kd,std::ofstream& Sim_data, double time,bool Save) {


    
    VertexData<Vector3> Force(*mesh);
    Force=buildFlowOperator(h,V_bar,nu,c0,P0,KA,KB,Kd);
    
    
    // This force is the gradient basically
    double alpha=1e-3;

    // I have the forces for almost everything i just need the bead.

    size_t vindex;
    size_t Nvert=mesh->nVertices();


    // I want to print the Volume, the area, VOl_E, Area_E, Bending_E 
    double V = geometry->totalVolume();
    double A = geometry->totalArea();
    double A_bar=4*PI*pow(3*V_bar/(4*PI*nu),2.0/3.0);
    double H_bar=sqrt(4*PI/A_bar)*c0/2.0; //Coment this with another comment
    double D_P=-1*P0*(V-V_bar)/V_bar/V_bar;
    double lambda=KA*(A-A_bar )/A_bar;    


    double E_Vol=E_Pressure(D_P,V,V_bar);
    double E_Sur=E_Surface(KA,A,A_bar);
    double E_Ben=E_Bending(H_bar,KB);
    double backtrackstep;


    // if(system_time< 50){
    //   backtrackstep=h;
    // }
    // else{
    backtrackstep=Backtracking(Force,D_P,V_bar,A_bar,KA,KB,H_bar);
    // }

    if(Save || backtrackstep<1){
    Sim_data << V_bar<<" "<< A_bar<<" "<< time <<" "<< V<<" " << A<<" " << E_Vol << " " << E_Sur << " " << E_Ben << " "<< grad_norm<<" " << backtrackstep<<" \n";
    }
    system_time+=1;
    
    // for (Vertex v : mesh->vertices()) {
    //     vindex=v.getIndex();

    //     // Update= { Total_force[vindex],Total_force[vindex+Nvert],Total_force[vindex+2*Nvert]  };
    //     // Update=h*Force[vindex];
    //     Update=backtrackstep*Force[vindex];
    //     geometry->inputVertexPositions[v] =geometry->inputVertexPositions[v]+ Update ; // placeholder
    // }
    // Force=MeshData<Vertex,Vector3>::fromVector(Delta_x);
    // geometry->inputVertexPositions=geometry->inputVertexPositions+backtrackstep*Force;
    
  
  return backtrackstep;
}







void Mem3DG::Grad_Vol_dx(std::ofstream& Gradient_file,double P0, double V_bar,size_t index) const{
// I want to calculate the gradient of the volume
VertexData<Vector3> initial_pos(*mesh);
initial_pos= geometry->inputVertexPositions;
// VertexData<Vector3> Gradients(*mesh,0.0);
double V=geometry->totalVolume();
double E_vol=0.0;
double E_vol_back=0.0;
double E_vol_front=0.0;
double dr;
Vector3 grad{0.0,0.0,0.0};
double D_P=-P0*(V-V_bar)/V_bar/V_bar;
E_vol=E_Pressure(D_P,V,V_bar);



VertexData<Vector3> Calc_grad=OsmoticPressure(D_P);

Gradient_file<<Calc_grad[index].x <<" "<<Calc_grad[index].y<<" "<< Calc_grad[index].z<<" \n";

Vector<Vector3> Gradients(10);
dr=10.0;
for(size_t exponent=0;exponent<20;exponent++){
  dr=dr/10.0;
  // dr=pow(10,-1*exponent);

  // std::cout<<dr<<" ";
  
  geometry->inputVertexPositions[index]=initial_pos[index]+ Vector3{dr,0,0};
  geometry->refreshQuantities();
  V=geometry->totalVolume();
  D_P=-P0*(V-V_bar)/V_bar/V_bar;
  E_vol_front=E_Pressure(D_P,V,V_bar);
  geometry->inputVertexPositions[index]=initial_pos[index]- Vector3{dr,0,0};
  geometry->refreshQuantities();
  V=geometry->totalVolume();
  D_P=-P0*(V-V_bar)/V_bar/V_bar;
  E_vol_back=E_Pressure(D_P,V,V_bar);
  grad.x=(E_vol_front-E_vol_back)/(2*dr);

  geometry->inputVertexPositions[index]=initial_pos[index]- Vector3{0,dr,0};
  geometry->refreshQuantities();
  V=geometry->totalVolume();
  D_P=-P0*(V-V_bar)/V_bar/V_bar;
  E_vol_back=E_Pressure(D_P,V,V_bar);
  geometry->inputVertexPositions[index]=initial_pos[index]+ Vector3{0,dr,0};
  geometry->refreshQuantities();
  V=geometry->totalVolume();
  D_P=-P0*(V-V_bar)/V_bar/V_bar;
  E_vol_front=E_Pressure(D_P,V,V_bar);
  grad.y=(E_vol_front-E_vol_back)/(2*dr);

  geometry->inputVertexPositions[index]=initial_pos[index]- Vector3{0,0,dr};
  geometry->refreshQuantities();
  V=geometry->totalVolume();
  D_P=-P0*(V-V_bar)/V_bar/V_bar;
  E_vol_back=E_Pressure(D_P,V,V_bar);
  geometry->inputVertexPositions[index]=initial_pos[index]+ Vector3{0,0,dr};
  geometry->refreshQuantities();
  V=geometry->totalVolume();
  D_P=-P0*(V-V_bar)/V_bar/V_bar;
  E_vol_front=E_Pressure(D_P,V,V_bar);
  grad.z=(E_vol_front-E_vol_back)/(2*dr);


  Gradients[exponent]=grad;
  Gradient_file<< Gradients[exponent].x <<" " <<Gradients[exponent].y<<" "<<Gradients[exponent].z <<" \n";

}




// Gradient_file<<"\n";
// This function should just write
std::cout<<"\n";

return ;
}

// The exponent will be 1e-6

VertexData<Vector3> Mem3DG::Grad_Vol(std::ofstream& Gradient_file,double P0, double V_bar,bool Save) const{
// I want to calculate the gradient of the volume
VertexData<Vector3> initial_pos(*mesh);
VertexData<Vector3> Finite_grad(*mesh);
initial_pos= geometry->inputVertexPositions;
// VertexData<Vector3> Gradients(*mesh,0.0);
double V=geometry->totalVolume();
double E_vol=0.0;
double E_vol_back=0.0;
double E_vol_front=0.0;
double dr;
double D_P=-P0*(V-V_bar)/V_bar/V_bar;
double total_grad_finite=0;
double total_grad_theory=0;

Vector3 grad{0.0,0.0,0.0};
Vector3 grad_theory;
Vector3 difference;



E_vol=E_Pressure(D_P,V,V_bar);

VertexData<Vector3> Calc_grad=OsmoticPressure(D_P);

dr=1e-6;

size_t N_vert = mesh->nVertices();
// for(size_t index=0; index<N_vert; index++){
size_t index;
for(Vertex v : mesh->vertices()){
  index=v.getIndex();
  grad_theory=Calc_grad[v];
  geometry->inputVertexPositions[v]=initial_pos[v]+ Vector3{dr,0,0};
  geometry->refreshQuantities();
  V=geometry->totalVolume();
  D_P=-P0*(V-V_bar)/V_bar/V_bar;
  E_vol_front=E_Pressure(D_P,V,V_bar);
  geometry->inputVertexPositions[v]=initial_pos[v]- Vector3{dr,0,0};
  geometry->refreshQuantities();
  V=geometry->totalVolume();
  D_P=-P0*(V-V_bar)/V_bar/V_bar;
  E_vol_back=E_Pressure(D_P,V,V_bar);
  grad.x=(E_vol_front-E_vol_back)/(2*dr);

  geometry->inputVertexPositions[v]=initial_pos[v]- Vector3{0,dr,0};
  geometry->refreshQuantities();
  V=geometry->totalVolume();
  D_P=-P0*(V-V_bar)/V_bar/V_bar;
  E_vol_back=E_Pressure(D_P,V,V_bar);
  geometry->inputVertexPositions[v]=initial_pos[v]+ Vector3{0,dr,0};
  geometry->refreshQuantities();
  V=geometry->totalVolume();
  D_P=-P0*(V-V_bar)/V_bar/V_bar;
  E_vol_front=E_Pressure(D_P,V,V_bar);
  grad.y=(E_vol_front-E_vol_back)/(2*dr);

  geometry->inputVertexPositions[v]=initial_pos[v]- Vector3{0,0,dr};
  geometry->refreshQuantities();
  V=geometry->totalVolume();
  D_P=-P0*(V-V_bar)/V_bar/V_bar;
  E_vol_back=E_Pressure(D_P,V,V_bar);
  geometry->inputVertexPositions[v]=initial_pos[v]+ Vector3{0,0,dr};
  geometry->refreshQuantities();
  V=geometry->totalVolume();
  D_P=-P0*(V-V_bar)/V_bar/V_bar;
  E_vol_front=E_Pressure(D_P,V,V_bar);
  grad.z=(E_vol_front-E_vol_back)/(2*dr);
  if(Save){
  difference= grad+grad_theory;
  Gradient_file<< difference.x <<" "<<difference.y<<" "<< difference.z<<" "<<difference.norm()/grad.norm() <<" " << grad.norm()/grad_theory.norm() <<" \n" ;
  // difference= grad_theory;
  // Gradient_file<< difference.x <<" "<<difference.y<<" "<< difference.z<<" "<<grad_theory.norm()<<" \n" ;
  total_grad_theory+=grad_theory.norm2();
  total_grad_finite+=grad.norm2();
  }
  Finite_grad[v]=-1*grad;
  geometry->inputVertexPositions[v]=initial_pos[v];
  geometry->refreshQuantities();
}
if(Save){
  Gradient_file<< sqrt(total_grad_theory)<<" "<< sqrt(total_grad_finite)<<"\n";
}

return Finite_grad;
}



VertexData<Vector3> Mem3DG::Grad_Area(std::ofstream& Gradient_file, double A_bar,double KA,bool Save) const{
// I want to calculate the gradient of the volume
VertexData<Vector3> initial_pos(*mesh);
VertexData<Vector3> Finite_grad(*mesh);
initial_pos= geometry->inputVertexPositions;
// VertexData<Vector3> Gradients(*mesh,0.0);
double A=geometry->totalArea();


double E_area=0.0;
double E_tot_area=0.0;
double E_area_back=0.0;
double E_area_front=0.0;
double dr;
double lambda=KA*(A-A_bar )/A_bar;
double total_grad_finite=0;
double total_grad_theory=0;

Vector3 grad{0.0,0.0,0.0};
Vector3 grad_theory;
Vector3 difference;

E_area=E_Surface(KA,A,A_bar);


VertexData<Vector3> Calc_grad=SurfaceTension(lambda);

dr=1e-6;

size_t N_vert = mesh->nVertices();
// for(size_t index=0; index<N_vert; index++){
size_t index;
for(Vertex v :mesh->vertices()){
  // A=geometry->totalArea();
  
  
  // E_tot_area=E_Surface(KA,A,A_bar);
  // std::cout<<"THe difference in energy is "<<E_tot_area-E_area<<" \n";
  index=v.getIndex();
  grad_theory=Calc_grad[v];
  geometry->inputVertexPositions[v]=initial_pos[v]+ Vector3{dr,0,0};
  geometry->refreshQuantities();
  A=geometry->totalArea();
  lambda=KA*(A-A_bar )/A_bar;
  E_area_front=E_Surface(KA,A,A_bar);
  geometry->inputVertexPositions[v]=initial_pos[v]- Vector3{dr,0,0};
  geometry->refreshQuantities();
  A=geometry->totalArea();
  lambda=KA*(A-A_bar )/A_bar;
  E_area_back=E_Surface(KA,A,A_bar);
  grad.x=(E_area_front-E_area_back)/(2*dr);

  geometry->inputVertexPositions[v]=initial_pos[v]- Vector3{0,dr,0};
  geometry->refreshQuantities();
  A=geometry->totalArea();
  lambda=KA*(A-A_bar )/A_bar;
  E_area_back=E_Surface(KA,A,A_bar);
  geometry->inputVertexPositions[v]=initial_pos[v]+ Vector3{0,dr,0};
  geometry->refreshQuantities();
  A=geometry->totalArea();
  lambda=KA*(A-A_bar )/A_bar;
  E_area_front=E_Surface(KA,A,A_bar);
  grad.y=(E_area_front-E_area_back)/(2*dr);

  geometry->inputVertexPositions[v]=initial_pos[v]- Vector3{0,0,dr};
  geometry->refreshQuantities();
  A=geometry->totalArea();
  lambda=KA*(A-A_bar )/A_bar;
  E_area_back=E_Surface(KA,A,A_bar);
  geometry->inputVertexPositions[v]=initial_pos[v]+ Vector3{0,0,dr};
  geometry->refreshQuantities();
  A=geometry->totalArea();
  lambda=KA*(A-A_bar )/A_bar;
  E_area_front=E_Surface(KA,A,A_bar);
  grad.z=(E_area_front-E_area_back)/(2*dr);
  if(Save){
  difference= grad+grad_theory;
  Gradient_file<< difference.x/grad.norm() <<" "<<difference.y/grad.norm()<<" "<< difference.z/grad.norm()<<" "<<difference.norm()/grad.norm()<<" " << grad.norm()/grad_theory.norm() <<" \n" ;
  // difference= grad_theory;
  // Gradient_file<< difference.x <<" "<<difference.y<<" "<< difference.z<<" "<<grad_theory.norm()<<" \n" ;
  total_grad_theory+=grad_theory.norm2();
  total_grad_finite+=grad.norm2();
  }
  Finite_grad[v]=-1*grad;
  geometry->inputVertexPositions[v]=initial_pos[v];
  geometry->refreshQuantities();
}
if(Save){
  Gradient_file<< sqrt(total_grad_theory)<<" "<< sqrt(total_grad_finite)<<"\n";
}

return Finite_grad ;
}




VertexData<Vector3> Mem3DG::Grad_Bending(std::ofstream& Gradient_file, double H_bar,double KB, bool Save) {
// I want to calculate the gradient of the volume
VertexData<Vector3> initial_pos(*mesh);
VertexData<Vector3> Finite_grad(*mesh);
initial_pos= geometry->inputVertexPositions;
// VertexData<Vector3> Gradients(*mesh,0.0);
double E_bending=0.0;
double E_bending_back=0.0;
double E_bending_front=0.0;
double dr;


Vector3 grad{0.0,0.0,0.0};
Vector3 grad_theory;
Vector3 difference;

E_bending=E_Bending(H_bar,KB);

double norm_whole_finite=0;
double norm_whole_theory=0;


VertexData<Vector3> Calc_grad=Bending(H_bar,KB);

dr=1e-6;

size_t N_vert = mesh->nVertices();



for(size_t index=0; index<N_vert; index++){
  grad_theory=Calc_grad[index];
  geometry->inputVertexPositions[index]=initial_pos[index]+ Vector3{dr,0,0};
  geometry->refreshQuantities();

  
  E_bending_front=E_Bending(H_bar,KB);
  geometry->inputVertexPositions[index]=initial_pos[index]- Vector3{dr,0,0};
  geometry->refreshQuantities();

  
  E_bending_back=E_Bending(H_bar,KB);
  grad.x=(E_bending_front-E_bending_back)/(2*dr);

  geometry->inputVertexPositions[index]=initial_pos[index]- Vector3{0,dr,0};
  geometry->refreshQuantities();

  
  E_bending_back=E_Bending(H_bar,KB);
  geometry->inputVertexPositions[index]=initial_pos[index]+ Vector3{0,dr,0};
  geometry->refreshQuantities();

  
  E_bending_front=E_Bending(H_bar,KB);
  grad.y=(E_bending_front-E_bending_back)/(2*dr);

  geometry->inputVertexPositions[index]=initial_pos[index]- Vector3{0,0,dr};
  geometry->refreshQuantities();

  
  E_bending_back=E_Bending(H_bar,KB);
  geometry->inputVertexPositions[index]=initial_pos[index]+ Vector3{0,0,dr};
  geometry->refreshQuantities();

  
  E_bending_front=E_Bending(H_bar,KB);
  grad.z=(E_bending_front-E_bending_back)/(2*dr);
  if(Save){
  difference= grad+grad_theory;
  Gradient_file<< difference.x <<" "<<difference.y<<" "<< difference.z<<" "<<difference.norm()/grad.norm() <<" " <<grad.norm()/grad_theory.norm() << " \n" ;
  // difference= grad_theory;
  // Gradient_file<< difference.x <<" "<<difference.y<<" "<< difference.z<<" "<<grad_theory.norm()<<" \n" ;
  norm_whole_finite+=grad.norm2();
  norm_whole_theory+=grad_theory.norm2();
  }
  Finite_grad[index]=-1*grad;
  geometry->inputVertexPositions[index]=initial_pos[index];
  geometry->refreshQuantities();
}

norm_whole_finite=sqrt(norm_whole_finite);
norm_whole_theory=sqrt(norm_whole_theory);
if(Save){
Gradient_file<< norm_whole_finite<<" "<< norm_whole_finite <<" \n";
}



return Finite_grad;
}













void Mem3DG::Grad_Bending_2(std::ofstream& Gradient_file, double H_bar,double KB) {
// I want to calculate the gradient of the volume
VertexData<Vector3> initial_pos(*mesh);
initial_pos= geometry->inputVertexPositions;
// VertexData<Vector3> Gradients(*mesh,0.0);
double E_bending=0.0;
double E_bending_back=0.0;
double E_bending_front=0.0;
double dr;


Vector3 grad{0.0,0.0,0.0};
Vector3 grad_theory;
Vector3 difference;

E_bending=E_Bending(H_bar,KB);

double norm_diff=0;



VertexData<Vector3> Calc_grad=Bending(H_bar,KB);

dr=1e-6;

size_t N_vert = mesh->nVertices();



for(size_t index=0; index<N_vert; index++){
  grad_theory=Calc_grad[index];
  geometry->inputVertexPositions[index]=initial_pos[index]+ Vector3{dr,0,0};
  geometry->refreshQuantities();

  
  E_bending_front=E_Bending(H_bar,KB);
  geometry->inputVertexPositions[index]=initial_pos[index]- Vector3{dr,0,0};
  geometry->refreshQuantities();

  
  E_bending_back=E_Bending(H_bar,KB);
  grad.x=(E_bending_front-E_bending_back)/(2*dr);

  geometry->inputVertexPositions[index]=initial_pos[index]- Vector3{0,dr,0};
  geometry->refreshQuantities();

  
  E_bending_back=E_Bending(H_bar,KB);
  geometry->inputVertexPositions[index]=initial_pos[index]+ Vector3{0,dr,0};
  geometry->refreshQuantities();

  
  E_bending_front=E_Bending(H_bar,KB);
  grad.y=(E_bending_front-E_bending_back)/(2*dr);

  geometry->inputVertexPositions[index]=initial_pos[index]- Vector3{0,0,dr};
  geometry->refreshQuantities();

  
  E_bending_back=E_Bending(H_bar,KB);
  geometry->inputVertexPositions[index]=initial_pos[index]+ Vector3{0,0,dr};
  geometry->refreshQuantities();

  
  E_bending_front=E_Bending(H_bar,KB);
  grad.z=(E_bending_front-E_bending_back)/(2*dr);

  difference= grad+grad_theory;

  // Gradient_file<< difference.x <<" "<<difference.y<<" "<< difference.z<<" \n" ;
  // difference= grad_theory;
  Gradient_file<< difference.x <<" "<<difference.y<<" "<< difference.z<<" "<<difference.norm()/grad_theory.norm()<<" \n" ;
  // norm_whole_finite+=grad.norm2();
  // norm_whole_theory+=grad_theory.norm2();

}

// norm_whole_finite=sqrt(norm_whole_finite);
// norm_whole_theory=sqrt(norm_whole_theory);

// Gradient_file<< norm_whole_finite<<" "<< norm_whole_finite <<" \n";




return ;
}

void Mem3DG::Bending_test(std::ofstream& Analysis_file, double H0,double KB) {

  // THe idea here would be to do the same as with the other function but store things 

    size_t neigh_index;
    size_t N_vert=mesh->nVertices();

    Vector3 Hij;
    Vector3 Kij;
    Vector3 Sij_1;
    Vector3 Sij_2;
    Vector3 Sij_22;

    Vector3 F1={0,0,0};
    Vector3 F2={0,0,0};
    Vector3 F3={0,0,0};
    Vector3 F4={0,0,0};

    Vector3 Position_1;
    Vector3 Position_2;

    VertexData<Vector3> Force(*mesh);

    VertexData<double> Scalar_MC(*mesh,0.0);
    double factor;
    size_t index1;
    for(Vertex v1 : mesh->vertices()) {
        index1=v1.getIndex();
        Scalar_MC[index1]=geometry->scalarMeanCurvature(v1)/geometry->barycentricDualArea(v1);
        }   
    
    
    size_t index;
    for(Vertex v: mesh->vertices()){
        F1={0,0,0};
        F2={0,0,0};
        F3={0,0,0};
        F4={0,0,0};

        index=v.getIndex();
        Position_1=geometry->inputVertexPositions[v];
        for(Halfedge he: v.outgoingHalfedges()){

            neigh_index=he.tipVertex().getIndex();
            Position_2=geometry->inputVertexPositions[neigh_index];
      
      
            Kij=computeHalfedgeGaussianCurvatureVector(he);
            // factor=-1*(Scalar_MC[index]-(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0))-1*(Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
            factor=-1*(Scalar_MC[index]-H0)-1*(Scalar_MC[neigh_index]-H0);
            F1=F1+factor*Kij;

            Hij=2*computeHalfedgeMeanCurvatureVector(he);
            // factor=(1/3.0)*(Scalar_MC[index]*Scalar_MC[index] -(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0)*(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0))+(2.0/3.0)*(Scalar_MC[neigh_index]*Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0)*(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
            factor=(1/3.0)*(Scalar_MC[index]*Scalar_MC[index] -H0*H0)+(2.0/3.0)*(Scalar_MC[neigh_index]*Scalar_MC[neigh_index]-H0*H0);
            F2=F2+factor*Hij;

            Sij_1 =  geometry->edgeLength(he.edge()) * dihedralAngleGradient(he,he.vertex());
            // factor=-1*(Scalar_MC[index]-(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0));
            factor=-1*(Scalar_MC[index]-H0);
            F3=F3+factor*Sij_1;

            Sij_2=(geometry->edgeLength(he.twin().edge())*dihedralAngleGradient(he.twin(),he.vertex())+ geometry->edgeLength(he.next().edge())*dihedralAngleGradient(he.next(),he.vertex()) + geometry->edgeLength(he.twin().next().next().edge())*dihedralAngleGradient(he.twin().next().next(), he.vertex()));

            // Sij_2=-1*( geometry->cotan(he.next().next())*geometry->faceNormal(he.face()) + geometry->cotan(he.twin())*geometry->faceNormal(he.twin().face()));

            // std::cout<<"Norm new "<< Sij_2.norm()<<"Norm old "<< Sij_22.norm()<<"\n";
            // std::cout<<"Norm ratio new/old = "<<Sij_2.norm()/Sij_22.norm()<<" \n";
            // std::cout<<"Relative orientation "<<dot(Sij_2/Sij_2.norm(),Sij_22/Sij_22.norm())<<" \n\n";


            // factor= -1*(Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
            factor= -1*(Scalar_MC[neigh_index]-H0);
            
            F4=F4+factor*Sij_2;

    }
    Analysis_file<<F1.x <<" "<< F1.y<<" "<< F1. z<<" "<< F2.x<<" "<< F2.y<<" "<< F2.z<<" "<< F3.x<<" "<< F3.y<<" "<< F3.z<<" "<< F4.x<<" "<< F4.y<<" "<< F4.z<<" \n";





    // Force[index]=F1+F2+F3+F4;

    }




  return ;
}




double Mem3DG::integrate_finite(double h, double V_bar, double nu, double c0,double P0,double KA,double KB, double Kd,std::ofstream& Sim_data, double time,std::ofstream& Gradient_file_vol,std::ofstream& Gradient_file_area,std::ofstream& Gradient_file_bending,bool Save ) {



    // Vector<double> Total_force=buildFlowOperator(h,V_bar,nu,c0,P0,KA,KB,Kd);
    VertexData<Vector3> Force(*mesh);
    Vector3 Update;
    size_t vindex;
    size_t Nvert=mesh->nVertices();


    // I want to print the Volume, the area, VOl_E, Area_E, Bending_E 
    double V = geometry->totalVolume();
    double A = geometry->totalArea();
    double A_bar=4*PI*pow(3*V_bar/(4*PI*nu),2.0/3.0);
    double H_bar=sqrt(4*PI/A_bar)*c0/2.0; //Coment this with another comment
    double D_P=-1*P0*(V-V_bar)/V_bar/V_bar;
    double lambda=KA*(A-A_bar )/A_bar;    
    // std::cout<<"The forces are calculated here\n";
    Force=Grad_Vol(Gradient_file_vol,P0,V_bar,Save)+ Grad_Area(Gradient_file_area,A_bar,KA,Save)+Grad_Bending(Gradient_file_bending,H_bar,KB,Save);
    // lambda=KA;
    // KB=0;


    double E_Vol=E_Pressure(D_P,V,V_bar);
    double E_Sur=E_Surface(KA,A,A_bar);
    double E_Ben=E_Bending(H_bar,KB);
    double backtrackstep;


    // std::cout<<"Now is the backtracking\n";
    
    backtrackstep=Backtracking(Force,D_P,V_bar,A_bar,KA,KB,H_bar);
    Sim_data << V_bar<<" "<< A_bar<<" "<< time <<" "<< V<<" " << A<<" " << E_Vol << " " << E_Sur << " " << E_Ben << " "<< grad_norm<<" " <<backtrackstep <<" "<<"\n";

    system_time+=1;
    
    for (Vertex v : mesh->vertices()) {
        vindex=v.getIndex();
        Update=backtrackstep*Force[vindex];
        geometry->inputVertexPositions[v] =geometry->inputVertexPositions[v]+ Update ; // placeholder
    }
    
  
  return backtrackstep;
}




VertexData<Vector3> Mem3DG::Grad_Bead(std::ofstream& Gradient_file,bool Save,bool Projection) {

std::cout<<"The interaction to consider is " <<Bead_1.interaction <<"\n";
// I want to calculate the gradient of the volume
VertexData<Vector3> initial_pos(*mesh);
VertexData<Vector3> Finite_grad(*mesh);
initial_pos= geometry->inputVertexPositions;

// VertexData<Vector3> Gradients(*mesh,0.0);
// double V=geometry->totalVolume();
double E_bead=0.0;
double E_bead_back=0.0;
double E_bead_front=0.0;
double dr;
double r_dist;
// double D_P=-P0*(V-V_bar)/V_bar/V_bar;
double total_grad_finite=0;
double total_grad_theory=0;

Vector3 grad{0.0,0.0,0.0};
Vector3 grad_theory;
Vector3 difference;
Vector3 Area_grad;
Vector3 r;
E_bead=Bead_1.Energy();

// E_vol=E_Pressure(D_P,V,V_bar);

VertexData<Vector3> Calc_grad=Bead_1.Gradient();
VertexData<Vector3> Grad_area=SurfaceGrad();
dr=1e-7;

size_t N_vert = mesh->nVertices();
// for(size_t index=0; index<N_vert; index++){
size_t index;
for(Vertex v : mesh->vertices()){
  index=v.getIndex();
  grad_theory=Calc_grad[v];
 
  geometry->inputVertexPositions[v]=initial_pos[v]+ Vector3{dr,0,0};
  geometry->refreshQuantities();
  E_bead_front=Bead_1.Energy();

  geometry->inputVertexPositions[v]=initial_pos[v]- Vector3{dr,0,0};
  geometry->refreshQuantities();
  E_bead_back=Bead_1.Energy();
  
  grad.x=(E_bead_front-E_bead_back)/(2*dr);

  geometry->inputVertexPositions[v]=initial_pos[v]- Vector3{0,dr,0};
  geometry->refreshQuantities();
  E_bead_back=Bead_1.Energy();
  
  geometry->inputVertexPositions[v]=initial_pos[v]+ Vector3{0,dr,0};
  geometry->refreshQuantities();
  E_bead_front=Bead_1.Energy();
  
  grad.y=(E_bead_front-E_bead_back)/(2*dr);

  geometry->inputVertexPositions[v]=initial_pos[v]- Vector3{0,0,dr};
  geometry->refreshQuantities();
  E_bead_back=Bead_1.Energy();

  geometry->inputVertexPositions[v]=initial_pos[v]+ Vector3{0,0,dr};
  geometry->refreshQuantities();
  E_bead_front=Bead_1.Energy();
  
  grad.z=(E_bead_front-E_bead_back)/(2*dr);
  



  if(Save){
  difference= grad+grad_theory;
  
  // I can do this i have grad and grad_theory so i can actually compare them 
  

  // I want to know a little more abt this direction.

    // r= Bead_1.Pos- geometry->inputVertexPositions[v];
    // r_dist=r.norm();
    // r= r.unit();
    // Area_grad=Grad_area[v].unit();
    // // Area_grad= geometry->vertexNormalMeanCurvature(v).unit();
    
    // double E_v=4*1.0*(pow(1.0/r_dist,12)-pow(1.0/r_dist,6));
    // Vector3 F2=E_v *-1*geometry->vertexNormalMeanCurvature(v);

  
  Gradient_file<< difference.x <<" "<<difference.y<<" "<< difference.z<<" "<<difference.norm()/grad.norm() <<" "<<difference.norm() << " " << grad.norm()/grad_theory.norm()<<" \n";//<< dot(r.unit(),difference.unit())<<" " << dot(Area_grad,difference.unit())<<" "<< dot(Area_grad,r.unit())<<" \n";//<< dot(r,HN) <<" \n" ;
  // difference= grad_theory;
  // Gradient_file<< difference.x <<" "<<difference.y<<" "<< difference.z<<" "<<grad_theory.norm()<<" \n" ;
  total_grad_theory+=grad_theory.norm2();
  total_grad_finite+=grad.norm2();
  }
  Finite_grad[v]=-1*grad;
  geometry->inputVertexPositions[v]=initial_pos[v];
  geometry->refreshQuantities();
}
if(Save){
  Gradient_file<< sqrt(total_grad_theory)<<" "<< sqrt(total_grad_finite)<<" "<<sqrt(total_grad_theory/total_grad_finite) << " \n";
}

return Finite_grad;
}






VertexData<Vector3> Mem3DG::Grad_tot_Area(std::ofstream& Gradient_file,bool Save) const{
// I want to calculate the gradient of the volume
VertexData<Vector3> initial_pos(*mesh);
VertexData<Vector3> Finite_grad(*mesh);
initial_pos= geometry->inputVertexPositions;
// VertexData<Vector3> Gradients(*mesh,0.0);
double A=geometry->totalArea();


double E_area=0.0;
double E_tot_area=0.0;
double E_area_back=0.0;
double E_area_front=0.0;
double dr;
// double lambda=KA*(A-A_bar )/A_bar;
double total_grad_finite=0;
double total_grad_theory=0;

Vector3 grad{0.0,0.0,0.0};
Vector3 grad_theory;
Vector3 difference;

E_area=A;


VertexData<Vector3> Calc_grad=SurfaceGrad();

dr=1e-7;

size_t N_vert = mesh->nVertices();
// for(size_t index=0; index<N_vert; index++){
size_t index;
for(Vertex v :mesh->vertices()){
  // A=geometry->totalArea();
  
  
  // E_tot_area=E_Surface(KA,A,A_bar);
  // std::cout<<"THe difference in energy is "<<E_tot_area-E_area<<" \n";
  index=v.getIndex();
  grad_theory=Calc_grad[v];
  geometry->inputVertexPositions[v]=initial_pos[v]+ Vector3{dr,0,0};
  geometry->refreshQuantities();
  A=geometry->totalArea();
  E_area_front=A;
  geometry->inputVertexPositions[v]=initial_pos[v]- Vector3{dr,0,0};
  geometry->refreshQuantities();
  A=geometry->totalArea();
  E_area_back=A;
  grad.x=(E_area_front-E_area_back)/(2*dr);

  geometry->inputVertexPositions[v]=initial_pos[v]- Vector3{0,dr,0};
  geometry->refreshQuantities();
  A=geometry->totalArea();

  E_area_back=A;
  geometry->inputVertexPositions[v]=initial_pos[v]+ Vector3{0,dr,0};
  geometry->refreshQuantities();
  A=geometry->totalArea();
  E_area_front=A;
  grad.y=(E_area_front-E_area_back)/(2*dr);

  geometry->inputVertexPositions[v]=initial_pos[v]- Vector3{0,0,dr};
  geometry->refreshQuantities();
  A=geometry->totalArea();
  E_area_back=A;
  geometry->inputVertexPositions[v]=initial_pos[v]+ Vector3{0,0,dr};
  geometry->refreshQuantities();
  A=geometry->totalArea();
  E_area_front=A;
  grad.z=(E_area_front-E_area_back)/(2*dr);
  if(Save){
  difference= grad+grad_theory;
  Gradient_file<< difference.x/grad.norm() <<" "<<difference.y/grad.norm()<<" "<< difference.z/grad.norm()<<" "<<difference.norm()/grad.norm()<<" " << grad.norm()/grad_theory.norm() <<" \n" ;
  // difference= grad_theory;
  // Gradient_file<< difference.x <<" "<<difference.y<<" "<< difference.z<<" "<<grad_theory.norm()<<" \n" ;
  total_grad_theory+=grad_theory.norm2();
  total_grad_finite+=grad.norm2();
  }
  Finite_grad[v]=-1*grad;
  geometry->inputVertexPositions[v]=initial_pos[v];
  geometry->refreshQuantities();
}
if(Save){
  Gradient_file<< sqrt(total_grad_theory)<<" "<< sqrt(total_grad_finite)<< " "<< sqrt(total_grad_theory/total_grad_finite)<<" \n";
}

return Finite_grad ;
}



void Mem3DG::Grad_Bead_dx(std::ofstream& Gradient_file,bool Save){
// I want to calculate the gradient of the volume
VertexData<Vector3> initial_pos(*mesh);
initial_pos= geometry->inputVertexPositions;
// VertexData<Vector3> Gradients(*mesh,0.0);
// double V=geometry->totalVolume();
double E_bead=0.0;
double E_bead_back=0.0;
double E_bead_front=0.0;
double dr;
Vector3 grad{0.0,0.0,0.0};
E_bead=Bead_1.Energy();
int index=7333;

VertexData<Vector3> Calc_grad=Bead_1.Gradient();

Gradient_file<<Calc_grad[index].x <<" "<<Calc_grad[index].y<<" "<< Calc_grad[index].z<<" \n";

Vector<Vector3> Gradients(20);
dr=10.0;
for(size_t exponent=0;exponent<20;exponent++){
  dr=dr/10.0;
  // dr=pow(10,-1*exponent);

  // std::cout<<dr<<" ";
  
  geometry->inputVertexPositions[index]=initial_pos[index]+ Vector3{dr,0,0};
  geometry->refreshQuantities();
  // V=geometry->totalVolume();
  // D_P=-P0*(V-V_bar)/V_bar/V_bar;
  E_bead_front=Bead_1.Energy();
  geometry->inputVertexPositions[index]=initial_pos[index]- Vector3{dr,0,0};
  geometry->refreshQuantities();
 
  E_bead_back=Bead_1.Energy();
  grad.x=(E_bead_front-E_bead_back)/(2*dr);

  geometry->inputVertexPositions[index]=initial_pos[index]- Vector3{0,dr,0};
  geometry->refreshQuantities();
  // V=geometry->totalVolume();
  // D_P=-P0*(V-V_bar)/V_bar/V_bar;
  E_bead_back=Bead_1.Energy();
  geometry->inputVertexPositions[index]=initial_pos[index]+ Vector3{0,dr,0};
  geometry->refreshQuantities();
  // V=geometry->totalVolume();
  // D_P=-P0*(V-V_bar)/V_bar/V_bar;
  E_bead_front=Bead_1.Energy();
  grad.y=(E_bead_front-E_bead_back)/(2*dr);

  geometry->inputVertexPositions[index]=initial_pos[index]- Vector3{0,0,dr};
  geometry->refreshQuantities();
  E_bead_back=Bead_1.Energy();
  geometry->inputVertexPositions[index]=initial_pos[index]+ Vector3{0,0,dr};
  geometry->refreshQuantities();
  // V=geometry->totalVolume();
  // D_P=-P0*(V-V_bar)/V_bar/V_bar;
  E_bead_front=Bead_1.Energy();
  grad.z=(E_bead_front-E_bead_back)/(2*dr);


  Gradients[exponent]=grad;
  Gradient_file<< Gradients[exponent].x <<" " <<Gradients[exponent].y<<" "<<Gradients[exponent].z <<" \n";

}


}



bool Mem3DG::Area_sanity_check(){
  bool all_positive=true;
  for(Vertex v : mesh->vertices()){
    if(geometry->circumcentricDualArea(v)<0){
      all_positive=false;
    }

  }

  return all_positive;
}