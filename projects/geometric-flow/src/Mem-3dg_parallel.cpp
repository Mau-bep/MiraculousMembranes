// Implement member functions for MeanCurvatureFlow class.
#include "Mem-3dg_parallel.h"
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
Mem3DG_p::Mem3DG_p(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
    // velocity = VertexData<Vector3>(*mesh,{0,0,0});
    Old_norm2=0;
    Current_norm2=0;
    // H_Vector_0 = VertexData<double>(*mesh,0.0);
    // dH_Vector = VertexData<double> (*mesh,0.0);
    system_time=0.0;
    grad_norm=0.0;

}

/* 
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
VertexData<Vector3> Mem3DG_p::buildFlowOperator(double h, double V_bar, double nu, double c0,double P0,double KA,double KB, double Kd) const {

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


VertexData<Vector3> Mem3DG_p::OsmoticPressure(double D_P) const {

    // You have the face normals
    
    // size_t index;
    Vector3 Normal;
    size_t N_vert=mesh->nVertices();

    VertexData<Vector3> Force(*mesh);

    #pragma omp parallel for private(Normal)
    // for(Vertex v : mesh->vertices()) {
      for(size_t index=0;index<N_vert;index++){
        Vertex v = mesh->vertex(index);
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


VertexData<Vector3> Mem3DG_p::SurfaceTension(double lambda) const {

    size_t index;
    Vector3 Normal;
    size_t N_vert=mesh->nVertices();
    VertexData<Vector3> Force(*mesh);
    

    #pragma omp parallel for
    for(size_t index=0;index<N_vert;index++){
          Vertex v = mesh->vertex(index);
          Force[index]=2*geometry->vertexNormalMeanCurvature(v);
        }

    return -1*lambda*Force;
}
Vector3 Mem3DG_p::computeHalfedgeMeanCurvatureVector(Halfedge he) const {
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


Vector3 Mem3DG_p::computeHalfedgeGaussianCurvatureVector(Halfedge he) const {
  Vector3 gaussVec{0, 0, 0};
  if (!he.edge().isBoundary()) {
    // gc::Vector3 eji{} = -vecFromHalfedge(he, *vpg);
    gaussVec = 0.5 * geometry->dihedralAngle(he)*( -1* geometry->inputVertexPositions[he.next().vertex()]+geometry->inputVertexPositions[he.vertex()]).unit();
  }
  else{
    std::cout<<"This mean gaussian curvature doesnt work";
  }
  return gaussVec;
}




Vector3 Mem3DG_p::dihedralAngleGradient(Halfedge he, Vertex v) const {
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






VertexData<Vector3> Mem3DG_p::Bending(double H0,double KB) const {

    
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
    // Vector<double> Vertex_area(N_vert);
    double factor;
    // omp_set_num_threads(4);  
    // size_t index1;
    #pragma omp parallel for
    for(size_t index1=0;index1<N_vert;index1++){
    // for(Vertex v1 : mesh->vertices()) {
        Vertex v1 = mesh->vertex(index1);
        index1=v1.getIndex();
        Scalar_MC[index1]=geometry->scalarMeanCurvature(v1)/geometry->barycentricDualArea(v1);
        // Vertex_area.coeffRef(index)=geometry->barycentricDualArea(v);
        }   
    
    // std::cout<< Scalar_MC<< "\n";
    // std::cout<<"The spontaneous is"<< H0 <<"\n";


    
    auto start = chrono::steady_clock::now();

    #pragma omp parallel for private(F1,F2,F3,F4,Position_1,neigh_index,Position_2,Kij,factor,Hij,Sij_1,Sij_2)
    for(size_t index=0;index<N_vert;index++){
      //  std::cout << "Thread id: " << omp_get_thread_num() << " loop id: " << index <<  std::endl;
      // std::cout<<index<<"\n";
    // size_t index;
    // for(Vertex v: mesh->vertices()){
        F1={0,0,0};
        F2={0,0,0};
        F3={0,0,0};
        F4={0,0,0};
        // F5={0,0,0};
        // F6={0,0,0};
        
        Vertex v =mesh->vertex(index);
        index=v.getIndex();
        Position_1=geometry->inputVertexPositions[v];
        for(Halfedge he: v.outgoingHalfedges()){

            
            neigh_index=he.tipVertex().getIndex();
            Position_2=geometry->inputVertexPositions[neigh_index];
            
            
            Kij=computeHalfedgeGaussianCurvatureVector(he);
            factor=-1*(Scalar_MC[index]-H0)-1*(Scalar_MC[neigh_index]-H0);
            F1=F1+factor*Kij;


            // Hij=0.5*(geometry->cotan(he.twin())+geometry->cotan(he))*e_ji;
            // Lets try the code one 
            Hij=2*computeHalfedgeMeanCurvatureVector(he);

            factor=(1/3.0)*(Scalar_MC[index]*Scalar_MC[index] -H0*H0)+(2.0/3.0)*(Scalar_MC[neigh_index]*Scalar_MC[neigh_index]-H0*H0);

            F2=F2+factor*Hij;

            // Sij_1=  0.5*(geometry->cotan(he.twin().next().next())* geometry->faceNormal(he.twin().face()) + geometry->cotan(he.next()) * geometry->faceNormal(he.face()));
            
            // THis says the implementation, but uses outgoing edges
            // Sij_1=  0.5*(geometry->cotan(he.next().next())* geometry->faceNormal(he.face()) + geometry->cotan(he.twin()) * geometry->faceNormal(he.twin().face()));

            Sij_1 =  geometry->edgeLength(he.edge()) * dihedralAngleGradient(he,he.vertex());



            factor=-1*(Scalar_MC[index]-H0);
            F3=F3+factor*Sij_1;

            Sij_2=(geometry->edgeLength(he.twin().edge())*dihedralAngleGradient(he.twin(),he.vertex())+ geometry->edgeLength(he.next().edge())*dihedralAngleGradient(he.next(),he.vertex()) + geometry->edgeLength(he.twin().next().next().edge())*dihedralAngleGradient(he.twin().next().next(), he.vertex()));



            factor= -1*(Scalar_MC[neigh_index]-H0);
            F4=F4+factor*Sij_2;




    }
    

    
    Force[index]=F1+F2+F3+F4;

    

    }
    auto end = chrono::steady_clock::now();
    

    // std::cout<<"THis loop took "<<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" miliseconds\n";
    

    return KB*Force;
}


// 

double Mem3DG_p::E_Pressure(double P0,double V, double V_bar) const {


    // double V = geometry->totalVolume();
    
    return -1*0.5*P0*(V-V_bar);
}

double Mem3DG_p::E_Surface(double KA,double A, double A_bar) const {


    // return 0.5*KA*A*A;
  
    return 0.5*KA*(A-A_bar)*(A-A_bar)/A_bar;
}


double Mem3DG_p::E_Bending(double H0,double KB) const{
    // size_t index;
    double Eb=0;
    double H;
    size_t N_vert=mesh->nVertices();
    #pragma omp parallel for private(H) reduction(+:Eb)
    for(size_t index=0;index<N_vert;index++){
    // for(Vertex v : mesh->vertices()) {
        Vertex v= mesh->vertex(index);
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
double Mem3DG_p::Backtracking(VertexData<Vector3> Force,double D_P,double V_bar,double A_bar,double KA,double KB,double H_bar) {
double c1=5e-4;
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
size_t N_vert=mesh->nVertices();

geometry->inputVertexPositions+=alpha * Force;
// std::cout<< geometry->inputVertexPositions[0]<<"and the other "<< initial_pos[0]<<"\n";

#pragma omp parallel for reduction(+:Projection)
for(size_t index=0;index<N_vert;index++){
  Vertex v = mesh->vertex(index);
// for(Vertex v : mesh->vertices()) {
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
    alpha=-1.0;
    break;
  }


  alpha*=rho;
  if(alpha<1e-7){
    std::cout<<"THe timestep got small so the simulation will end \n";
    alpha=-1.0;
    break;
  }


  #pragma omp parallel for 
  for(size_t index2=0; index2<N_vert;index2++){
  // for(Vertex vi : mesh->vertices()){
    Vertex vi= mesh->vertex(index2);
    geometry->inputVertexPositions[vi.getIndex()]= initial_pos[vi.getIndex()]+alpha*Force[vi.getIndex()];
  }

  geometry->refreshQuantities();
 
  
  A=geometry->totalArea();
  V=geometry->totalVolume();
  E_Vol=E_Pressure(D_P,V,V_bar);
  E_Sur=E_Surface(KA,A,A_bar);
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
double Mem3DG_p::integrate(double h, double V_bar, double nu, double c0,double P0,double KA,double KB, double Kd,std::ofstream& Sim_data, double time ) {



    // Vector<double> Total_force=buildFlowOperator(h,V_bar,nu,c0,P0,KA,KB,Kd);
    VertexData<Vector3> Force(*mesh);
    Force=buildFlowOperator(h,V_bar,nu,c0,P0,KA,KB,Kd);
    // std::cout<<"WHat is the first force"<<Force[0]<<"\n";
    // std::cout<<"is this true?"<<"\n";
    // if(time==0.0){

    //   std::cout<<"am i always here?\n";
    //   velocity=Force;
    //   // std::cout<<"WHat is the first velocity"<<velocity[0]<<"\n";
    //   for(Vertex v : mesh->vertices()){
    //       Old_norm2+=Force[v.getIndex()].norm2();

    //   }
    // }
    // else { 
    
    //   for(Vertex v : mesh->vertices()){
    //       Current_norm2+=Force[v.getIndex()].norm2();

    //   }
    //   velocity*=Current_norm2/Old_norm2;
    //   velocity+=Force;
    //   // std::cout<<"The current norm2 is "<<Current_norm2<<"\n";
    //   // std::cout<<"The old norm2 is "<<Old_norm2<<"\n";
    //   Old_norm2=Current_norm2;
    //   // Current_norm2=Old_norm2;
    // }




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

    // lambda=KA;
    // KB=0;


    double E_Vol=E_Pressure(D_P,V,V_bar);
    double E_Sur=E_Surface(KA,A,A_bar);
    double E_Ben=E_Bending(H_bar,KB);
    double backtrackstep;

    
    backtrackstep=Backtracking(Force,D_P,V_bar,A_bar,KA,KB,H_bar);
    
    Sim_data << V_bar<<" "<< A_bar<<" "<< time <<" "<< V<<" " << A<<" " << E_Vol << " " << E_Sur << " " << E_Ben << " "<<grad_norm<<" \n";

    

    system_time+=backtrackstep;
    
    #pragma omp parallel for private(Update)
    for(size_t vindex=0;vindex<Nvert;vindex++){

    // for (Vertex v : mesh->vertices()) {
        // vindex=v.getIndex();
        Vertex v=mesh->vertex(vindex);

        // Update= { Total_force[vindex],Total_force[vindex+Nvert],Total_force[vindex+2*Nvert]  };
        // Update=h*Force[vindex];
        Update=backtrackstep*Force[vindex];
        geometry->inputVertexPositions[v] =geometry->inputVertexPositions[v]+ Update ; // placeholder
    }
    
  
  return backtrackstep;
}