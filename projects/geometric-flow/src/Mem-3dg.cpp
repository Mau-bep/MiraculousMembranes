// Implement member functions for MeanCurvatureFlow class.
#include "Mem-3dg.h"
#include "Energy_Handler.h"
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
    Area_evol_steps=100000;
    stop_increasing = false;
    small_TS = false;
    recentering = true;
    boundary = false;
    Field="None";
    Field_vals.resize(0);

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
    // Beads.push_back(input_Bead);
    pulling=false;
    Save_SS = false;
    stop_increasing = false;
    small_TS = false;
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
    stop_increasing = false;
    small_TS = false;
}


void Mem3DG::Add_bead(Bead *bead){


  Beads.push_back(bead);


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
    double D_P=-P0*(V-V_bar)/(V_bar*V_bar);
    
    double A_bar=4*PI*pow(3*V_bar/(4*PI*nu),2.0/3.0);
    double H_bar=sqrt(4*PI/A_bar)*c0/2.0; //Coment this with another comment
    double A=geometry->totalArea();
    double lambda=KA*(A-A_bar )/(A_bar*A_bar);
    
    
    // This lines are for the bunny i need to delete them later
    // lambda=KA;
    // KB=0;
    // return (SurfaceTension(lambda)+OsmoticPressure(D_P));


    return (KB*Bending(H_bar)+D_P*OsmoticPressure()+lambda*SurfaceTension());
    // return (Bending(H_bar,KB)+SurfaceTension(lambda));
    
    // 
    // +SurfaceTension(lambda)
}
VertexData<Vector3> Mem3DG::buildFlowOperator(double V_bar, double P0,double KA,double KB, double h)  {

    // Lets get our target area and curvature
    
    double V=geometry->totalVolume();
    double D_P=-P0*(V-V_bar)/V_bar/V_bar;
    
    // double A_bar=4*PI*pow(3*V_bar/(4*PI*nu),2.0/3.0);
    double H_bar=0.0; //Coment this with another comment
    double A=geometry->totalArea();
    double lambda=KA;
    
    
    // This lines are for the bunny i need to delete them later
    // lambda=KA;
    // KB=0;
    // return (SurfaceTension(lambda)+OsmoticPressure(D_P));
    // return (Bending(H_bar,KB)+lambda*SurfaceTension());

    return (KB*Bending(H_bar)+D_P*OsmoticPressure()+lambda*SurfaceTension());
    // return (Bending(H_bar,KB)+SurfaceTension(lambda));
    
    // 
    // +SurfaceTension(lambda)
}

VertexData<Vector3> Mem3DG::buildFlowOperator(double h, double V_bar,double P0,double KA)  {

    // Lets get our target area and curvature
    
    double V=geometry->totalVolume();
    double D_P=-P0*(V-V_bar)/V_bar/V_bar;
 
    
  

    return (D_P*OsmoticPressure()+KA*SurfaceTension());
    // return (Bending(H_bar,KB)+SurfaceTension(lambda));
    
    // 
    // +SurfaceTension(lambda)
}


VertexData<Vector3> Mem3DG::OsmoticPressure() const {

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
        Force[v.getIndex()]=Normal/3.0;
    }

    // std::cout<< "THe osmotic pressure force in magnitude is: "<< D_P*sqrt(Force.transpose()*Force) <<"\n";
    return Force;
}


VertexData<Vector3> Mem3DG::SurfaceTension() const {

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
    return -1*Force;
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
    gaussVec = 0.5 * geometry->dihedralAngle(he)*( -1* geometry->inputVertexPositions[he.next().vertex()]+geometry->inputVertexPositions[he.vertex()] ).unit();
    // std::cout<< "Dihedral angle "<<0.5 * geometry->dihedralAngle(he)<<"\n";
    // std::cout<<" Unit vector of an edge"<< ( -1* geometry->inputVertexPositions[he.next().vertex()]+geometry->inputVertexPositions[he.vertex()] ).unit()<<"\n";
    // std::cout<<"This mean gaussian curvature shouldnt work";
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





VertexData<Vector3> Mem3DG::Bending(double H0) const {

    
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

 

 
    return Force;
}


// 
double Mem3DG::E_Volume_constraint(double KV, double V, double V_bar) const {
  return 0.5*KV*(V-V_bar)*(V-V_bar)/(V_bar*V_bar);
}


double Mem3DG::E_Pressure(double P0,double V, double V_bar) const {


    // double V = geometry->totalVolume();
      
    return -1*0.5*P0*(V-V_bar);
}

double Mem3DG::E_Area_constraint(double KA, double A, double A_bar) const{

  return 0.5*KA*(A-A_bar)*(A-A_bar)/(A_bar*A_bar);
}
double Mem3DG::E_Surface(double KA,double A, double A_bar) const {

    // return 0.5*KA*A*A;
    return 0.5*KA*(A-A_bar)*(A-A_bar)/(A_bar*A_bar);
}


double Mem3DG::E_Bending(double H0,double KB) const{
    size_t index;
    double Eb=0;
    double H;
    double r_eff2;
    Vector3 Pos;
    for(Vertex v : mesh->vertices()) {
        //boundary_fix
        if(v.isBoundary()) continue;
        index=v.getIndex();
        // Scalar_MC.coeffRef(index)
        Pos = geometry->inputVertexPositions[v];
        r_eff2 = Pos.z*Pos.z + Pos.y*Pos.y ;      
        if(r_eff2 > 1.6 && boundary ) continue;
        
        H=abs(geometry->scalarMeanCurvature(v)/geometry->barycentricDualArea(v));
        
        if(std::isnan(H)){
          continue;
          std::cout<<"Dual area: "<< geometry->barycentricDualArea(v);
          std::cout<<"Scalar mean Curv"<< geometry->scalarMeanCurvature(v);
          std::cout<<"One of the H is not a number\n";
        }

        
        

        Eb+=KB*H*H*geometry->barycentricDualArea(v);
        
        }   
    
    return Eb;
}





SparseMatrix<double> Mem3DG::H2_operator(bool CM, bool Vol_const, bool Area_const) {


    // Ok so i have the gradient i want to calculate

    // std::cout<<"We are doing the sobolev operator ! \n";
    // std::cout<<"THe constraints are CM: "<< CM << " Volume " <<  Vol_const <<" and area " << Area_const <<" \n";
    SparseMatrix<double> L = geometry->laplaceMatrix();
    SparseMatrix<double> M = geometry->massMatrix();
    SparseMatrix<double> Inv_M(mesh->nVertices(),mesh->nVertices());
    for(size_t index= 0; index<mesh->nVertices();index++)
    {
        Inv_M.coeffRef(index,index)= 1.0/(M.coeff(index,index));
       
    }

    SparseMatrix<double> J = L.transpose()*Inv_M*L;
    // I need the constraint gradientsx
    VertexData<Vector3> grad_sur = SurfaceGrad();
    VertexData<Vector3> grad_vol = OsmoticPressure();
    // VertexData<Vector3> grad_ben = Bending(0.0);

    // i CAN MAYBE MULTIPLY BY THE SURFACE AREAS 
    

    // std::cout<<"Grad ben \n";

    // std::cout<<"We are debugginggg \n";
    size_t N_vert=mesh->nVertices();


    int Num_constraints = 0;
    if(CM ) Num_constraints+=3;
    if(Vol_const) Num_constraints+=1;
    if(Area_const) Num_constraints+=1;

    // std::cout<<"THe number of constraints is "<< Num_constraints<<" \n";

    SparseMatrix<double> S(N_vert*3+Num_constraints,N_vert*3+Num_constraints);

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    int highest_row = 0;
    int highest_col = 0;
    // We add the constraints first

    int Constraint_number = 0;


    for(size_t index = 0; index<mesh->nVertices();index++)
    {

        Constraint_number = 0;
    //  Area constraint 
        if(Area_const){
          // std::cout<<"Area constraint on \n";
        tripletList.push_back(T(3*index,3*N_vert+Constraint_number, grad_sur[index].x ) );
        tripletList.push_back(T(3*index+1, 3*N_vert+Constraint_number, grad_sur[index].y ) );
        tripletList.push_back(T(3*index+2, 3*N_vert+Constraint_number, grad_sur[index].z ) );
        tripletList.push_back(T(3*N_vert+Constraint_number,3*index, grad_sur[index].x ) );
        tripletList.push_back(T(3*N_vert+Constraint_number, 3*index+1, grad_sur[index].y ) );
        tripletList.push_back(T(3*N_vert+Constraint_number,3*index+2, grad_sur[index].z ) );
        
        Constraint_number++;
        }
        
    // Volume Constraint
        if(Vol_const){
          // std::cout<<"Vol const on\n";
        tripletList.push_back(T(3*index, 3*N_vert+Constraint_number, grad_vol[index].x ) );
        tripletList.push_back(T(3*index+1, 3*N_vert+Constraint_number, grad_vol[index].y ) );
        tripletList.push_back(T(3*index+2, 3*N_vert+Constraint_number, grad_vol[index].z ) );
        tripletList.push_back(T(3*N_vert+Constraint_number, 3*index, grad_vol[index].x ) );
        tripletList.push_back(T(3*N_vert+Constraint_number, 3*index+1, grad_vol[index].y ) );
        tripletList.push_back(T(3*N_vert+Constraint_number, 3*index+2, grad_vol[index].z ) );
        
        Constraint_number++;
        }
    // Position Constraint I
        
        if(CM){
          // std::cout<<"POs constraint on \n";
        tripletList.push_back(T(3*index, 3*N_vert +Constraint_number, 1 ) );
        tripletList.push_back(T(3*N_vert +Constraint_number, 3*index, 1 ) );
        Constraint_number++;
        
        tripletList.push_back(T(3*index + 1, 3*N_vert +Constraint_number, 1 ) );
        tripletList.push_back(T(3*N_vert +Constraint_number, 3*index +1, 1 ) );
        Constraint_number++;
        
        tripletList.push_back(T(3*index + 2 , 3*N_vert+Constraint_number, 1 ) );
        tripletList.push_back(T(3*N_vert+Constraint_number, 3*index + 2, 1 ) );
        Constraint_number++;
        }


    }

    //  std::cout<<"Before setting from tripplets\n";
    // std::cout<<"THe number of vertices is " << N_vert <<"\n";
    // std::cout<<"THe highest expected col is " << 3*N_vert+2 <<" \n";
    // std::cout<<"THe calculated col is " << highest_col <<" \n";
    // std::cout<<"THe highest expected row is " << 3*N_vert+2 <<" \n";
    // std::cout<<"THe calculated row is " << highest_row <<" \n";
    // Now we iterate over the Laplacian
    int row;
    int col;
    double value;

    for( long int k = 0; k < J.outerSize(); ++k ) {
        for( SparseMatrix<double>::InnerIterator it(J,k); it; ++it ) {
            value = it.value();
            row = it.row();
            col = it.col();
            tripletList.push_back(T(3*row,3*col,value));
            tripletList.push_back(T(3*row+1,3*col+1,value));
            tripletList.push_back(T(3*row+2,3*col+2,value));
            // if( 3*row > highest_row) highest_row = 3*row;
            // if( 3*col > highest_col) highest_col = 3*col;
        
        }
    }

    // std::cout<<"Before setting from tripplets\n";
    // std::cout<<"THe number of vertices is " << N_vert <<"\n";
    // std::cout<<"THe highest expected col is " << 3*N_vert+2 <<" \n";
    // std::cout<<"THe calculated col is " << highest_col <<" \n";
    // std::cout<<"THe highest expected row is " << 3*N_vert+2 <<" \n";
    // std::cout<<"THe calculated row is " << highest_row <<" \n";
    S.setFromTriplets(tripletList.begin(),tripletList.end());
    // std::cout<<"THe matrix has been formed\n";

    return S;


    // std::cout<<"AFTER setting from tripplets\n";
    // // We need to create the vector of the RHS

    // Vector<double> RHS(N_vert*3+Num_constraints);
    // int highest_idx = 0;
    // for(size_t index; index<mesh->nVertices();index++)
    // {
    //     RHS.coeffRef(3*index)=grad_ben[index].x;
    //     RHS.coeffRef(3*index+1)=grad_ben[index].y;
    //     RHS.coeffRef(3*index+2)=grad_ben[index].z;
    //     if( 3*index +2 > highest_idx) highest_idx = 3*index+2;
    // }

    // // std::cout<<"The highest index is " << highest_idx <<"\n";
    // // std::cout<<"CREATING RHS\n";
    // for(size_t index = 0; index<Num_constraints;index++)
    // {
    //     RHS.coeffRef(3*N_vert+index)=0;
    // }
    // // RHS.coeffRef(3*N_vert)=0;
    // // RHS.coeffRef(3*N_vert+1)=0;
    // // RHS.coeffRef(3*N_vert+2)=0;
    // // RHS.coeffRef(3*N_vert+3)=0;
    // // RHS.coeffRef(3*N_vert+4)=0;

    // // std::cout<<"RHS IMPLEMENTED\n";

    // // std::cout<<"The RHS is " << RHS << "\n";
    // // We have The RHS and the matrix, its tiiiime 
    // // std::cout<<"Lets solve \n";
    // Eigen::SparseLU<SparseMatrix<double>> solver;
    // // std::cout<<"Solver defined \n";
    
    // solver.compute(S);



    
    // Vector<double> result = solver.solve(RHS);
    // // std::cout<<"Result retrieved\n";
    // VertexData<Vector3> Final_Force(*mesh);
    // for(size_t index; index<mesh->nVertices();index++)
    // {
    //     Final_Force[index]=Vector3{result.coeff(3*index),result.coeff(3*index+1),result.coeff(3*index+2)};
    // } 

    // // std::cout<<"The two lambda for the constraints are "<< result[N_vert] << " and "<< result[N_vert+1] <<"\n";
    // // return Final_Force;
    



}
SparseMatrix<double> Mem3DG::H1_operator(bool CM, bool Vol_const, bool Area_const) {


    // Ok so i have the gradient i want to calculate

    // std::cout<<"We are doing the sobolev operator ! \n";
    // std::cout<<"THe constraints are CM: "<< CM << " Volume " <<  Vol_const <<" and area " << Area_const <<" \n";
    SparseMatrix<double> L = geometry->laplaceMatrix();
    // SparseMatrix<double> M = geometry->massMatrix();
    // SparseMatrix<double> Inv_M(mesh->nVertices(),mesh->nVertices());
    // for(size_t index= 0; index<mesh->nVertices();index++)
    // {
    //     Inv_M.coeffRef(index,index)= 1.0/(M.coeff(index,index));
       
    // }

    SparseMatrix<double> J = L;
    // I need the constraint gradientsx
    VertexData<Vector3> grad_sur = SurfaceGrad();
    VertexData<Vector3> grad_vol = OsmoticPressure();
    // VertexData<Vector3> grad_ben = Bending(0.0);

    // i CAN MAYBE MULTIPLY BY THE SURFACE AREAS 
    

    // std::cout<<"Grad ben \n";

    // std::cout<<"We are debugginggg \n";
    size_t N_vert=mesh->nVertices();


    int Num_constraints = 0;
    if(CM ) Num_constraints+=3;
    if(Vol_const) Num_constraints+=1;
    if(Area_const) Num_constraints+=1;

    // std::cout<<"THe number of constraints is "<< Num_constraints<<" \n";

    SparseMatrix<double> S(N_vert*3+Num_constraints,N_vert*3+Num_constraints);

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    int highest_row = 0;
    int highest_col = 0;
    // We add the constraints first

    int Constraint_number = 0;


    for(size_t index = 0; index<mesh->nVertices();index++)
    {

        Constraint_number = 0;
    //  Area constraint 
        if(Area_const){
          // std::cout<<"Area constraint on \n";
        tripletList.push_back(T(3*index,3*N_vert+Constraint_number, grad_sur[index].x ) );
        tripletList.push_back(T(3*index+1, 3*N_vert+Constraint_number, grad_sur[index].y ) );
        tripletList.push_back(T(3*index+2, 3*N_vert+Constraint_number, grad_sur[index].z ) );
        tripletList.push_back(T(3*N_vert+Constraint_number,3*index, grad_sur[index].x ) );
        tripletList.push_back(T(3*N_vert+Constraint_number, 3*index+1, grad_sur[index].y ) );
        tripletList.push_back(T(3*N_vert+Constraint_number,3*index+2, grad_sur[index].z ) );
        
        Constraint_number++;
        }
        
    // Volume Constraint
        if(Vol_const){
          // std::cout<<"Vol const on\n";
        tripletList.push_back(T(3*index, 3*N_vert+Constraint_number, grad_vol[index].x ) );
        tripletList.push_back(T(3*index+1, 3*N_vert+Constraint_number, grad_vol[index].y ) );
        tripletList.push_back(T(3*index+2, 3*N_vert+Constraint_number, grad_vol[index].z ) );
        tripletList.push_back(T(3*N_vert+Constraint_number, 3*index, grad_vol[index].x ) );
        tripletList.push_back(T(3*N_vert+Constraint_number, 3*index+1, grad_vol[index].y ) );
        tripletList.push_back(T(3*N_vert+Constraint_number, 3*index+2, grad_vol[index].z ) );
        
        Constraint_number++;
        }
    // Position Constraint I
        
        if(CM){
          // std::cout<<"POs constraint on \n";
        tripletList.push_back(T(3*index, 3*N_vert +Constraint_number, 1 ) );
        tripletList.push_back(T(3*N_vert +Constraint_number, 3*index, 1 ) );
        Constraint_number++;
        
        tripletList.push_back(T(3*index + 1, 3*N_vert +Constraint_number, 1 ) );
        tripletList.push_back(T(3*N_vert +Constraint_number, 3*index +1, 1 ) );
        Constraint_number++;
        
        tripletList.push_back(T(3*index + 2 , 3*N_vert+Constraint_number, 1 ) );
        tripletList.push_back(T(3*N_vert+Constraint_number, 3*index + 2, 1 ) );
        Constraint_number++;
        }


    }

    //  std::cout<<"Before setting from tripplets\n";
    // std::cout<<"THe number of vertices is " << N_vert <<"\n";
    // std::cout<<"THe highest expected col is " << 3*N_vert+2 <<" \n";
    // std::cout<<"THe calculated col is " << highest_col <<" \n";
    // std::cout<<"THe highest expected row is " << 3*N_vert+2 <<" \n";
    // std::cout<<"THe calculated row is " << highest_row <<" \n";
    // Now we iterate over the Laplacian
    int row;
    int col;
    double value;

    for( long int  k = 0; k < J.outerSize(); ++k ) {
        for( SparseMatrix<double>::InnerIterator it(J,k); it; ++it ) {
            value = it.value();
            row = it.row();
            col = it.col();
            tripletList.push_back(T(3*row,3*col,value));
            tripletList.push_back(T(3*row+1,3*col+1,value));
            tripletList.push_back(T(3*row+2,3*col+2,value));
            // if( 3*row > highest_row) highest_row = 3*row;
            // if( 3*col > highest_col) highest_col = 3*col;
        
        }
    }

    S.setFromTriplets(tripletList.begin(),tripletList.end());

    return S;


}



VertexData<Vector3> Mem3DG::Linear_force_field(double x0,double slope) const{
  // Can i make it volume preserving?

  VertexData<Vector3> FField(*mesh,Vector3({0,0,0}));
  VertexData<Vector3> FField_correction(*mesh,Vector3({0,0,0}));
  Vector3 Normal;
  double mag = 0;
  double tot_mag = 0;
  for( Vertex v : mesh->vertices()){
    // I want to iterate over all vertices
    Normal = geometry->vertexNormalAreaWeighted(v);
    mag = abs(geometry->inputVertexPositions[v].x-x0)*slope;
    FField[v] = mag*Normal;
    tot_mag += mag;
    FField_correction[v] = Normal;
  }

  // int index;
  // for(Vertex v : mesh->vertices()) {
  //    // do science here
  //       index=v.getIndex();
  //       Normal={0,0,0};
  //       for(Face f : v.adjacentFaces()) {
  //           Normal+=geometry->faceArea(f)*geometry->faceNormal(f);

  //       }
  //       // Force[v.getIndex()]=D_P*Normal/3.0;
  //       mag = abs(geometry->inputVertexPositions[v].x-x0)*slope;
  //       tot_mag += mag;
  //       FField[v.getIndex()]=mag*Normal/3.0;
  //   }


  tot_mag = -1*tot_mag/mesh->nVertices();

  // Now that i have the area average of this magnitude we can add it 
  for(Vertex v : mesh->vertices()){
    FField_correction[v]*=tot_mag;
  }


  // FField_correction = OsmoticPressure(tot_mag);
  // return FField;
  return FField+FField_correction;
  
}







/*
Performs backtracking
Input:

Returns:

void Mem3DG::Smooth_vertices(){






  return;
}



*/

void Mem3DG::Smooth_vertices(){


  size_t N_vert = mesh->nVertices();
  SparseMatrix<double> L(N_vert,N_vert);

  L = geometry->uniformlaplacianMatrix();

  // Ok i have the matrix now what
  Vector<double> xpos = Vector<double>(N_vert);
  Vector<double> ypos = Vector<double>(N_vert);
  Vector<double> zpos =Vector<double>(N_vert);
  Vector3 Pos;
  for(size_t i = 0; i < N_vert; i++){
    Pos = geometry->inputVertexPositions[i];
    xpos[i] = Pos.x;
    ypos[i] = Pos.y;
    zpos[i] = Pos.z;

  }
  xpos = xpos + L*xpos;
  ypos = ypos + L*ypos;
  zpos = zpos + L*zpos;

  for(size_t i = 0; i <N_vert; i++){
    geometry->inputVertexPositions[i] = Vector3({xpos[i],ypos[i],zpos[i]});
  }


  return;
}


/**
 * @brief Performs the backtracking algorithm for the Mem3DG class.
 *
 * This function is responsible for performing the backtracking algorithm in the Mem3DG class. It takes in the force vector, a vector of energy names, and a matrix of energy constants as input. It updates the position of the vertices and beads based on the force, and calculates the energy values for each energy term. It then performs backtracking to find the optimal step size that minimizes the energy.
 *
 * @param Force The force vector applied to the vertices.
 * @param Energies A vector of energy names.
 * @param Energy_constants A matrix of energy constants.
 * @return The optimal step size for the backtracking algorithm.
 */
double Mem3DG::Backtracking(){

  double c1 = 1e-4;
  double rho = 0.5;
  double alpha = 1;
  // alpha = 5e-4;
  double position_Projeection = 0;
  double X_pos;

  // if(system_time>10) std::cout<<"Printing this as system time is " << system_time <<" \n";
  // Sim_handler.Calculate_energies(previousE);
  double previousE = 0;
  Sim_handler->Calculate_energies(&previousE);
  // std::cout<<"Current e is  " << previousE << "\n";
  for(size_t i = 0; i < Sim_handler->Energies.size(); i++) {
    // previousE += Sim_handler.Energy_values[i];
    if(isnan(Sim_handler->Energy_values[i])) std::cout<<"Energy " << Sim_handler->Energies[i] << " is nan\n";
  }

  // std::cout<<"THe previous energy is "<< previousE<<" \n";
  double NewE;
  VertexData<Vector3> initial_pos(*mesh);
  if(recentering){
    if(Field != "None"){
      // std::cout<<"THere is a field\n";
      Vector3 leftmost = Vector3({1e10,1e10,1e10});
      for(Vertex v: mesh->vertices()){
      if(geometry->inputVertexPositions[v].x < leftmost.x) leftmost = geometry->inputVertexPositions[v];
      }
      // I have the leftmost vertex, the position of this vertex should be P 
      leftmost = Vector3({-1.0,0.0,0.0})-leftmost;
      VertexData<Vector3> Displacement(*mesh,leftmost);
      geometry->inputVertexPositions+=Displacement;

    }
    else{
    // Ok so if there is a field the way we normalize is different.
    // 
      Vector3 CoM = geometry->centerOfMass();

      geometry->normalize(Vector3({0.0,0.0,0.0}),false);
    }
  }
  Vector3 CoM = geometry->centerOfMass();
    
  initial_pos = geometry->inputVertexPositions;

  std::vector<Vector3> Bead_init;

  for( size_t i = 0 ; i < Beads.size() ; i++) Bead_init.push_back( Beads[i]->Pos);

  double Projection = 0;
  Vector3 center;  


  // We start the evolution
  // We move the vertices
  // std::cout<<" checking for nan grads\n";
  for(Vertex v : mesh->vertices()){
    if(isnan(Sim_handler->Current_grad[v].x || Sim_handler->Current_grad[v].y || Sim_handler->Current_grad[v].z)) std::cout<<" Is this force nan at vertex but norm2 is " << Sim_handler->Current_grad[v.getIndex()].norm2() <<"\n";
  
  }

  // std::cout<<"doing the stepping \n";
  geometry->inputVertexPositions+=alpha * Sim_handler->Current_grad;
  bool nanflag = false;
  for(Vertex v : mesh->vertices()){
    if(isnan( geometry->inputVertexPositions[v].norm2()))  nanflag = true;
  }

  if(nanflag) std::cout<<"At least one vertex has nan position\n";

  geometry->refreshQuantities();
  center = geometry->centerOfMass();
  Vector3 Vertex_pos;
  
  // We move the beads;
  for(size_t i = 0 ; i < Beads.size(); i++) Beads[i]->Move_bead(alpha,Vector3({0,0,0})); 

  Projection = Sim_handler->Gradient_norms[Sim_handler->Energies.size()];

  size_t bead_count = 0;
  NewE = 0.0;
  // std::cout<<" calculating energies again\n";
  Sim_handler->Calculate_energies(&NewE);
  // std::cout<<"New E is " << NewE << "\n";
  // std::cout<<"The projection is " << Projection << "\n";
  size_t counter = 0;
  
  bool displacement_cond = true;

  // std::cout<<"THe projection is " << Projection <<" \n";
  // std::cout<<"The number of beads is " << Beads.size() << " \n";
  
  if(Projection < 1e-7) {
    small_TS = true;
      std::cout<<"The energy diff is quite small and so is the gradient\n";
      std::cout<<"The energy diff is"<< abs(NewE-previousE)/previousE<<"\n";
      std::cout<<"The projection is"<< Projection<<"\n";
      return -1.0;
  }

  while(true) {
    displacement_cond = true;
    
    for(size_t i = 0 ; i< Beads.size() ; i++) displacement_cond = displacement_cond && Beads[i]->Total_force.norm()*alpha<0.1*Beads[i]->sigma;
    
      if(NewE <= previousE - c1 * alpha * Projection && displacement_cond  && fabs(NewE-previousE)<5e1 ) {
    
        if(fabs(NewE-previousE) > 5e1){
          
      
          std::cout<<"The energies are ";
          for(size_t i = 0; i < Sim_handler->Energies.size(); i++) std::cout<< Sim_handler->Energies[i] << " is " << Sim_handler->Energy_values[i] << " ";
          std::cout<<" \n";
          std::cout<<"The projection is " << Projection <<" \n";

          double Max_projection = 0.0;
          int maxproj_index = 0;

          double maxDisplacement = 0.0;
          std::cout<<"Finding breaking point\n";
          for (Vertex v : mesh->vertices()) {
            if(Sim_handler->Current_grad[v].norm() > Max_projection) {
              Max_projection = Sim_handler->Current_grad[v].norm();
              maxproj_index = v.getIndex();
              }
            double displacement = (geometry->inputVertexPositions[v] - initial_pos[v]).norm();
            if (displacement > maxDisplacement) {
              maxDisplacement = displacement;
            }
          }
          std::cout<<"The max displacement is " << maxDisplacement <<" \n";
          std::cout<<"The value of alpha is " << alpha <<" \n";
          std::cout<<"We will recalculate the energies, lets go back one step for now\n";
          geometry->inputVertexPositions = initial_pos;
          geometry->refreshQuantities();
          mesh->compress();
          // I want something else 
        
        alpha = 0.0;

        // Lets troubleshoot this hehe
        std::cout<<"The previous energy was" << previousE <<" \n";
        std::cout<<"The projection of the bigges vertex is " << Max_projection <<" \n";
        std::cout<<"This vertex is located at " << geometry->inputVertexPositions[maxproj_index] <<" \n";
        std::cout<<"This vertex in init pos is  at " << initial_pos[maxproj_index] <<" \n";
      
        // Lets explore the sorroundings
        Vertex v = mesh->vertex(maxproj_index);
        for(Face f : v.adjacentFaces()){
          std::cout<<"The adjacent faces are " << f.getIndex() <<" \n";
          std::cout<<"With area " << geometry->faceArea(f) <<" \n";
        }
        for(Halfedge he : v.outgoingHalfedges()){
          std::cout<<"The adjacent halfedges are " << he.getIndex() <<" \n";
          std::cout<<"With cotan " << geometry->cotan(he) <<" \n";
        }

        }
      break;
    }
    

    if(std::isnan(NewE)){
      alpha = -1.0;
      break;

    }
    
    alpha *=rho;
    if( (abs((NewE-previousE)/previousE) < 1e-7 && Projection < 0.5) || Projection < 1e-5){
      small_TS = true;
      std::cout<<"The energy diff is quite small and so is the gradient\n";
      std::cout<<"The energy diff is"<< abs(NewE-previousE)/previousE<<"\n";
      std::cout<<"The projection is"<< Projection<<"\n";
      return -1.0;
    }

    if(alpha<1e-10){
      std::cout<<"THe timestep got small so the simulation would end \n";
      std::cout<<"THe timestep is "<< alpha <<" \n";
      std::cout<<"The energy diff is"<< abs(NewE-previousE)<<"\n";
      std::cout<<"THe relative energy diff  is"<<abs((NewE-previousE)/previousE)<<"\n";
      std::cout<<"The projection is"<< Projection<<"\n";
      std::cout<<"The projection is too big "<< (Projection>1.0e8) <<" \n";
      if(Projection>1.0e8){
      // return alpha;
      std::cout<<"The gradient got crazy\n";
      std::cout<<"The projections is "<< Projection<<"\n";
      geometry->inputVertexPositions = initial_pos;
      return -1;
      }
      small_TS = true;

      break;

    }
    else if(small_TS) small_TS = false;
    // std::cout<<"System time is" << system_time <<" \n";
    if(alpha>0) {geometry->inputVertexPositions = initial_pos + alpha*Sim_handler->Current_grad;
    
    for(size_t i = 0; i<Beads.size(); i++){
      Beads[i]->Reset_bead(Bead_init[i]);
      Beads[i]->Move_bead(alpha, Vector3({0,0,0}));
    }
    }
    else {
    for(size_t i = 0; i<Beads.size(); i++){
      Beads[i]->Reset_bead(Bead_init[i]);
      // Beads[i]->Move_bead(alpha, Vector3({0,0,0}));
    }
    geometry->inputVertexPositions = initial_pos;
    }
    
    geometry->refreshQuantities();

    bead_count = 0;

    NewE = 0.0;
    Sim_handler->Calculate_energies(&NewE);
   
  
  }


  nanflag = false;

  for(Vertex v : mesh->vertices()) if(isnan(geometry->inputVertexPositions[v].x+ geometry->inputVertexPositions[v].y  + geometry->inputVertexPositions[v].z )) nanflag = true;  

  if(nanflag) std::cout<< "After backtracking one vertex is nan :( also the value of alpha is"<< alpha << " \n";
  if(alpha<=0.0){ 
  // std::cout<<"Repositioning\n";
  geometry->inputVertexPositions = initial_pos;
  }
  if(recentering) {

    if(alpha<=0.0){
      std::cout<<"NotRecentering after crisis\n";
    }
    else{
    // std::cout<<"rENORMALIZING\n";
    CoM = geometry->centerOfMass();

    if(Field != "None"){
      // std::cout<<"THere is a field\n";
      Vector3 leftmost = Vector3({1e10,1e10,1e10});
      for(Vertex v: mesh->vertices()){
      if(geometry->inputVertexPositions[v].x < leftmost.x) leftmost = geometry->inputVertexPositions[v];
      }
      // I have the leftmost vertex, the position of this vertex should be P 
      leftmost = Vector3({-1.0,0.0,0.0})-leftmost;
      CoM = leftmost;
      VertexData<Vector3> Displacement(*mesh,leftmost);
      geometry->inputVertexPositions+=Displacement;

    }
    else{
    geometry->normalize(Vector3({0.0,0.0,0.0}),false);
    }
  // CoM = geometry->centerOfMass();
    
  for(size_t i = 0; i < Beads.size(); i++){

    // Here i need to move the bead
    // if(Beads[i]->state!= "froze"){
      Beads[i]->Pos -= CoM;
    // }

  }
  // }
  // 
    }
  geometry->refreshQuantities();
  }


  // std::cout<<"The difference in energy is " << fabs(NewE-previousE) <<"(: \n";

  return alpha;
}


/**
 * @brief Performs the backtracking algorithm for the Mem3DG class using a Force given.
 *
 * This function is responsible for performing the backtracking algorithm in the Mem3DG class. It takes in the force vector, a vector of energy names, and a matrix of energy constants as input. It updates the position of the vertices and beads based on the force, and calculates the energy values for each energy term. It then performs backtracking to find the optimal step size that minimizes the energy.
 *
 * @param Force The force vector applied to the vertices.
 * @param Energies A vector of energy names.
 * @param Energy_constants A matrix of energy constants.
 * @return The optimal step size for the backtracking algorithm.
 */
double Mem3DG::Backtracking_BFGS(VertexData<Vector3> Force){

  double c1 = 1e-4;
  double rho = 0.5;
  double alpha = 1;
  // alpha = 5e-4;
  double position_Projeection = 0;
  double X_pos;

  // if(system_time>10) std::cout<<"Printing this as system time is " << system_time <<" \n";
  double previousE = 0;
  Sim_handler->Calculate_energies(&previousE);
  // std::cout<<"Current e is  " << previousE << "\n";

  // std::cout<<"Sim hanlder call\n";
  for(size_t i = 0; i < Sim_handler->Energies.size(); i++) {
    // previousE += Sim_handler.Energy_values[i];
    if(isnan(Sim_handler->Energy_values[i])) std::cout<<"Energy " << Sim_handler->Energies[i] << " is nan\n";
  }

  // std::cout<<"THe previous energy is "<< previousE<<" \n";
  double NewE;
  VertexData<Vector3> initial_pos(*mesh);
  if(recentering){
  geometry->normalize(Vector3({0.0,0.0,0.0}),false);
  }
  Vector3 CoM = geometry->centerOfMass();
    
  initial_pos = geometry->inputVertexPositions;

  std::vector<Vector3> Bead_init;

  // std::cout<<"mOVING BEADS\n"; 
  for( size_t i = 0 ; i < Beads.size() ; i++) Bead_init.push_back( Beads[i]->Pos);

  double Projection = 0;
  Vector3 center;  


  // We start the evolution
  // We move the vertices
  // std::cout<<" checking for nan grads\n";
  // std::cout<<"some nan calcs\n";
  for(Vertex v : mesh->vertices()){
    if(isnan(Force[v].x || Force[v].y || Force[v].z)) std::cout<<" Is this force nan at vertex but norm2 is " << Force[v.getIndex()].norm2() <<"\n";
  
  }

  // std::cout<<"doing the stepping \n";
  geometry->inputVertexPositions+=alpha * Force;
  bool nanflag = false;
  for(Vertex v : mesh->vertices()){
    if(isnan( geometry->inputVertexPositions[v].norm2()))  nanflag = true;
  }

  if(nanflag) std::cout<<"At least one vertex has nan position\n";

  geometry->refreshQuantities();
  center = geometry->centerOfMass();
  Vector3 Vertex_pos;
  
  // We move the beads;
  // std::cout<<"Moving beads\n";
  for(size_t i = 0 ; i < Beads.size(); i++) Beads[i]->Move_bead(alpha,Vector3({0,0,0})); 
  Projection = 0.0;
  // std::cout<<"Calculating projection\n";
  for(Vertex v: mesh->vertices()){
    Projection += Force[v].norm2();
  }

  // Projection = Sim_handler->Gradient_norms[Sim_handler->Energies.size()];

  size_t bead_count = 0;
  NewE = 0.0;
  // std::cout<<" calculating energies again\n";
  Sim_handler->Calculate_energies(&NewE);
  // std::cout<<"New E is " << NewE << "\n";
  // std::cout<<"The projection is " << Projection << "\n";
  size_t counter = 0;
  
  bool displacement_cond = true;

  // std::cout<<"THe projection is " << Projection <<" \n";
  // std::cout<<"The number of beads is " << Beads.size() << " \n";
  
  if(Projection < 1e-7) {
    small_TS = true;
      std::cout<<"The energy diff is quite small and so is the gradient\n";
      std::cout<<"The energy diff is"<< abs(NewE-previousE)/previousE<<"\n";
      std::cout<<"The projection is"<< Projection<<"\n";
      return -1.0;
  }

  while(true) {
    displacement_cond = true;
    
    for(size_t i = 0 ; i< Beads.size() ; i++) displacement_cond = displacement_cond && Beads[i]->Total_force.norm()*alpha<0.1*Beads[i]->sigma;
    
      if(NewE <= previousE - c1 * alpha * Projection && displacement_cond  && fabs(NewE-previousE)<5e1 ) {
    
        if(fabs(NewE-previousE) > 5e1){
          
      
          std::cout<<"The energies are ";
          for(size_t i = 0; i < Sim_handler->Energies.size(); i++) std::cout<< Sim_handler->Energies[i] << " is " << Sim_handler->Energy_values[i] << " ";
          std::cout<<" \n";
          std::cout<<"The projection is " << Projection <<" \n";

          double Max_projection = 0.0;
          int maxproj_index = 0;

          double maxDisplacement = 0.0;
          std::cout<<"Finding breaking point\n";
          for (Vertex v : mesh->vertices()) {
            if(Sim_handler->Current_grad[v].norm() > Max_projection) {
              Max_projection = Sim_handler->Current_grad[v].norm();
              maxproj_index = v.getIndex();
              }
            double displacement = (geometry->inputVertexPositions[v] - initial_pos[v]).norm();
            if (displacement > maxDisplacement) {
              maxDisplacement = displacement;
            }
          }
          std::cout<<"The max displacement is " << maxDisplacement <<" \n";
          std::cout<<"The value of alpha is " << alpha <<" \n";
          std::cout<<"We will recalculate the energies, lets go back one step for now\n";
          geometry->inputVertexPositions = initial_pos;
          geometry->refreshQuantities();
          mesh->compress();
          // I want something else 
        
        alpha = 0.0;

        // Lets troubleshoot this hehe
        std::cout<<"The previous energy was" << previousE <<" \n";
        std::cout<<"The projection of the bigges vertex is " << Max_projection <<" \n";
        std::cout<<"This vertex is located at " << geometry->inputVertexPositions[maxproj_index] <<" \n";
        std::cout<<"This vertex in init pos is  at " << initial_pos[maxproj_index] <<" \n";
      
        // Lets explore the sorroundings
        Vertex v = mesh->vertex(maxproj_index);
        for(Face f : v.adjacentFaces()){
          std::cout<<"The adjacent faces are " << f.getIndex() <<" \n";
          std::cout<<"With area " << geometry->faceArea(f) <<" \n";
        }
        for(Halfedge he : v.outgoingHalfedges()){
          std::cout<<"The adjacent halfedges are " << he.getIndex() <<" \n";
          std::cout<<"With cotan " << geometry->cotan(he) <<" \n";
        }

        }
      break;
    }
    

    if(std::isnan(NewE)){
      alpha = -1.0;
      break;

    }
    
    alpha *=rho;
    if( (abs((NewE-previousE)/previousE) < 1e-7 && Projection < 0.5) || Projection < 1e-5){
      small_TS = true;
      std::cout<<"The energy diff is quite small and so is the gradient\n";
      std::cout<<"The energy diff is"<< abs(NewE-previousE)/previousE<<"\n";
      std::cout<<"The projection is"<< Projection<<"\n";
      return -1.0;
    }

    if(alpha<1e-10){
      std::cout<<"THe timestep got small so the simulation would end \n";
      std::cout<<"THe timestep is "<< alpha <<" \n";
      std::cout<<"The energy diff is"<< abs(NewE-previousE)<<"\n";
      std::cout<<"THe relative energy diff  is"<<abs((NewE-previousE)/previousE)<<"\n";
      std::cout<<"The projection is"<< Projection<<"\n";
      std::cout<<"The projection is too big "<< (Projection>1.0e8) <<" \n";
      if(Projection>1.0e8){
      // return alpha;
      std::cout<<"The gradient got crazy\n";
      std::cout<<"The projections is "<< Projection<<"\n";
      geometry->inputVertexPositions = initial_pos;
      return -1;
      }
      small_TS = true;

      break;

    }
    else if(small_TS) small_TS = false;
    // std::cout<<"System time is" << system_time <<" \n";
    if(alpha>0) {geometry->inputVertexPositions = initial_pos + alpha*Force;
    
    for(size_t i = 0; i<Beads.size(); i++){
      Beads[i]->Reset_bead(Bead_init[i]);
      Beads[i]->Move_bead(alpha, Vector3({0,0,0}));
    }
    }
    else {
    for(size_t i = 0; i<Beads.size(); i++){
      Beads[i]->Reset_bead(Bead_init[i]);
      // Beads[i]->Move_bead(alpha, Vector3({0,0,0}));
    }
    geometry->inputVertexPositions = initial_pos;
    }
    
    geometry->refreshQuantities();

    bead_count = 0;

    NewE = 0.0;
    Sim_handler->Calculate_energies(&NewE);
   
  
  }

  std::cout<<"The energy final is" << NewE <<" \n";

  nanflag = false;

  for(Vertex v : mesh->vertices()) if(isnan(geometry->inputVertexPositions[v].x+ geometry->inputVertexPositions[v].y  + geometry->inputVertexPositions[v].z )) nanflag = true;  

  if(nanflag) std::cout<< "After backtracking one vertex is nan :( also the value of alpha is"<< alpha << " \n";
  if(alpha<=0.0){ 
  // std::cout<<"Repositioning\n";
  geometry->inputVertexPositions = initial_pos;
  }
  if(recentering) {

    if(alpha<=0.0){
      std::cout<<"NotRecentering after crisis\n";
    }
    else{
    // std::cout<<"rENORMALIZING\n";
  CoM = geometry->centerOfMass();
  // if(recentering){
  geometry->normalize(Vector3({0.0,0.0,0.0}),false);
  
  // CoM = geometry->centerOfMass();
    
  for(size_t i = 0; i < Beads.size(); i++){

    // Here i need to move the bead
    // if(Beads[i]->state!= "froze"){
      Beads[i]->Pos -= CoM;
    // }

  }
  
  // }
  // 
    }
  geometry->refreshQuantities();
  }


  // std::cout<<"The difference in energy is " << fabs(NewE-previousE) <<"(: \n";

  return alpha;
}


double Mem3DG::Backtracking(VertexData<Vector3> Force,double P0,double V_bar,double A_bar,double KA,double KB,double H_bar,bool bead, bool pulling) {

// std::cout<<"Backtracking\n";
bool other_pulling=false;
bool two_bead_pulling = true;
double Bead_force_prev=-12;
double Bead_force_new=-15;
double c1=1e-4;
double rho=0.5;
double alpha=1e-3;
double positionProjection = 0;
double X_pos;

// A=geometry->totalArea();
// V=geometry->totalVolume();
// E_Vol = E_Pressure(D_P,V,V_bar);
// E_Sur = E_Surface(KA,A,A_bar);
// E_Ben = E_Bending(H_bar,KB);
// E_Bead = Bead_1.Energy();
double previousE=E_Vol+E_Sur+E_Ben+E_Bead;
double NewE;
VertexData<Vector3> initial_pos(*mesh);
if(recentering){
geometry->normalize(Vector3({0.0,0.0,0.0}),false);
}
initial_pos= geometry->inputVertexPositions;

std::vector<Vector3> Bead_init;

for( size_t i = 0 ; i < Beads.size() ; i++) Bead_init.push_back( Beads[i]->Pos);

// std::cout<<"THe current energy is "<<previousE <<"Is this awful?\n";
// Zeroth iteration
double Projection=0;
Vector3 center; 
    



geometry->inputVertexPositions+=alpha * Force;
geometry->refreshQuantities();

center = geometry->centerOfMass();
Vector3 Vertex_pos;



if(!pulling){
  for(size_t i = 0 ; i < Beads.size(); i++) Beads[i]->Move_bead(alpha,Vector3({0,0,0})); 
  // this->Bead_1.Move_bead(alpha,Vector3({0,0,0}));
}



Total_force = Vector3({0,0,0});
for(Vertex v : mesh->vertices()) {
  Projection+=Force[v.getIndex()].norm2();
  Total_force+=Force[v.getIndex()];
  } 


grad_norm=Projection;


A=geometry->totalArea();
V=geometry->totalVolume();
double D_P = -1*P0*(V-V_bar)/(V_bar*V_bar);
E_Vol=E_Pressure(D_P,V,V_bar);
E_Sur=E_Surface(KA,A,A_bar);
E_Ben=E_Bending(H_bar,KB);
E_Bead = 0;
for(size_t i =0 ; i < Beads.size(); i++) E_Bead+=Beads[i]->Energy();
// E_Bead=Bead_1.Energy();
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

bool displacement_cond = true;

// std::cout<<"starting while\n";
while(true){
  // if(true){
  displacement_cond = true;
  for(size_t i = 0 ; i < Beads.size() ; i++) displacement_cond = displacement_cond && Beads[i]->Total_force.norm()*alpha<0.1;
  // std::cout<<"Displacement cond is"<< displacement_cond <<"\n";
  if( NewE<= previousE - c1 * alpha * Projection && ( displacement_cond ) && abs(NewE-previousE)<10  ) {
    // std::cout<<Bead_1.Total_force.norm()*alpha<<" Displacement of the bead\n";
    break;

    }
  // if(abs(NewE-previousE)>100 && Projection>1e6){
  //   std::cout<<"The energy diff is"<< abs(NewE-previousE)<<"\n";
  //   std::cout<<"THe relative energy diff  is"<<abs((NewE-previousE)/previousE)<<"\n";
  //   std::cout<<"The projection is"<< Projection<<"\n";
  //   std::cout<<"Peak in energy variation, will stop out of safety\n";
  //   alpha=-1;
  //   break;
  // }
  
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
  
    // if(pulling){
      
    //   double current_force=Bead_1.Total_force.norm2();
    //   std::cout<<"The force relative dif is "<< abs((current_force-Bead_1.prev_force)/current_force) <<" \n";
    //   std::cout<<"The bead prev force is"<< Bead_1.prev_force<<"\n";
    //   std::cout<<current_force<<"\n";
    //   if(abs((current_force-Bead_1.prev_force)/current_force)<1e-4){
    //   // std::cout<<"\t \t Reseting bead position\n";
    //   // geometry->normalize(Vector3({0.0,0.0,0.0}),false);
    //   // geometry->refreshQuantities();
    
    //   // X_pos=0.0;
    //   // for(Vertex v : mesh->vertices()){
    //   //   Vertex_pos=geometry->inputVertexPositions[v];
    //   //   if(Vertex_pos.x>X_pos){
    //   //     X_pos=Vertex_pos.x;
    //   //   }  
    //   // }
    //   // std::cout<<" The difference in distance after a step is "<< abs(Bead_1.Pos.x -(X_pos+1.4)) <<" \n";
    //   // if(abs(Bead_1.Pos.x -(X_pos+1.4))<1e-4 ){ //|| Bead_1.Pos.x>X_pos+1.4
    //       // Bead_1.strength = Bead_1.strength+0.1;
    //       // std::cout<<"Increasing strength\n";
    //     // }
    //   // Bead_1.Reset_bead(Vector3(Bead_1.Pos+Vector3({alpha,0,0})));
    //   // Bead_1.Reset_bead(Vector3({X_pos+1.4,0.0,0.0}));
    //   if(X_pos>40.0){
    //     return -1.0;
    //   }
    // // return alpha;
    // break;
    // }
    // }


  if(alpha<1e-10){


    std::cout<<"THe timestep got small so the simulation would end \n";
    std::cout<<"THe timestep is "<< alpha <<" \n";
    std::cout<<"The energy diff is"<< abs(NewE-previousE)<<"\n";
    std::cout<<"THe relative energy diff  is"<<abs((NewE-previousE)/previousE)<<"\n";
    std::cout<<"The projection is"<< Projection<<"\n";
    std::cout<<"The projection is too big "<< (Projection>1.0e10) <<" \n";
    if(Projection>1.0e10){
    // return alpha;
    std::cout<<"The gradient got crazy\n";
    std::cout<<"The projections is "<< Projection<<"\n";
    geometry->inputVertexPositions = initial_pos;
    return -1;
    }

    if(!pulling){
    // if(small_TS==true){
    //   alpha=-1.0;
      
    // }
    if(system_time>999 ){
    small_TS = true;
    std::cout<<"small timestep\n";
    break;
    }
    

    }
    
    // if(Projection>1e10){
    //   std::cout<<"Not moving forward\n";
    //   alpha=0.0;
      // return -1;
    // }
    // continue;
    break;
    

  }
  



  

  geometry->inputVertexPositions = initial_pos+alpha*Force;

  if(pulling && two_bead_pulling){
    // I want to move the beads
      for( size_t i = 0; i < Beads.size(); i++){
        
        if(Beads[i]->state == "default" || Beads[i]->state == "froze") {
          // std::cout<<"TWo here\n";
          Beads[i]->Move_bead(alpha,Vector3({0,0,0}));
          Beads[i]->Reset_bead(Bead_init[i]);
        }
      }



  }
  // geometry->normalize(Vector3({0.0,0.0,0.0}),false);
  
  // center = geometry->centerOfMass();

  if(!pulling){
    for(size_t i = 0 ; i < Beads.size() ; i++){
      Beads[i]->Reset_bead(Bead_init[i]);
      Beads[i]->Move_bead(alpha, Vector3({0,0,0}));
    
    }
  // std::cout<<"pulling is false\n";
  }
  geometry->refreshQuantities();

  
  A=geometry->totalArea();
  V=geometry->totalVolume();
  D_P = -1*P0*(V-V_bar)/(V_bar*V_bar);
  E_Vol=E_Pressure(D_P,V,V_bar);
  E_Sur=E_Surface(KA,A,A_bar);
  E_Ben=E_Bending(H_bar,KB);
  E_Bead = 0;
  for(size_t i =0 ; i < Beads.size(); i++) E_Bead+=Beads[i]->Energy();
  NewE=E_Vol+E_Sur+E_Ben+E_Bead;
  // std::cout<<"THe old energy is "<< previousE <<"\n";
  // std::cout<<"Alpha is "<< alpha<<"and the new energy is"<< NewE << "\n";
  // std::cout<<"The projection is :"<<Projection<<"\n";
  // // std::cout<<"THe energy changed to"<<NewE<<"\n";
  // std::cout<< "Volume E"<<E_Vol <<"Surface E" << E_Sur <<"\n";



}


if(pulling && other_pulling && !two_bead_pulling){
  std::cout<<"pulling is true\n";
      
      // 
      double current_force = Beads[0]->Total_force.norm2();
      // double current_force=Bead_1.Total_force.norm2();
      geometry->normalize(Vector3({0.0,0.0,0.0}),false);
      geometry->refreshQuantities();
      X_pos=0.0;
      for(Vertex v : mesh->vertices()){
        Vertex_pos=geometry->inputVertexPositions[v];
        if(Vertex_pos.x>X_pos){
          X_pos=Vertex_pos.x;
          }  
        }
        // std::cout<<"The system time is"<<system_time<<"\n";
      // if(system_time >2e4){
        std::cout<<"The force difference is "<<abs((current_force-Bead_1.prev_force)/current_force)<<"\n";
      // }
      if(abs((current_force-Beads[0]->prev_force)/current_force)<1e-5 || alpha<1e-10){
        // std::cout<<" The difference in distance after a step is "<< abs(Bead_1.Pos.x -(X_pos+1.4)) <<" \n";
        
        std::cout<<"Some form of steady state was reached\n";
        if(Beads[0]->Pos.x -X_pos<1.4){
          // Just reset the bead
          std::cout<<"\t\t Moving bead\n"; 
          std::cout<<"It is moving "<< abs(Beads[0]->Pos.x -(X_pos+1.4)) <<" \n";
          
          if(abs(Beads[0]->Pos.x-(X_pos+1.4)) <1e-6){
            
            Save_SS = true;
            std::cout<<"\t\t Also increasing interaction strength\n";
            Beads[0]->strength = Beads[0]->strength + 0.1;
          }
          Beads[0]->Reset_bead(Vector3({X_pos+1.4,0,0}));
        
        }
        else{
          // if(abs(Bead_1.Pos.x -(X_pos+1.4))<1e-4 ){
          Save_SS = true;
          std::cout<<"\t\tIncreasing interaction strength because it is going back\n";
          Beads[0]->strength = Beads[0]->strength+0.1;
          }
          
          Beads[0]->prev_E_stationary = E_Bead;

      if(X_pos>40){
        return -1.0;
      }
    //  Bead_1.Reset_bead(Vector3(Bead_1.Pos+Vector3({alpha,0,0})));
    return alpha;
    }
    
    // Bead_1.Reset_bead(Vector3({X_pos+1.4,0,0}));



}

if(pulling && !other_pulling && !two_bead_pulling){
    // std::cout<<"pulling is true\n";
      
      double current_force=Bead_1.Total_force.norm2();
      // std::cout<<"\t \t Reseting bead position\n";
      geometry->normalize(Vector3({0.0,0.0,0.0}),false);
      geometry->refreshQuantities();
      X_pos=0.0;
      for(Vertex v : mesh->vertices()){
        Vertex_pos=geometry->inputVertexPositions[v];
        if(Vertex_pos.x>X_pos){
          X_pos=Vertex_pos.x;
          }  
        }
      if(abs((current_force-Bead_1.prev_force)/current_force)<1e-6 || alpha<1e-10){
        // std::cout<<" The difference in distance after a step is "<< abs(Bead_1.Pos.x -(X_pos+1.4)) <<" \n";
        std::cout<<"It is moving "<< abs(Bead_1.Pos.x -(X_pos+1.4)) <<" \n";
        // if(abs(Bead_1.Pos.x -(X_pos+1.4))<1e-4 ){
          std::cout<<"Increasing strength\n";
          Save_SS = true;
          std::cout<<"The energy difference is "<< abs(E_Bead-Bead_1.prev_E_stationary)<<" \n"; 
          if(abs(E_Bead-Bead_1.prev_E_stationary)<1e-3 && !stop_increasing && Bead_1.strength>10.0){
            stop_increasing = true;
            std::cout<<"\t \t Stopped increasing potential energy\n";
          }
          if(!stop_increasing) Bead_1.strength = Bead_1.strength+0.1;

          Bead_1.prev_E_stationary = E_Bead;

        // }
      Bead_1.Reset_bead(Vector3({X_pos+1.4,0,0}));
      if(X_pos>40){
        return -1.0;
      }
    //  Bead_1.Reset_bead(Vector3(Bead_1.Pos+Vector3({alpha,0,0})));
    return alpha;
    }
    
    Bead_1.Reset_bead(Vector3({X_pos+1.4,0,0}));


    }
  if(pulling && two_bead_pulling){

    for( size_t i = 0 ; i < Beads.size(); i++){        
      if(Beads[i]->state == "manual") {
        Beads[i]->Reset_bead(Bead_init[i]);
        Beads[i]->Move_bead(alpha,Vector3({0,0,0}));
        
      }
      
  }
  }

// std::cout<<"finished while\n";
if(alpha<0.0){
    
    geometry->inputVertexPositions = initial_pos;
}
geometry->normalize(Vector3({0.0,0.0,0.0}),false);
geometry->refreshQuantities();


return alpha;

}


double Mem3DG::Backtracking_field(VertexData<Vector3> Force,double D_P,double V_bar,double A_bar,double KA,double KB, double H_bar) {

// std::cout<<"Backtracking\n";
bool other_pulling=true;
// double Bead_force_prev=-12;
// double Bead_force_new=-15;
double c1=1e-4;
double rho=0.5;
double alpha=1e-3;
double positionProjection = 0;
double X_pos;

// A=geometry->totalArea();
// V=geometry->totalVolume();
// E_Vol = E_Pressure(D_P,V,V_bar);
// E_Sur = E_Surface(KA,A,A_bar);
// E_Ben = E_Bending(H_bar,KB);
// E_Bead = Bead_1.Energy();
double previousE=E_Vol+E_Sur+E_Ben;
double NewE;
VertexData<Vector3> initial_pos(*mesh);

geometry->normalize(Vector3({0.0,0.0,0.0}),false);

initial_pos= geometry->inputVertexPositions;
// Vector3 Bead_init = this->Bead_1.Pos;
// std::cout<<"THe current energy is "<<previousE <<"Is this awful?\n";
// Zeroth iteration
double Projection=0;
Vector3 center; 
    



geometry->inputVertexPositions+=alpha * Force;
geometry->refreshQuantities();

center = geometry->centerOfMass();
Vector3 Vertex_pos;


Total_force = Vector3({0,0,0});
for(Vertex v : mesh->vertices()) {
  Projection+=Force[v.getIndex()].norm2();
  Total_force+=Force[v.getIndex()];
  } 


grad_norm=Projection;


A=geometry->totalArea();
V=geometry->totalVolume();
E_Vol=E_Pressure(D_P,V,V_bar);
E_Sur=E_Surface(KA,A,A_bar);
E_Ben=E_Bending(H_bar,KB);
// E_Bead=Bead_1.Energy();
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

// std::cout<<"starting while\n";
while(true){
  // if(true){
  
  if( NewE<= previousE - c1 * alpha * Projection &&  abs(NewE-previousE)<10  ) {
    // std::cout<<Bead_1.Total_force.norm()*alpha<<" Displacement of the bead\n";
    break;

    }
  // if(abs(NewE-previousE)>100 && Projection>1e6){
  //   std::cout<<"The energy diff is"<< abs(NewE-previousE)<<"\n";
  //   std::cout<<"THe relative energy diff  is"<<abs((NewE-previousE)/previousE)<<"\n";
  //   std::cout<<"The projection is"<< Projection<<"\n";
  //   std::cout<<"Peak in energy variation, will stop out of safety\n";
  //   alpha=-1;
  //   break;
  // }
  
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
  
    // if(pulling){
      
    //   double current_force=Bead_1.Total_force.norm2();
    //   std::cout<<"The force relative dif is "<< abs((current_force-Bead_1.prev_force)/current_force) <<" \n";
    //   std::cout<<"The bead prev force is"<< Bead_1.prev_force<<"\n";
    //   std::cout<<current_force<<"\n";
    //   if(abs((current_force-Bead_1.prev_force)/current_force)<1e-4){
    //   // std::cout<<"\t \t Reseting bead position\n";
    //   // geometry->normalize(Vector3({0.0,0.0,0.0}),false);
    //   // geometry->refreshQuantities();
    
    //   // X_pos=0.0;
    //   // for(Vertex v : mesh->vertices()){
    //   //   Vertex_pos=geometry->inputVertexPositions[v];
    //   //   if(Vertex_pos.x>X_pos){
    //   //     X_pos=Vertex_pos.x;
    //   //   }  
    //   // }
    //   // std::cout<<" The difference in distance after a step is "<< abs(Bead_1.Pos.x -(X_pos+1.4)) <<" \n";
    //   // if(abs(Bead_1.Pos.x -(X_pos+1.4))<1e-4 ){ //|| Bead_1.Pos.x>X_pos+1.4
    //       // Bead_1.strength = Bead_1.strength+0.1;
    //       // std::cout<<"Increasing strength\n";
    //     // }
    //   // Bead_1.Reset_bead(Vector3(Bead_1.Pos+Vector3({alpha,0,0})));
    //   // Bead_1.Reset_bead(Vector3({X_pos+1.4,0.0,0.0}));
    //   if(X_pos>40.0){
    //     return -1.0;
    //   }
    // // return alpha;
    // break;
    // }
    // }


  if(alpha<1e-10){


    std::cout<<"THe timestep got small so the simulation would end \n";
    std::cout<<"THe timestep is "<< alpha <<" \n";
    std::cout<<"The energy diff is"<< abs(NewE-previousE)<<"\n";
    std::cout<<"THe relative energy diff  is"<<abs((NewE-previousE)/previousE)<<"\n";
    std::cout<<"The projection is"<< Projection<<"\n";
    

    // if(Projection>1e10){
    //   std::cout<<"Not moving forward\n";
    //   alpha=0.0;
      // return -1;
    // }
    // continue;
    break;
    

  }

  // geometry->inputVertexPositions = initial_pos+alpha*Force;
  // geometry->normalize(Vector3({0.0,0.0,0.0}),false);
  // geometry->refreshQuantities();
  // center = geometry->centerOfMass();

  

  
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





// std::cout<<"finished while\n";
if(alpha<0.0){
    
    geometry->inputVertexPositions = initial_pos;
}
// geometry->normalize(Vector3({0.0,0.0,0.0}),false);
// geometry->refreshQuantities();


return alpha;

}

// THis is for the membrane flow

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

// double previousE=E_Vol+E_Sur+E_Ben;
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
// NewE=E_Sur+E_Ben;
// if(std::isnan(E_Vol)){
//   std::cout<<"E vol is nan\n";
// }
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
  // if(std::isnan(E_Vol)){
  // std::cout<<"E vol is nan\n";
  //   }
if(std::isnan(E_Sur)){
  std::cout<<"E sur is nan\n";
    }
if(std::isnan(E_Ben)){
  std::cout<<"E ben is nan\n";
    }

  

  alpha*=rho;
  if(alpha<1e-8){
    // std::cout<<"THe timestep got small\n";
    if(system_time<2*Area_evol_steps){
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
  // NewE=E_Sur+E_Ben;
  NewE=E_Vol+E_Sur+E_Ben;
  
  if(std::isnan(NewE)){
    std::cout<<"The energy got Nan\n";
    
    alpha=-1.0;
    break;
  }




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


// THis is backtracking for volume preserving mean curvature flow
double Mem3DG::Backtracking(VertexData<Vector3> Force,double D_P,double V_bar,double KA) {
double c1=1e-4;
double rho=0.7;
double alpha=1e-3;
double positionProjection = 0;
double A=geometry->totalArea();
double V=geometry->totalVolume();
double E_Vol=E_Pressure(D_P,V,V_bar);
double E_Sur=A*KA;

// double previousE=E_Vol+E_Sur+E_Ben;
double previousE=E_Vol+E_Sur;
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
E_Sur=A*KA;
NewE=E_Vol+E_Sur;
// NewE=E_Sur+E_Ben;
// if(std::isnan(E_Vol)){
//   std::cout<<"E vol is nan\n";
// }
if(std::isnan(E_Sur)){
  std::cout<<"E sur is nan\n";
}



size_t counter=0;
while(true){
  // if(true){
  // std::cout<<"THe new energy is "<<NewE <<"\n";
  if( NewE<= previousE - c1*alpha*Projection ) {
    break;

    }
  // if(std::isnan(E_Vol)){
  // std::cout<<"E vol is nan\n";
  //   }
if(std::isnan(E_Sur)){
  std::cout<<"E sur is nan\n";
    }

  if(std::isnan(NewE)){
    std::cout<<"The energy got Nan\n";
    
    alpha=-1.0;
    break;
  }


  alpha*=rho;
  if(alpha<1e-8){
    // std::cout<<"THe timestep got small\n";
    if(system_time<2*Area_evol_steps){
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
  E_Sur=A*KA;
  NewE=E_Vol+E_Sur;  
  

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


// void Mem3DG::Get_Energies(std::vector<std::string>Energies, std::vector<std::vector<double>> Energy_constants, double* NewE){

//   double V_bar;

//   *NewE = 0.0;
//   // *NewE = 1.0;

//   int bead_count = 0;
//   for(size_t i = 0; i < Energies.size(); i++){
//   // We will do the calculations here
//     if(Energies[i]=="Volume_constraint"){
//       double V = geometry->totalVolume();
//       V_bar = Energy_constants[i][1];
//       // double D_?P = -1*Energy_constants[i][0]*(V-V_bar)/(V_bar*V_bar);

//       // We calculate the energy
//       Energy_vals[i] = E_Volume_constraint(Energy_constants[i][0],V,V_bar);
//       *NewE += Energy_vals[i];
//       continue;
//     }
//     if(Energies[i]=="Area_constraint"){
      
//       double A = geometry->totalArea();
//       double A_bar = Energy_constants[i][1];
//       double KA = Energy_constants[i][0];

//       Energy_vals[i] = E_Area_constraint(KA,A,A_bar);
//       *NewE += Energy_vals[i];
//       continue;
//     }
//     if(Energies[i]=="Bending" || Energies[i] == "H1_Bending" || Energies[i]=="H2_Bending"){
//       // std::cout<<"Calculating BENERGY 1\n";
//       double KB = Energy_constants[i][0];
//       double H0 = 0.0;

//       Energy_vals[i] = E_Bending(H0, KB);
//       *NewE += Energy_vals[i];
//       continue;

//     }

//     if(Energies[i]=="Bead" || Energies[i]=="H1_Bead" || Energies[i]=="H2_Bead"){

//       // The energy
//       Energy_vals[i] = Beads[bead_count]->Energy(); 

//       bead_count +=1;
//       *NewE += Energy_vals[i];
//       continue;


//     }
//     if(Energies[i]=="Surface_tension"|| Energies[i] == "H1_Surface_tension" || Energies[i] == "H2_Surface_tension"){
//       // What if the energy is surface tension 
//       double sigma = Energy_constants[i][0];
//       // I need to add an energy term here;
//       // std::cout<<"Calculating old energy\n";
//       A  = geometry->totalArea();
//       Energy_vals[i] = sigma*A; 
//       *NewE += Energy_vals[i];
//       continue;

//     }

//   }

//   return;
// }


double Mem3DG::integrate(std::ofstream& Sim_data , double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data){

  auto start = chrono::steady_clock::now();
  auto end = chrono::steady_clock::now();
  auto construction_start = chrono::steady_clock::now();
  auto construction_end = chrono::steady_clock::now();
  auto solve_start = chrono::steady_clock::now();
  auto solve_end = chrono::steady_clock::now();

  
  // std::cout<<"Got to integrating\n";

  double time_construct = 0;
  double time_solve = 0;
  double time_compute = 0;
  double time_gradients = 0;
  double time_backtracking = 0;


  size_t bead_count = 0;
  // std::cout<<"Bead data\n";
  if(Bead_data_filenames.size()!=0 && Save_output_data){
    std::ofstream Bead_data;
    for(size_t i = 0; i < Beads.size(); i++){
      
      Bead_data = std::ofstream(Bead_data_filenames[i],std::ios_base::app);
      Bead_data<<Beads[i]->Pos.x <<" "<< Beads[i]->Pos.y << " "<< Beads[i]->Pos.z<< " "<< Beads[i]->Total_force.x <<" "<< Beads[i]->Total_force.y << " "<< Beads[i]->Total_force.z <<" \n";
      // std::cout<<Bead_1.Pos.x << " "<< Bead_1.Pos.y << " "<< Bead_1.Pos.z<<" \n";
      // std::cout<<"The total force is "<<Bead_1.Total_force <<"\n";
      Bead_data.close();
      }

  }


  // I need to get the force and their gradients


  // std::vector<double> Gradient_norms;
  // double grad_value;
  // VertexData<Vector3> Force(*mesh,Vector3({0.0,0.0,0.0}));
  // VertexData<Vector3> Force_temp(*mesh);
  // if(Energy_vals.size()==0) Energy_vals.resize(Energies.size());
  
  // double A_bar = 4.0*3.1415926535;
  // double V_bar = 4.0*3.14115926535/3.0;
  // double V;
  // double A = 0;

  // Lets calculate the areas of every vertex

  // std::cout<<"Barycentric area\n";
  Vector<double> Barycentric_area(mesh->nVertices());

  for(size_t index = 0; index < mesh->nVertices(); index++) {
    Barycentric_area[index] = geometry->barycentricDualArea(mesh->vertex(index));
    A += Barycentric_area[index];
    Barycentric_area[index] = 1.0;
  }


  bool H2_grad = false;
  bool H1_grad = false;
  size_t N_vert = mesh->nVertices();
  int N_constraints = 4;
  SparseMatrix<double> H2_mat;//(N_vert*3+N_constraints,N_vert*3+N_constraints);
  SparseMatrix<double> H1_mat;
  Eigen::SparseLU<SparseMatrix<double>> solverH2;
  Eigen::SparseLU<SparseMatrix<double>> solverH1;
  // int N_vert = mesh->nVertices();

  // We will calculate the gradients here
  start = chrono::steady_clock::now();
  // std::cout<<"Calling simulation handler for gradiebts\n";
  // std::cout<<"Calling sim handler\n";
  Sim_handler->Calculate_gradient();

  // std::cout<<"Got the gradient\n";
  // std::cout<<"Sucesfully calculated\n";
  // Sim_handler->Do_nothing();
  // Sim_handler->mesh->nVertices();
  // Sim_handler->Do_nothing();
  
  // Sim_handler->Add_Bead(Beads[0]);

  end = chrono::steady_clock::now();

  time_gradients = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  
  
  // THis is the time it takes for the gradients 



  // std::cout<<" \n";
  //  std::cout<<"The energy vals are ";
  // for(size_t i = 0; i < Energy_vals.size(); i++){
  //   std::cout<< Energy_vals[i] <<" ";

  // }
  // std::cout<<" \n";
  // Ok now i have the force
  double alpha = 1e-3;
  double backtrackstep;
  Total_force = Vector3({0.0, 0.0, 0.0});
  double Grad_tot_norm  = 0;
 
  V = geometry->totalVolume();
  double r_eff = 0;


  // F_dist.close();
  // std::cout<<"Moving to backtracking\n";
  
  start = chrono::steady_clock::now();
  backtrackstep = Backtracking();
  end = chrono::steady_clock::now();
  time_backtracking = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

  // OK NOW I HAVE TO save the timings

  if(Field !="None"){
    VertexData<Vector3> Field(*mesh,Vector3({0.0,0.0,0.0}));
    Field = Linear_force_field(Field_vals[0],Field_vals[1]);
    geometry->inputVertexPositions+=backtrackstep*Field;
  }
  // Ok so we added the force field

  std::ofstream Timings = std::ofstream(Bead_data_filenames[Beads.size()],std::ios_base::app);
  Timings << time_gradients <<" "<< time_backtracking <<" "<< time_construct <<" "<< time_compute <<" "<< time_solve <<" \n";
  Timings.close();
  // std::cout<<"The backtrackstep is " << backtrackstep << " \n";
  // After backtracking i have to save.

  // 
  if(Save_output_data || backtrackstep <0 ){
  double tot_E=0;
  Sim_data << time <<" "<< V<<" " << A<<" ";
  for(size_t i = 0; i < Sim_handler->Energies.size(); i++){
    // std::cout<<"Printing " << Energies[i] << " ";

    Sim_data << Sim_handler->Energy_values[i] << " ";
    // std::cout<<"the val is " << Energy_vals[i] << " ";
    tot_E += Sim_handler->Energy_values[i];
  }
  // std::cout<<" \n";
  Sim_data<< tot_E <<" ";
  for(size_t i = 0; i < Sim_handler->Gradient_norms.size(); i++){
    Sim_data << Sim_handler->Gradient_norms[i]<< " ";
  }
  // Sim_data << Grad_tot_norm << " ";
  Sim_data<< backtrackstep<<" \n";
  
    }

  




  return backtrackstep;
}

Eigen::MatrixXd Mem3DG::integrate_BFGS(std::ofstream& Sim_data , double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data, Eigen::MatrixXd Hessian){

  // We need to store the Hessian somewhere, maybe it should be received as a pointer
  auto start = chrono::steady_clock::now();
  auto end = chrono::steady_clock::now();
  auto construction_start = chrono::steady_clock::now();
  auto construction_end = chrono::steady_clock::now();
  auto solve_start = chrono::steady_clock::now();
  auto solve_end = chrono::steady_clock::now();

  

  double time_construct = 0;
  double time_solve = 0;
  double time_compute = 0;
  double time_gradients = 0;
  double time_backtracking = 0;


  size_t bead_count = 0;
  // std::cout<<"Bead data\n";
  if(Bead_data_filenames.size()!=0 && Save_output_data){
    std::ofstream Bead_data;
    for(size_t i = 0; i < Beads.size(); i++){
      
      Bead_data = std::ofstream(Bead_data_filenames[i],std::ios_base::app);
      Bead_data<<Beads[i]->Pos.x <<" "<< Beads[i]->Pos.y << " "<< Beads[i]->Pos.z<< " "<< Beads[i]->Total_force.x <<" "<< Beads[i]->Total_force.y << " "<< Beads[i]->Total_force.z <<" \n";
      // std::cout<<Bead_1.Pos.x << " "<< Bead_1.Pos.y << " "<< Bead_1.Pos.z<<" \n";
      // std::cout<<"The total force is "<<Bead_1.Total_force <<"\n";
      Bead_data.close();
      }

  }
  // std::cout<<"BEAD DATA DONE\n";



  // I need to get the force and their gradients


  // std::vector<double> Gradient_norms;
  // double grad_value;
  // VertexData<Vector3> Force(*mesh,Vector3({0.0,0.0,0.0}));
  // VertexData<Vector3> Force_temp(*mesh);
  // if(Energy_vals.size()==0) Energy_vals.resize(Energies.size());
  
  // double A_bar = 4.0*3.1415926535;
  // double V_bar = 4.0*3.14115926535/3.0;
  // double V;
  // double A = 0;

  start = chrono::steady_clock::now();
  // std::cout<<"Calling simulation handler for gradiebts\n";

  // backtrackstep = Sim_handler->Backtracking(BFGS);s
  // I need to check here that 
  // std::cout<<"Calculating gradient\n";


  Sim_handler->Calculate_gradient();
  // std::cout<<"Gradient calculated\n";

  
  end = chrono::steady_clock::now();

  time_gradients = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  
  
  // THis is the time it takes for the gradients 



  // std::cout<<" \n";
  //  std::cout<<"The energy vals are ";
  // for(size_t i = 0; i < Energy_vals.size(); i++){
  //   std::cout<< Energy_vals[i] <<" ";

  // }
  // std::cout<<" \n";
  // Ok now i have the force
  double alpha = 1e-3;
  double backtrackstep;
  Total_force = Vector3({0.0, 0.0, 0.0});
  double Grad_tot_norm  = 0;
 
  V = geometry->totalVolume();
  double r_eff = 0;


  // F_dist.close();
  // std::cout<<"Moving to backtracking\n";
  
  start = chrono::steady_clock::now();
  
  // We want to know how this is transformed 

  // std::cout<<"Started converting stuff\n";
  Eigen::VectorXd Grad_vec = Eigen::VectorXd::Zero(3*mesh->nVertices());

  for(Vertex v : mesh->vertices()){
    Grad_vec[v.getIndex()] = Sim_handler->Current_grad[v].x;
    Grad_vec[v.getIndex()+mesh->nVertices()] = Sim_handler->Current_grad[v].y;
    Grad_vec[v.getIndex()+2*mesh->nVertices()] = Sim_handler->Current_grad[v].z;
  }

  // std::cout<<"Finished loading grad_vec\n";

  // std::cout<<"The size of gravec is " << Grad_vec.size() <<" \n";

  // std::cout<<"The gradvec at some position is " << Grad_vec[10] << " \n";
  // std::cout<<"The t gradvec at same position is " << Grad_vec.transpose()[10] <<" \n";
  Grad_vec = Hessian * Grad_vec;

  // std::cout<<"The gradvec at some position is " << Grad_vec[10] << " \n";
  // std::cout<<"The t gradvec at same position is " << Grad_vec.transpose()[10] <<" \n";
  
  // std::cout<<"Works?\n";
  // Eigen::RowVectorXd T_Grad_vec1 = Grad_vec.transpose();

  // std::cout<<"The size of gravecafter is " << Grad_vec.size() <<" \n";


  // std::cout<<"Vertex product\n";
  VertexData<Vector3> Force(*mesh,Vector3({0.0,0.0,0.0}));
  for(Vertex v : mesh->vertices()){
    // Vector3 Some_Force = Vector3({Grad_vec[v.getIndex()],Grad_vec[v.getIndex()+mesh->nVertices()],Grad_vec[v.getIndex()+2*mesh->nVertices()]});
    // std::cout<<"Force norm2 is " << Some_Force.norm2() << " and currentgrad norm2 is " << Sim_handler->Current_grad[v].norm2() << " \n";
    Force[v].x = Grad_vec[v.getIndex()];
    Force[v].y = Grad_vec[v.getIndex()+mesh->nVertices()];
    Force[v].z = Grad_vec[v.getIndex()+2*mesh->nVertices()];
  }
  // std::cout<<"Force converted\n now backtracking \n";

  backtrackstep = Backtracking_BFGS( Force);

  // std::cout<<"Backtrackstep is " << backtrackstep <<" \n";


  // Here we just moved the positions, lets see what we need to do
  // std::cout<<"Backtracking down boots\n New gradient\n";
  Sim_handler->Calculate_gradient();

  
  // std::cout<<"done with it, now magic quantities\n";
  VertexData<Vector3> yk = Sim_handler->Current_grad-Sim_handler->Previous_grad;

  // std::cout<<"YK DONE\n";
  // Deleting the scalar product
  Grad_vec = backtrackstep* Grad_vec; //So this is sk
  
  // std::cout<<"The size of grade vec after is also now " << Grad_vec.size() <<" \n";
  // std::cout<<"THis works\n";

  Eigen::VectorXd yk_vec = Eigen::VectorXd::Zero(3*mesh->nVertices());
  for(Vertex v : mesh->vertices()){
    yk_vec[v.getIndex()] = yk[v].x;
    yk_vec[v.getIndex()+mesh->nVertices()] = yk[v].y;
    yk_vec[v.getIndex()+2*mesh->nVertices()] = yk[v].z;
  }
  // std::cout<<"Yk vec at 10 is" << yk_vec[10] << " \n";

  // std::cout<<"YK VEC DONE\n";
  // I have yk ready, sk and Hk

  // Lets see whos size is wrong
  // Eigen::DenseMatri
  // std::cout<<"THe size of Grad_vec is " << Grad_vec.size() << " \n";
  // std::cout<<"THe size of ykvec is " << yk_vec.size() << " \n";
  // std::cout<<"The size of Hessian is " << Hessian.rows() << "  and " << Hessian.cols() <<"\n";

  // std::cout<<"THe size of "<< (Grad_vec* Grad_vec.transpose()).size() << " \n";
  // std::cout<<"THe size of "<< (yk_vec.transpose()*Grad_vec).size() << " \n";
  

  // std::cout<<"The rows are " << (Grad_vec * Grad_vec.transpose()).coeff(0,0) << " and the columns are " << (Grad_vec * Grad_vec.transpose()).cols() << "\n";
  // He ssian = Hessian + (Grad_vec.transpose()* yk_vec + yk_vec.transpose()*Hessian * yk_vec)*(Grad_vec*Grad_vec.transpose())/(Grad_vec.transpose()*yk_vec*yk_vec.transpose()*Grad_vec) - (Hessian * yk_vec*Grad_vec.transpose() + Grad_vec*yk_vec.transpose()*Hessian)/(yk_vec.transpose()*Grad_vec); 


  // Eigen::VectorXd Simple_vec(4);
  // Simple_vec << 1,2,3,4;
  // Eigen::VectorXd other_simple_vec(4);
  // other_simple_vec << 1,2,3,4;

  // std::cout<<"Simple vec is " << Simple_vec << " \n";
  // std::cout<<"First they are different\n";
  // std::cout<<Simple_vec * other_simple_vec.transpose() << " \n";
  // std::cout<<"Second they are the same";
  // std::cout<< Simple_vec * Simple_vec.transpose() << " \n";

  // Eigen::RowVectorXd trans_simple_vec = Simple_vec.transpose();
  // std::cout<<"Simple vec transposed easy " << trans_simple_vec <<" \n";
  Eigen::RowVectorXd T_Grad_vec = Grad_vec.transpose();
  // std::cout<<"This row vec exists\n";
  // std::cout<<T_Grad_vec << " \n";

  // DenseMatrix<double> Hessian_temp = (Grad_vec*Grad_vec.transpose());
  // std::cout<<"Lets output the vector " << Grad_vec <<" \n";
  
  // std::cout<< Grad_vec.transpose().size() << " \n";
  // std::cout<<"The gradvec size is " << Grad_vec.transpose().size() << " \n";  

  // Eigen::MatrixXd Hessian_temp = Grad_vec * T_Grad_vec1;


  Hessian = Hessian + (Grad_vec.dot(yk_vec) + yk_vec.dot(Hessian * yk_vec) )*(Grad_vec*T_Grad_vec)/(Grad_vec.dot(yk_vec) *yk_vec.dot(Grad_vec) ) - (Hessian * yk_vec*Grad_vec.transpose() + Grad_vec*yk_vec.transpose()*Hessian)/(yk_vec.dot(Grad_vec));
  // std::cout<<"Big op done\n";


  end = chrono::steady_clock::now();
  time_backtracking = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

  // OK NOW I HAVE TO save the timings

  std::ofstream Timings = std::ofstream(Bead_data_filenames[Beads.size()],std::ios_base::app);
  Timings << time_gradients <<" "<< time_backtracking <<" "<< time_construct <<" "<< time_compute <<" "<< time_solve <<" \n";
  Timings.close();
  // std::cout<<"The backtrackstep is " << backtrackstep << " \n";
  // After backtracking i have to save.

  // 
  if(Save_output_data || backtrackstep <0 ){
  double tot_E=0;
  Sim_data << time <<" "<< V<<" " << A<<" ";
  for(size_t i = 0; i < Sim_handler->Energies.size(); i++){
    // std::cout<<"Printing " << Energies[i] << " ";

    Sim_data << Sim_handler->Energy_values[i] << " ";
    // std::cout<<"the val is " << Energy_vals[i] << " ";
    tot_E += Sim_handler->Energy_values[i];
  }
  // std::cout<<" \n";
  Sim_data<< tot_E <<" ";
  for(size_t i = 0; i < Sim_handler->Gradient_norms.size(); i++){
    Sim_data << Sim_handler->Gradient_norms[i]<< " ";
  }
  // Sim_data << Grad_tot_norm << " ";
  Sim_data<< backtrackstep<<" \n";
  
    }

  




  return Hessian;
}


double Mem3DG::integrate_implicit(std::vector<std::string> Energies,  std::vector<std::vector<double>> Energy_constants, std::ofstream& Sim_data , double time, std::vector<std::string>Bead_data_filenames, bool Save_output_data){


  // Okok what now

  

  return 0.0;
}




/*
 * Performs integration with the bead
 *
 * Input: The timestep <h>.
 * Returns:
 */
double Mem3DG::integrate(double h, double V_bar, double nu, double c0,double P0,double KA,double KB, double Kd,std::ofstream& Sim_data, double time, bool bead,std::vector<std::string> Bead_data_filenames,bool Save_output_data,bool pulling) {
//, Beads Bead_1 
    std::ofstream Bead_data;
    // std::cout<<"1_2\n";
    Bead *Active_bead;
    if(bead ){
      // std::cout<<"This is being called\n";
      // std::cout<<"There are "<< Beads.size() <<" beads\n";
      for(size_t i = 0; i < Beads.size(); i++){
      
      Bead_data = std::ofstream(Bead_data_filenames[i],std::ios_base::app);
      Bead_data<<Beads[i]->Pos.x <<" "<< Beads[i]->Pos.y << " "<< Beads[i]->Pos.z<< " "<< Beads[i]->Total_force.x <<" "<< Beads[i]->Total_force.y << " "<< Beads[i]->Total_force.z <<" \n";
      // std::cout<<Bead_1.Pos.x << " "<< Bead_1.Pos.y << " "<< Bead_1.Pos.z<<" \n";
      // std::cout<<"The total force is "<<Bead_1.Total_force <<"\n";
      Bead_data.close();
      }
    }

    // std::cout<<"2_2\n";
    // Vector<double> Total_force=buildFlowOperator(h,V_bar,nu,c0,P0,KA,KB,Kd);

    
    VertexData<Vector3> Force(*mesh);
    
    Force=buildFlowOperator(h,V_bar,nu,c0,P0,KA,KB,Kd);//+Bead_1.Gradient();


    // std::cout<<"3_2\n";
    // std::cout<<"\t \tThe size of Beads is"<< Beads.size()<<"\n";
    VertexData<Vector3> Bead_force;
    for(size_t i = 0; i < Beads.size(); i++){
      // std::cout<<"Iterating over bead number "<<i+1<<"\n";
      Bead_force = Beads[i]->Gradient();
      Beads[i]->Bead_interactions();
      Force+=Bead_force;
    }
    // VertexData<Vector3> Bead_force = Bead_1.Gradient();

    // VertexData<Vector3> Bead_force=Project_force(Bead_1.Gradient());
    
    // std::cout<<"4_2\n";
    // Force+=Bead_force;


    // This force is the grdient basically
    double alpha=h;

    // I have the forces for almost everything i just need the bead.

    size_t vindex;
    size_t Nvert=mesh->nVertices();


    // I want to print the Volume, the area, VOl_E, Area_E, Bending_E 
    V = geometry->totalVolume();
    A = geometry->totalArea();
    double A_bar=4*PI*pow(3*V_bar/(4*PI*nu),2.0/3.0);
    double H_bar=sqrt(4*PI/A_bar)*c0/2.0; //Coment this with another comment
    double D_P=-1*P0*(V-V_bar)/(V_bar*V_bar);
    double lambda=KA*(A-A_bar )/(A_bar*A_bar);    



    E_Vol = E_Pressure(D_P,V,V_bar);
    E_Sur = E_Surface(KA,A,A_bar);
    E_Ben = E_Bending(H_bar,KB);
    E_Bead= 0.0;

    for( size_t i = 0; i < Beads.size(); i++){

    E_Bead += Beads[i]->Energy();
    }

    double backtrackstep;

    // std::cout<<" The position of vertex 120 is "<<geometry->inputVertexPositions[120].x<<" "<<geometry->inputVertexPositions[120].y<< " "<<geometry->inputVertexPositions[120].z<<" \n"; 

    // if(system_time< 50){
    //   backtrackstep=h;
    // }
    // else{
      // std::cout<<"Backtracking is being called\n";
    Total_force = Vector3({0,0,0});
    backtrackstep=Backtracking(Force,P0,V_bar,A_bar,KA,KB,H_bar,bead,pulling);
    // Total_force*=backtrackstep;
    // I need to find the stop increasing bt
    // }
    // std::cout<<" The position of vertex 120 is "<<geometry->inputVertexPositions[120].x<<" "<<geometry->inputVertexPositions[120].y<< " "<<geometry->inputVertexPositions[120].z<<" \n";
    
        
    // geometry->inputVertexPositions=geometry->inputVertexPositions+backtrackstep*Force;

    // Bead_1.Move_bead(backtrackstep,center);

    if((Save_output_data ) || backtrackstep<0  ){
    Sim_data << V_bar<<" "<< A_bar<<" "<< time <<" "<< V<<" " << A<<" " << E_Vol << " " << E_Sur << " " << E_Ben <<" " << E_Bead << " "<< grad_norm<<" " << backtrackstep << " \n";
    }
    system_time+=1;
    
    
  
  return backtrackstep;
} 

   

double Mem3DG::integrate_field(double h, double V_bar, double nu, double P0,double KA,double KB, double slope,double x0 ,std::ofstream& Sim_data, double time,bool Save) {
//, Beads Bead_1 


    // std::cout<<"2_2\n";
    // Vector<double> Total_force=buildFlowOperator(h,V_bar,nu,c0,P0,KA,KB,Kd);
    VertexData<Vector3> Force(*mesh);
    VertexData<Vector3> Field(*mesh);
    double c0 = 0.0;
    double Kd = 0.0;
    Force=buildFlowOperator(V_bar,P0,KA,KB,h);//+Bead_1.Gradient();
    


    // This force is the grdient basically
    double alpha=h;

    // I have the forces for almost everything i just need the bead.

    size_t vindex;
    size_t Nvert=mesh->nVertices();


    // I want to print the Volume, the area, VOl_E, Area_E, Bending_E 
    V = geometry->totalVolume();
    A = geometry->totalArea();
    double A_bar=4*PI*pow(3*V_bar/(4*PI*nu),2.0/3.0);
    double H_bar=sqrt(4*PI/A_bar)*c0/2.0; //Coment this with another comment
    double D_P=-1*P0*(V-V_bar)/(V_bar*V_bar);
    double lambda=KA*(A-A_bar )/(A_bar*A_bar);    



    E_Vol = E_Pressure(D_P,V,V_bar);
    E_Sur = E_Surface(KA,A,A_bar);
    E_Ben = E_Bending(H_bar,KB);
    E_Bead = Bead_1.Energy();
    
    // double backtrackstep;

    // std::cout<<" The position of vertex 120 is "<<geometry->inputVertexPositions[120].x<<" "<<geometry->inputVertexPositions[120].y<< " "<<geometry->inputVertexPositions[120].z<<" \n"; 

    // if(system_time< 50){
    //   backtrackstep=h;
    // }
    // else{
      // std::cout<<"Backtracking is being called\n";
    Total_force = Vector3({0,0,0});
    // backtrackstep=Backtracking(Force,D_P,V_bar,A_bar,KA,KB,H_bar,bead,pulling);
    // Total_force*=backtrackstep;
    // I need to find the stop increasing bt
    // }
    // std::cout<<" The position of vertex 120 is "<<geometry->inputVertexPositions[120].x<<" "<<geometry->inputVertexPositions[120].y<< " "<<geometry->inputVertexPositions[120].z<<" \n";
    // double backtrackstep = Backtracking_field(Force,D_P,V_bar,A_bar,KA,KB,H_bar);    
    double backtrackstep = 1e-4; //lowest 5e-7
    // I need to calculate the force field


    Field = Linear_force_field(x0,slope); 
    // I want to make it volume preserving



    Force += Field;
    geometry->inputVertexPositions=geometry->inputVertexPositions+backtrackstep*Force;
    
    
    
    
    // I want to find the leftmost particle
    
    Vector3 Leftmost({0,0,0});
    Vector3 Rightmost({0,0,0});
    
    int rightmost;
    int leftmost;
    for(Vertex v: mesh->vertices()){
      if(geometry->inputVertexPositions[v].x<Leftmost.x){
        Leftmost = geometry->inputVertexPositions[v];
        leftmost = v.getIndex();
      }
      if(geometry->inputVertexPositions[v].x>Rightmost.x){
        Rightmost = geometry->inputVertexPositions[v];
        rightmost = v.getIndex();
      }
    
    }
    Leftmost-= Vector3({-5,0,0});
    VertexData<Vector3> Displacement(*mesh,Leftmost);

    geometry->inputVertexPositions-=Displacement;

    // std::cout<<"The force at the leftmost is"<< Field[leftmost] <<" with magnitude "<< Field[leftmost].norm()<<"\n";
    // std::cout<<"The force at the leftmost is"<< Field[rightmost] <<" with magnitude "<< Field[rightmost].norm()<<"\n";

    // std::cout<<" The difference in force in theory " << slope*(geometry->inputVertexPositions[rightmost].x - geometry->inputVertexPositions[leftmost].x)<< " and the obtained is " << Field[rightmost].norm()-Field[leftmost].norm()<<"\n";


    // geometry->normalize(Vector3({0.0,0.0,0.0}),false);
    geometry->refreshQuantities();
    
    // Bead_1.Move_bead(backtrackstep,center);

    if(Save ){
    Sim_data << V_bar<<" "<< A_bar<<" "<< time <<" "<< V<<" " << A<<" " << E_Vol << " " << E_Sur << " " << E_Ben <<" " << E_Bead << " "<< grad_norm<<" " << backtrackstep << " \n";
    }
    system_time+=1;
    
    
  
  return backtrackstep;
}







// THis is the main integrate 
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
    if(Save || backtrackstep<0){
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


// THis integrate will do surface tension + volume preservation
double Mem3DG::integrate(double h, double V_bar,double P0,double KA,std::ofstream& Sim_data, double time,bool Save) {


    
    VertexData<Vector3> Force(*mesh);
    Force=buildFlowOperator(h,V_bar,P0,KA);
    
    
    // This force is the gradient basically
    double alpha=1e-2;

    // I have the forces for almost everything i just need the bead.

    size_t vindex;
    size_t Nvert=mesh->nVertices();


    // I want to print the Volume, the area, VOl_E, Area_E, Bending_E 
    double V = geometry->totalVolume();
    double A = geometry->totalArea();
    double D_P=-1*P0*(V-V_bar)/V_bar/V_bar;
    

    double E_Vol=E_Pressure(D_P,V,V_bar);
    double E_Sur=KA*A;
    double backtrackstep;

    backtrackstep=Backtracking(Force,D_P,V_bar,KA);
    // }
    if(Save || backtrackstep<0){
    Sim_data << V_bar<<" "<< time <<" "<< V<<" " << A<<" " << E_Vol << " " << E_Sur << " "<< grad_norm<<" " << backtrackstep<<" \n";
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
double D_P=-P0*(V-V_bar)/(V_bar*V_bar);
E_vol=E_Pressure(D_P,V,V_bar);



VertexData<Vector3> Calc_grad=D_P*OsmoticPressure();

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
// std::cout<<"\n";

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

VertexData<Vector3> Calc_grad=D_P*OsmoticPressure();

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


VertexData<Vector3> Calc_grad=lambda*SurfaceTension();

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


VertexData<Vector3> Calc_grad=KB*Bending(H_bar);

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



VertexData<Vector3> Calc_grad=KB*Bending(H_bar);

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
std::cout<<"The number of vertices is "<< N_vert <<" \n";
// for(size_t index=0; index<N_vert; index++){
size_t index;
for(Vertex v : mesh->vertices()){
  // std::cout<<"One vertex\n";
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



Eigen::Matrix2d Mem3DG::Face_sizing(Face f){
  Eigen::Matrix2d Sizing; //({0.0,0.0,0.0,0.0});
  // Eigen::Matrix2d Sizing{{0.0,0.0,0.0,0.0}};
  Sizing << 0.0, 0.0, 0.0,0.0;

  Eigen::Matrix3d Rotation;


  Vector3 Normal1 = geometry->faceNormal(f);
  Vector3 Axis1({0.0,1.0,0.0});
  
  // std::cout<<"THe norm of the normal is " << Normal1.norm2()<<" \n";
  Eigen::Vector3d Normal{Normal1.x,Normal1.y,Normal1.z};
  Eigen::Vector3d Axis{Axis1.x ,Axis1.y ,Axis1.z};

  Eigen::Vector3d uvw = Normal.cross(Axis);
  
  

  double rcos = Normal.dot(Axis);
  double rsin = uvw.norm();

  if(rsin>1e-7){
    uvw/=rsin;
  }

  // Eigen::Matrix3d V_x{{0, uvw(2), -1*uvw(1) },{-1*uvw(2), 0, uvw(0)},{uvw(1), -1*uvw(0), 0}};

  // Eigen::Matrix3d V_x{{0, -1*uvw(2), uvw(1) },{uvw(2), 0, -1*uvw(0)},{-1*uvw(1), uvw(0), 0}};
  // Eigen::Matrix3d V_x({0, -1*uvw(2), uvw(1) ,uvw(2), 0, -1*uvw(0),-1*uvw(1), uvw(0), 0});
  Eigen::Matrix3d V_x; 

  V_x << 0, -1*uvw(2), uvw(1) ,uvw(2), 0, -1*uvw(0),-1*uvw(1), uvw(0), 0;

  
  Rotation = rcos*Eigen::Matrix3d::Identity()+ rsin*V_x+ uvw*uvw.transpose()*(1-rcos);


  // We want to compare the outer proudct of two vectors
  // Eigen::Matrix3d outer{{uvw(0)*uvw(0),uvw(1)*uvw(0),uvw(2)*uvw(0) }, {uvw(0)*uvw(1),uvw(1)*uvw(1),uvw(2)*uvw(1)}, {uvw(0)*uvw(2),uvw(1)*uvw(2),uvw(2)*uvw(2)}};

  // std::cout<<"Outer is \n";
  // std::cout<<outer <<"\n";
  // std::cout<<"And the eigen one is "<<" \n";
  // std::cout<< uvw*uvw.transpose() <<"\n";

  Eigen::Vector3d trial({1.0,5.0,4.0});

  // std::cout<<"v_x is " << V_x << " \n";

  // std::cout<<"Trial has norm  " <<trial.norm() <<" \n";

  trial = Rotation*trial; 
// 
  // std::cout<<"Trial has now norm " << trial.norm() << "\n ";
  // std::cout<<"Trial now is " << trial <<" \n";




  Vector3 Pos;
  Eigen::Vector3d spare_vec;
  vector<Eigen::Vector2d> Positions(0);
  vector<int> index(0);

  for( Vertex v : f.adjacentVertices()){

    Pos = geometry->inputVertexPositions[v];
    spare_vec = Eigen::Vector3d{Pos.x, Pos.y, Pos.z}; 
    spare_vec = Rotation * spare_vec;
    index.push_back(v.getIndex());
    Positions.push_back({spare_vec(0),spare_vec(2)});
    // std::cout<<"THe y pos is " << spare_vec(1) <<"\t";

  }
  // std::cout<<"\n";

  // std::cout<<"The order is " << index[0] << " then " << index[1] << "then " << index[2]<<" \n";

  int counter = 0;
  for( Edge e : f.adjacentEdges()){

    // std::cout<<"The original edgelength is " << geometry->edgeLength(e) << " \n";
    

    Eigen::Vector2d e_mat = Positions[(counter+1)%3]-Positions[counter%3];
    
    // std::cout<<"THe projected length is " << e_mat.norm() << " \n";


    Eigen::Vector2d t_mat{ -1*e_mat(1), e_mat(0)};
    t_mat.normalize();


    double theta = geometry->dihedralAngle(e.halfedge());
    if(e.halfedge().face() != f){
      theta = geometry->dihedralAngle(e.halfedge().twin());
    }

    Sizing -= 1/2. * theta * e_mat.norm() * t_mat*t_mat.transpose();
    
    
    counter+=1;


  }

  Sizing/= geometry->faceArea(f);

  return Sizing;
}


FaceData<double> Mem3DG::Face_sizings(){

  FaceData<double> sizing(*mesh,0.0);
  double aspect_min = 0.2;
  double refine_angle = 0.7;

  // Lets create the face sizing
  Eigen::Matrix2d Face_sizing_mat ; //({0.0,0.0, 0.0,0.0});
  Face_sizing_mat << 0.0, 0.0, 0.0, 0.0;
  Eigen::EigenSolver<Eigen::Matrix2d> Solve;
  Eigen::Vector2d eigenvals;

  for(Face f : mesh->faces()){
    // For every face i need to do the sizing thingy
    Face_sizing_mat = Face_sizing(f) ;

    Face_sizing_mat = Face_sizing_mat.transpose()* Face_sizing_mat/(refine_angle * refine_angle);

    Solve.compute(Face_sizing_mat);

    eigenvals = Solve.eigenvalues().real();

    // 
    for(int i = 0; i <2  ;i++) eigenvals[i] = clamp(eigenvals[i],1.f/(0.2*0.2), 1.f/(0.001*0.001));

    double lmax = std::max(eigenvals[0]  ,eigenvals[1]);
    double lmin = lmax * aspect_min * aspect_min;

    for(int i = 0; i <2; i++) if(eigenvals[i] < lmin ) eigenvals[i] = lmin;

    
    sizing[f.getIndex()] = std::max( fabs(eigenvals[0]),fabs(eigenvals[1]));

  }




  return sizing;


}


VertexData<double> Mem3DG::Vert_sizing(FaceData<double> Face_sizings){

  VertexData<double> sizing(*mesh,0.0);

  // We need to create the vertsizing first.
  // WE have the sizing of every face
  for(Vertex v : mesh->vertices()){
    double sum = 0.0;
    for(Face f : v.adjacentFaces()){
      sum += Face_sizings[f]*geometry->faceArea(f)/3.;
    }
    sizing[v] = sum/geometry->vertexDualArea(v);
  }


  return sizing;
}

EdgeData<double> Mem3DG::Edge_sizing(VertexData<double> Vert_sizings){

  EdgeData<double> sizing(*mesh,0.0);

  Vertex v1;
  Vertex v2;
  for(Edge e: mesh->edges()){
    v1 = e.halfedge().vertex();
    v2 = e.halfedge().next().vertex();
    sizing[e] = geometry->edgeLength(e)*std::sqrt((Vert_sizings[v1]+Vert_sizings[v2])/2.0);
  }
  return sizing;
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