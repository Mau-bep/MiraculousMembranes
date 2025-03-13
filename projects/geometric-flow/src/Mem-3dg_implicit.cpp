// Implement member functions for IMem3DG class.
#include "Mem-3dg_implicit.h"
#include "Beads.h"
#include <fstream>
#include <chrono>
#include <Eigen/Core>

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
using namespace geometrycentral;
using namespace geometrycentral::surface;
IMem3DG::IMem3DG(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
    

}

/* 
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
SparseMatrix<double> IMem3DG::buildFlowOperator(const SparseMatrix<double>& M, double h,double nu,double V_bar,double P0, double KA,double KB) const {

    // VertexData<double> Scalar_MC(*mesh,0.0);
    size_t index;
    size_t nvert = mesh->nVertices();
    
    
    Vector<Vector3> Face_normals(mesh->nFaces());

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    
    SparseMatrix<double> Inv_M(nvert,nvert);
    for(size_t index; index<nvert;index++)
    {
        tripletList.push_back(T(index,index, 1.0/(M.coeff(index,index))  ) );
       
    }
    Inv_M.setFromTriplets(tripletList.begin(),tripletList.end());
    SparseMatrix<double> laplacian=geometry->laplaceMatrix();
    

    SparseMatrix<double> resulting_op(nvert,nvert);
    
    double A=geometry->totalArea();
    double A_bar=4*PI*pow(3*V_bar/(4*PI*nu),2.0/3.0);
    // double lambda=KA;
    double lambda=KA*(A-A_bar)/(A_bar*A_bar);
    // std::cout<<"Lambda is "<< lambda <<"\n";
    // resulting_op=(M+h*lambda*laplacian);
    // resulting_op=(M+h*lambda*laplacian+h*KB*laplacian*Inv_M*laplacian);
    
    resulting_op = (M+ h * lambda * laplacian + h * KB * (laplacian).transpose()*Inv_M*laplacian);
 
    return resulting_op;
    // return resulting_op; // placeholder
}





/*
 * Performs mean curvature flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void IMem3DG::integrate(double h,double nu,double V_bar,double P0, double KA,double KB) {

    // TODO
    // Note: Geometry Central has linear solvers: https://geometry-central.net/numerical/linear_solvers/
    

    // Lets recenter the mesh

    // geometry->normalize(Vector3({0.0,0.0,0.0}),false);



    size_t nvert=mesh->nVertices();
    SparseMatrix<double> M=geometry->massMatrix();
    
    SparseMatrix<double> Op_flow(nvert,nvert);
    Op_flow=buildFlowOperator(M,h,nu,V_bar,P0,KA,KB);
    Eigen::SparseLU<SparseMatrix<double> > solver;


    Vector<double> x_pos(nvert);
    Vector<double> y_pos(nvert);
    Vector<double> z_pos(nvert);
    
    for(Vertex v: mesh-> vertices())
        {
        Vector3 Position =geometry->inputVertexPositions[v.getIndex()];
        x_pos.coeffRef(v.getIndex())=Position.x;
        y_pos.coeffRef(v.getIndex())=Position.y;
        z_pos.coeffRef(v.getIndex())=Position.z;
        }


    // Vector<double> result = solver.solve(rhs);
    solver.compute(Op_flow);

    Vector<double> rhs_1=M*(x_pos);
   
    Vector<double> new_x= solver.solve(rhs_1);

    Vector<double> rhs_2=M*(y_pos);

    Vector<double> new_y= solver.solve(rhs_2);

    Vector<double> rhs_3=M*(z_pos );
    
    Vector<double> new_z= solver.solve(rhs_3);

    // Note: Update positions via geometry->inputVertexPositions
    Vector3 Update;
    size_t index_actual;
    for (Vertex v : mesh->vertices()) {

        Vector3 Position =geometry->inputVertexPositions[v.getIndex()];

        index_actual=v.getIndex();
       
        Update= { new_x[index_actual],new_y[index_actual],new_z[index_actual]  };
        geometry->inputVertexPositions[v] = Update ;
    }

    geometry->refreshQuantities();
    double Vol=geometry->totalVolume();
    
    double k =pow(V_bar/Vol,1/3.0);
    // std::cout<<"The value of k is" << k <<"\n";
    for(Vertex v: mesh->vertices()){
    geometry->inputVertexPositions[v]=k*geometry->inputVertexPositions[v];
        
    }
    geometry->refreshQuantities();
    // double NewVol=geometry->totalVolume();

    // std::cout<<"THe old volume was "<< Vol << "and the current volume is"<< NewVol<<" \n";
}



// We will add the Energy Gradient here 


Vector3 IMem3DG::computeHalfedgeMeanCurvatureVector(Halfedge he) const {
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



Vector3 IMem3DG::computeHalfedgeGaussianCurvatureVector(Halfedge he) const {
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


Vector3 IMem3DG::dihedralAngleGradient(Halfedge he, Vertex v) const {
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


VertexData<Vector3> IMem3DG::Bending(double H0) const {

    
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


VertexData<Vector3> IMem3DG::SurfaceGrad() const {

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

VertexData<Vector3> IMem3DG::OsmoticPressure() const {

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




VertexData<Vector3> IMem3DG::Sobolev_operator() const {


    // Ok so i have the gradient i want to calculate

    // std::cout<<"We are doing the sobolev operator ! \n";

    SparseMatrix<double> L = geometry->laplaceMatrix();
    SparseMatrix<double> M = geometry->massMatrix();
    SparseMatrix<double> Inv_M(mesh->nVertices(),mesh->nVertices());
    for(size_t index; index<mesh->nVertices();index++)
    {
        Inv_M.coeffRef(index,index)= 1.0/(M.coeff(index,index));
       
    }

    SparseMatrix<double> J = L.transpose()*Inv_M*L;
    // I need the constraint gradientsx
    VertexData<Vector3> grad_sur = SurfaceGrad();
    VertexData<Vector3> grad_vol = OsmoticPressure();
    VertexData<Vector3> grad_ben = Bending(0.0);

    // i CAN MAYBE MULTIPLY BY THE SURFACE AREAS 
    double area;
    for(size_t i = 0; i < mesh->nVertices(); i++){
        // area = geometry->barycentricDualArea(mesh->vertex(i));
        grad_ben[mesh->vertex(i)] = grad_ben[mesh->vertex(i)]*M.coeff(i,i);
    }


    // std::cout<<"Grad ben \n";

    // std::cout<<"We are debugginggg \n";
    size_t N_vert=mesh->nVertices();


    int Num_constraints = 4;
    SparseMatrix<double> S(N_vert*3+Num_constraints,N_vert*3+Num_constraints);

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    int highest_row = 0;
    int highest_col = 0;
    // We add the constraints first

    int Constraint_number = 0;


    for(size_t index; index<mesh->nVertices();index++)
    {
        Constraint_number = 0;
    // //  Area constraint 
    //     tripletList.push_back(T(3*index,3*N_vert+Constraint_number, grad_sur[index].x ) );
    //     tripletList.push_back(T(3*index+1, 3*N_vert+Constraint_number, grad_sur[index].y ) );
    //     tripletList.push_back(T(3*index+2, 3*N_vert+Constraint_number, grad_sur[index].z ) );
    //     tripletList.push_back(T(3*N_vert+Constraint_number,3*index, grad_sur[index].x ) );
    //     tripletList.push_back(T(3*N_vert+Constraint_number, 3*index+1, grad_sur[index].y ) );
    //     tripletList.push_back(T(3*N_vert+Constraint_number,3*index+2, grad_sur[index].z ) );
            // Constraint_number++;
    // Volume Constraint
        tripletList.push_back(T(3*index, 3*N_vert+Constraint_number, grad_vol[index].x ) );
        tripletList.push_back(T(3*index+1, 3*N_vert+Constraint_number, grad_vol[index].y ) );
        tripletList.push_back(T(3*index+2, 3*N_vert+Constraint_number, grad_vol[index].z ) );
        tripletList.push_back(T(3*N_vert+Constraint_number, 3*index, grad_vol[index].x ) );
        tripletList.push_back(T(3*N_vert+Constraint_number, 3*index+1, grad_vol[index].y ) );
        tripletList.push_back(T(3*N_vert+Constraint_number, 3*index+2, grad_vol[index].z ) );

    // Position Constraint 
        Constraint_number++;
        tripletList.push_back(T(3*index, 3*N_vert +Constraint_number, 1 ) );
        tripletList.push_back(T(3*N_vert +Constraint_number, 3*index, 1 ) );

        Constraint_number++;
        tripletList.push_back(T(3*index + 1, 3*N_vert +Constraint_number, 1 ) );
        tripletList.push_back(T(3*N_vert +Constraint_number, 3*index +1, 1 ) );

        Constraint_number++;
        tripletList.push_back(T(3*index + 2 , 3*N_vert+Constraint_number, 1 ) );
        tripletList.push_back(T(3*N_vert+Constraint_number, 3*index + 2, 1 ) );



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

    for( size_t k = 0; k < L.outerSize(); ++k ) {
        for( SparseMatrix<double>::InnerIterator it(L,k); it; ++it ) {
            value = it.value();
            row = it.row();
            col = it.col();
            tripletList.push_back(T(3*row,3*col,value));
            tripletList.push_back(T(3*row+1,3*col+1,value));
            tripletList.push_back(T(3*row+2,3*col+2,value));
            if( 3*row > highest_row) highest_row = 3*row;
            if( 3*col > highest_col) highest_col = 3*col;
        
        }
    }

    // std::cout<<"Before setting from tripplets\n";
    // std::cout<<"THe number of vertices is " << N_vert <<"\n";
    // std::cout<<"THe highest expected col is " << 3*N_vert+2 <<" \n";
    // std::cout<<"THe calculated col is " << highest_col <<" \n";
    // std::cout<<"THe highest expected row is " << 3*N_vert+2 <<" \n";
    // std::cout<<"THe calculated row is " << highest_row <<" \n";
    S.setFromTriplets(tripletList.begin(),tripletList.end());
    // Now u have a Sovolev operator

    // std::cout<<"AFTER setting from tripplets\n";
    // We need to create the vector of the RHS

    Vector<double> RHS(N_vert*3+Num_constraints);
    int highest_idx = 0;
    for(size_t index; index<mesh->nVertices();index++)
    {
        RHS.coeffRef(3*index)=grad_ben[index].x;
        RHS.coeffRef(3*index+1)=grad_ben[index].y;
        RHS.coeffRef(3*index+2)=grad_ben[index].z;
        // if( 3*index +2 > highest_idx) highest_idx = 3*index+2;
    }

    // std::cout<<"The highest index is " << highest_idx <<"\n";
    // std::cout<<"CREATING RHS\n";
    for(size_t index = 0; index<Num_constraints;index++)
    {
        RHS.coeffRef(3*N_vert+index)=0;
    }
    // RHS.coeffRef(3*N_vert)=0;
    // RHS.coeffRef(3*N_vert+1)=0;
    // RHS.coeffRef(3*N_vert+2)=0;
    // RHS.coeffRef(3*N_vert+3)=0;
    // RHS.coeffRef(3*N_vert+4)=0;

    // std::cout<<"RHS IMPLEMENTED\n";

    // std::cout<<"The RHS is " << RHS << "\n";
    // We have The RHS and the matrix, its tiiiime 
    // std::cout<<"Lets solve \n";
    Eigen::SparseLU<SparseMatrix<double> > solver;
    // std::cout<<"Solver defined \n";
    
    solver.compute(S);
    // std::cout<<"Solver computed \n";
    Vector<double> result = solver.solve(RHS);

    // std::cout<<"Result retrieved\n";
    VertexData<Vector3> Final_Force(*mesh);
    for(size_t index; index<mesh->nVertices();index++)
    {
        Final_Force[index]=Vector3{result.coeff(3*index),result.coeff(3*index+1),result.coeff(3*index+2)};
    } 

    // std::cout<<"The two lambda for the constraints are "<< result[N_vert] << " and "<< result[N_vert+1] <<"\n";
    return Final_Force;
    



}

double IMem3DG::E_Bending(double H0,double KB) const{
    size_t index;
    double Eb=0;
    double H;
    double r_eff2;
    Vector3 Pos;
    for(Vertex v : mesh->vertices()) {
        //boundary_fix
        if(v.isBoundary()) continue;
        
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

double IMem3DG::Backtracking(VertexData<Vector3> Force, double h ) {

    double ts = h;


    // std::cout<<"Bending energy\n";
    double E_bend = E_Bending(0.0,1.0);


    Vector3 Total_force;
    double previous_energy = E_bend;
    // std::cout<<"Previous energy is " << E_bend <<" \n";
    double NewE = 0.0;

    VertexData<Vector3> initial_positions(*mesh);
    // std::cout<<"Saving pos\n";
    initial_positions = geometry->inputVertexPositions;

    // std::cout<<"Evolving pos\n";
    geometry->inputVertexPositions+=h*Force;
    // std::cout<<"recentering\n";
    geometry->refreshQuantities();

    Total_force = Vector3({0.0,0.0,0.0});
    double Force_norm2 = 0.0;
    // std::cout<<"FORCEQUANTITIES\n";
    for(Vertex v : mesh->vertices()) {
        Total_force += Force[v];
        Force_norm2 += Force[v].norm2();
    }
    // std::cout<<"calculated\n";
    NewE = 0.0;
    // std::cout<<"Bending again\n";
    NewE = E_Bending(0.0,1.0);
    // std::cout<<"The energy is "<< NewE <<"\n";
    // if( std::fabs(NewE-previous_energy) <=1e-5) std::cout<<"The two energies are the same "<< NewE << " and "<< previous_energy <<"\n";  
    // std::cout<<"\n";
    while(true){
        // std::cout<<"eval\n";
        if(NewE <= previous_energy - 1e-4*ts*Force_norm2){
            // std::cout<<"Satisfied, next step\n";
            break;
        }

        if( fabs(NewE - previous_energy)/previous_energy < 1e-9 ){
            std::cout <<"Relative energy difference is small "<< fabs(NewE - previous_energy)/previous_energy <<" stopping sim\n"; 
            ts = -1;
            break;
        }
        // std::cout<<"Reinitialization\n";
        geometry->inputVertexPositions = initial_positions;
        ts = 0.5*ts;
        // std::cout<<"adding force\n";
        geometry->inputVertexPositions += ts*Force;
        // std::cout<<"refreshing\n";
        geometry->refreshQuantities();
        // std::cout<<"Bendinggg\n";
        NewE = E_Bending(0.0,1.0);
        // std::cout<<"check\n";
        if(ts <1e-7){
            std::cout<<"timestep very small, stopping this \n";
            ts = -1;
            break;
        }

    }
    // Ok we want to backtrack
    std::cout<<"Timestep is " << ts <<" \n";
    return ts;


}


double IMem3DG::integrate_Sob(double h,double Kb, std::ofstream& Sim_data,bool Save_output_data){

    // geometry->normalize(Vector3({0.0,0.0,0.0}),false);
    // geometry->refreshQuantities();

    VertexData<Vector3> Sobolev_gradient = this->Sobolev_operator();

    // I have the gradient
    // std::cout<<"Backtracking\n";
    double ts = Backtracking(Sobolev_gradient,h);
    // std::cout<<"Backtracking done\n";

    double V = geometry->totalVolume();
    double A = geometry->totalArea();
    double Bending = E_Bending(0.0,Kb);

    Sim_data << V << " " << A << " " << Bending << " " << ts << "\n";
    std::cout<<"THe bending energy is " << Bending << " \n";
    // geometry->inputVertexPositions += h*Sobolev_gradient;
    // geometry->normalize(Vector3({0.0,0.0,0.0}),false);
    // geometry->refreshQuantities();
    return ts;


}