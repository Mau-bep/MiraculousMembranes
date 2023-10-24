// Implement member functions for IMem3DG class.
#include "mem-3dg_implicit.h"
#include <fstream>
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

    VertexData<double> Scalar_MC(*mesh,0.0);
    size_t index;
    for(Vertex v : mesh->vertices()){
        index=v.getIndex();
        Scalar_MC[index]=geometry->scalarMeanCurvature(v)/geometry->barycentricDualArea(v);
    }
    Vector<Vector3> Face_normals(mesh->nFaces());
    for(Face f : mesh->faces()){
        Face_normals[f.getIndex()]=geometry->faceNormal(f);
    }


    SparseMatrix<double> laplace=geometry->laplaceMatrix3d();
    SparseMatrix<double> Normal=geometry->NormalFlowMat();
    SparseMatrix<double> MeanH= geometry->MeanHFlowMat3d(Scalar_MC,Face_normals);
    SparseMatrix<double> Schlafi = geometry->Schlafi3d(Scalar_MC,Face_normals);
    SparseMatrix<double> Gauss = geometry->GaussFlowMat3d(Scalar_MC);
    size_t nvert=mesh->nVertices();
    SparseMatrix<double> resulting_op(3*nvert,3*nvert);
    


    // typedef Eigen::Triplet<double> T;
    // std::vector<T> tripletList;
    // for(size_t index=0; index<3*mesh->nVertices();index++)
    // {
    //     tripletList.push_back(T(index,index,1.0) );
        
    // }
    // resulting_op.setFromTriplets(tripletList.begin(),tripletList.end());
    double V= geometry->totalVolume();
    double DP=-P0*(V-V_bar)/V_bar/V_bar;
    double A=geometry->totalArea();
    double A_bar=4*PI*pow(3*V_bar/(4*PI*nu),2.0/3.0);
    double lambda=KA*(A-A_bar)/A_bar;
    // lambda=KA;
    // std::cout<<"THe value of DP is "<<DP<<"\n";
    // resulting_op=(M);
    // resulting_op=(M+h*lambda*laplace+h*DP*M*Normal);// + h*KB*M*MeanH);
    resulting_op=(M+h*lambda*laplace+h*DP*M*Normal + M*h*KB*(Schlafi+Gauss+MeanH ));
    // resulting_op=(M+M*h*KB*(Schlafi+Gauss+Meanh));
    // resulting_op=resulting_op;
    // TODO
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
    

    SparseMatrix<double> M=geometry->massMatrix3d();


    SparseMatrix<double> Op_flow=buildFlowOperator(M,h,nu,V_bar,P0,KA,KB);
    Eigen::SparseLU<SparseMatrix<double> > solver;

    size_t nvert=mesh->nVertices();


    Vector<double> X_vect(3*nvert);
    
    for(Vertex v: mesh-> vertices())
        {
        Vector3 Position =geometry->inputVertexPositions[v.getIndex()];
        X_vect.coeffRef(v.getIndex())=Position.x;
        X_vect.coeffRef(v.getIndex()+nvert)=Position.y;
        X_vect.coeffRef(v.getIndex()+2*nvert)=Position.z;
        }


    // Vector<double> result = solver.solve(rhs);
    solver.compute(Op_flow);
    

    Vector<double> rhs=M*(X_vect);
   
    Vector<double> new_X= solver.solve(rhs);


    // Note: Update positions via geometry->inputVertexPositions
    Vector3 Update;
    size_t index_actual;
    for (Vertex v : mesh->vertices()) {

        Vector3 Position =geometry->inputVertexPositions[v.getIndex()];

        index_actual=v.getIndex();
        Update= { new_X[index_actual],new_X[index_actual+nvert],new_X[index_actual+2*nvert]  };
        // if(index_actual==0){
        //     std::cout << "The position was "<< ","<<Position[0]<< ","<<Position[1]<< ","<<Position[2] << ","<<"and now the position will be "<< Update.x<< ","<< Update.y <<  ","<<Update.z <<    "\n \n";

        // }
        
        geometry->inputVertexPositions[v] = Update ; // placeholder
    }
}