// Implement member functions for IMem3DG class.
#include "Mem-3dg_implicit.h"
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
    double lambda=KA*(A-A_bar)/A_bar;
    
    
    resulting_op=(M+h*lambda*laplacian+h*KB*laplacian*Inv_M*laplacian);
    
 
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