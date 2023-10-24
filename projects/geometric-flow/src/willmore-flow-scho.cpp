// Implement member functions for NOrmalFlow class.
#include "willmore-flow-scho.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
WillmoreFlowScho::WillmoreFlowScho(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {
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
SparseMatrix<double> WillmoreFlowScho::buildFlowOperator(const SparseMatrix<double>& M, double h) const {

    SparseMatrix<double> FlowMat=geometry->FullKWillmore();
    size_t nvert=mesh->nVertices();
    SparseMatrix<double> resulting_op(nvert,nvert);

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

    for(size_t index; index<nvert;index++)
    {
        tripletList.push_back(T(index,index,1.0) );
    }
    resulting_op.setFromTriplets(tripletList.begin(),tripletList.end());
    

    return (resulting_op+h*FlowMat); // placeholder
}




/*
 * Performs willmore flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void WillmoreFlowScho::integrate(double h) {

    
    SparseMatrix<double> M=geometry->massMatrix();
    SparseMatrix<double> Op_flow=buildFlowOperator(M,h);
    
    size_t nvert=mesh->nVertices();
    
    
    Eigen::SparseLU<SparseMatrix<double> > solver;
    

    Vector<double> x_pos(mesh->nVertices());
    Vector<double> y_pos(mesh->nVertices());
    Vector<double> z_pos(mesh->nVertices());
    
    for(Vertex v: mesh-> vertices())
        {
        Vector3 Position =geometry->inputVertexPositions[v.getIndex()];
    
        x_pos.coeffRef(v.getIndex())=Position.x;
        y_pos.coeffRef(v.getIndex())=Position.y;
        z_pos.coeffRef(v.getIndex())=Position.z;
        }

    solver.compute(Op_flow);
    
  
    Vector<double> rhs_1=(x_pos);
        
    Vector<double> new_x= solver.solve(rhs_1);

    Vector<double> rhs_2=(y_pos);

    Vector<double> new_y= solver.solve(rhs_2);

    Vector<double> rhs_3=(z_pos );
    Vector<double> new_z= solver.solve(rhs_3);
    



    // rhs_1= (FlowMat-identityMatrix<double>(nvert)  +h*M*FlowMat)*new_x + x_pos;
    // rhs_2= (FlowMat-identityMatrix<double>(nvert)  +h*M*FlowMat)*new_y + y_pos;
    // rhs_3= (FlowMat-identityMatrix<double>(nvert)  +h*M*FlowMat)*new_z + z_pos;
    
    // new_x=solve2.solve(rhs_1);
    // new_y=solve2.solve(rhs_2);
    // new_z=solve2.solve(rhs_3);
    
//    DenseMatrix<double> InverseFlowMat= FlowMat.toDense().inverse();


    // size_t n_iter=0;
    // for (size_t n_iter; n_iter<3; n_iter++)
    // {
    //  new_x=( identityMatrix<double>(nvert)-InverseFlowMat+ h*InverseFlowMat*M*FlowMat )*new_x + InverseFlowMat*x_pos;
    //  new_y=( identityMatrix<double>(nvert)-InverseFlowMat+ h*InverseFlowMat*M*FlowMat )*new_y + InverseFlowMat*y_pos;
    //  new_z=( identityMatrix<double>(nvert)-InverseFlowMat+ h*InverseFlowMat*M*FlowMat )*new_z + InverseFlowMat*z_pos;
    

    // }




    // Note: Update positions via geometry->inputVertexPositions
    Vector3 Update;
    size_t index_actual;
    for (Vertex v : mesh->vertices()) {

        Vector3 Position =geometry->inputVertexPositions[v.getIndex()];

        index_actual=v.getIndex();
        
        if(index_actual==0){
            std::cout << "The position was "<< ","<<Position[0]<< ","<<Position[1]<< ","<<Position[2] << ","<<"and now the position will be "<< new_x[index_actual]<< ","<< new_y[index_actual] <<  ","<<new_z[index_actual] <<    "\n \n";

        }

        Update= { new_x[index_actual],new_y[index_actual],new_z[index_actual]  };
        geometry->inputVertexPositions[v] = Update ; // placeholder
    }

}