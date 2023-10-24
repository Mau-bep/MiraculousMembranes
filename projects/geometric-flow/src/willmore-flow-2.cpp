// Implement member functions for NOrmalFlow class.
#include "willmore-flow-2.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
WillmoreFlow2::WillmoreFlow2(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {
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
SparseMatrix<double> WillmoreFlow2::buildFlowOperator(const SparseMatrix<double>& M, double h) const {

    // SparseMatrix<double> FlowMat=geometry->FullQWillmore();
    SparseMatrix<double> Laplacian=geometry->laplaceMatrix();

   
    size_t nvert=mesh->nVertices();
    SparseMatrix<double> Inv_M(nvert,nvert);
    SparseMatrix<double> resulting_op(nvert,nvert);

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

    for(size_t index; index<nvert;index++)
    {
        tripletList.push_back(T(index,index, 1.0/(M.coeff(index,index))  ) );
    }

    Inv_M.setFromTriplets(tripletList.begin(),tripletList.end());
    SparseMatrix<double> FlowOp = Laplacian.transpose()*Inv_M*Laplacian;
    
    // resulting_op=FlowOp;
    // TODO
    return FlowOp; // placeholder
}




/*
 * Performs willmore flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void WillmoreFlow2::integrate(double h) {

    
    SparseMatrix<double> M=geometry->massMatrix();
    // SparseMatrix<double> Q = geometry->FullQWillmore();
    SparseMatrix<double> Only_flow_Op=buildFlowOperator(M,h);
    SparseMatrix<double> Op_flow= M+h*Only_flow_Op;
    

    


    size_t nvert=mesh->nVertices();
    double WillE=0.0;
    
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
    WillE=x_pos.transpose()*(Only_flow_Op*x_pos);
    WillE= WillE + y_pos.transpose()*(Only_flow_Op*y_pos);
    WillE= WillE + z_pos.transpose()*(Only_flow_Op*z_pos);
    std::cout << "The Willmore energy is " << WillE << "THe area is "<< geometry->totalArea()<< "and the Volume is " << geometry->totalVolume()<< "\n";
    
  
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
        
        if(index_actual==0){
            std::cout << "The position was "<< ","<<Position[0]<< ","<<Position[1]<< ","<<Position[2] << ","<<"and now the position will be "<< new_x[index_actual]<< ","<< new_y[index_actual] <<  ","<<new_z[index_actual] <<    "\n \n";

        }

        Update= { new_x[index_actual],new_y[index_actual],new_z[index_actual]  };
        geometry->inputVertexPositions[v] = Update ; // placeholder
    }
    WillE=new_x.transpose()*(Only_flow_Op*new_x);
    WillE= WillE + new_y.transpose()*(Only_flow_Op*new_y);
    WillE= WillE + new_z.transpose()*(Only_flow_Op*new_z);
    std::cout<< "After the iteration we get the following results\n";
    
    std::cout << "The Willmore energy is " << WillE << "THe area is "<< geometry->totalArea()<< "and the Volume is " << geometry->totalVolume()<< "\n\n\n";
    


}