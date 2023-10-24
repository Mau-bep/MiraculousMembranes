// Implement member functions for NOrmalFlow class.
#include "willmore-flow.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
WillmoreFlow::WillmoreFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {
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
SparseMatrix<double> WillmoreFlow::buildFlowOperator(const SparseMatrix<double>& Q,const SparseMatrix<double>& M, double h) const {

    // SparseMatrix<double> FlowMat=geometry->FullQWillmore();
    size_t nvert=mesh->nVertices();
    SparseMatrix<double> resulting_op(nvert,nvert);

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

    for(size_t index; index<nvert;index++)
    {
        tripletList.push_back(T(index,index,1.0) );
    }
    resulting_op.setFromTriplets(tripletList.begin(),tripletList.end());
    
    



    resulting_op=(M-h*Q);
    // TODO
    return resulting_op; // placeholder
}




/*
 * Performs willmore flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void WillmoreFlow::integrate(double h) {

    
    SparseMatrix<double> M=geometry->massMatrix();
    SparseMatrix<double> Q = geometry->FullQWillmore();
    SparseMatrix<double> Op_flow=buildFlowOperator(Q,M,h);
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
    
    WillE= x_pos.transpose()*(Q*x_pos) ;
    WillE= WillE + y_pos.transpose()*(Q*y_pos);
    WillE= WillE + z_pos.transpose()*(Q*z_pos);
    
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
        
        // if(index_actual==0){
        //     std::cout << "The position was "<< ","<<Position[0]<< ","<<Position[1]<< ","<<Position[2] << ","<<"and now the position will be "<< new_x[index_actual]<< ","<< new_y[index_actual] <<  ","<<new_z[index_actual] <<    "\n \n";

        // }

        Update= { new_x[index_actual],new_y[index_actual],new_z[index_actual]  };
        geometry->inputVertexPositions[v] = Update ; // placeholder
    }
    for(Vertex v: mesh-> vertices())
        {
        Vector3 Position =geometry->inputVertexPositions[v.getIndex()];
    
        x_pos.coeffRef(v.getIndex())=Position.x;
        y_pos.coeffRef(v.getIndex())=Position.y;
        z_pos.coeffRef(v.getIndex())=Position.z;
        }


    WillE= new_x.transpose()*(Q*new_x) ;
    WillE= WillE + new_y.transpose()*(Q*new_y);
    WillE= WillE + new_z.transpose()*(Q*new_z);
    std::cout<< "After the iteration we get the following results\n";
    std::cout << "The Willmore energy is " << WillE << "THe area is "<< geometry->totalArea()<< "and the Volume is " << geometry->totalVolume()<< "\n\n\n";
    



}