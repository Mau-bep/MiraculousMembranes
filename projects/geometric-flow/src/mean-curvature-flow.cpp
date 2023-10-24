// Implement member functions for MeanCurvatureFlow class.
#include "mean-curvature-flow.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
MeanCurvatureFlow::MeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

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
SparseMatrix<double> MeanCurvatureFlow::buildFlowOperator(const SparseMatrix<double>& M, double h) const {

    SparseMatrix<double> laplace=geometry->laplaceMatrix();
    
    size_t nvert=mesh->nVertices();
    SparseMatrix<double> resulting_op(nvert,nvert);
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for(size_t index; index<mesh->nVertices();index++)
    {
        tripletList.push_back(T(index,index,1.0) );
    }
    resulting_op.setFromTriplets(tripletList.begin(),tripletList.end());


    resulting_op=(M+h*laplace);
    resulting_op=resulting_op;
    // TODO
    return resulting_op; // placeholder
}





/*
 * Performs mean curvature flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void MeanCurvatureFlow::integrate(double h) {

    // TODO
    // Note: Geometry Central has linear solvers: https://geometry-central.net/numerical/linear_solvers/
    
    SparseMatrix<double> M=geometry->massMatrix();
    SparseMatrix<double> Op_flow=buildFlowOperator(M,h);
    Eigen::SparseLU<SparseMatrix<double> > solver;
    

    Vector<double> x_pos(mesh->nVertices());
    Vector<double> y_pos(mesh->nVertices());
    Vector<double> z_pos(mesh->nVertices());
    z_pos=x_pos+y_pos;
    // // rhs=other_m*( );
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
        
        if(index_actual==0){
            std::cout << "The position was "<< ","<<Position[0]<< ","<<Position[1]<< ","<<Position[2] << ","<<"and now the position will be "<< new_x[index_actual]<< ","<< new_y[index_actual] <<  ","<<new_z[index_actual] <<    "\n \n";

        }

        Update= { new_x[index_actual],new_y[index_actual],new_z[index_actual]  };
        geometry->inputVertexPositions[v] = Update ; // placeholder
    }
}