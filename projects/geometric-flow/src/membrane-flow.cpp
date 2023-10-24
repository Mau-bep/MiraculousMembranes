// Implement member functions for MeanCurvatureFlow class.
#include "membrane-flow.h"
#include "mean-curvature-flow.h"
#include "modified-mean-curvature-flow.h"
#include "normal-flow.h"
#include "gauss-curvature-flow.h"
#include "willmore-flow.h"

#include <Eigen/Core>
#include <Eigen/Cholesky>


/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
MembraneFlow::MembraneFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

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
SparseMatrix<double> MembraneFlow::buildFlowOperator(const SparseMatrix<double>& M,double sigma, double kappa, double H0 ,double h) const {

    SparseMatrix<double> Laplace = geometry->laplaceMatrix();
    SparseMatrix<double> Gauss = geometry->GaussFlowMat();

    size_t nvert=mesh->nVertices();

    SparseMatrix<double> Inv_M(nvert,nvert);
    SparseMatrix<double> resulting_op(nvert,nvert);

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    std::vector<T> tripletList2;

    for(size_t index; index<nvert;index++)
    {
        tripletList.push_back(T(index,index, 1.0/(M.coeff(index,index))  ) );
        tripletList2.push_back(T(index,index,1.0) );
    }

    Inv_M.setFromTriplets(tripletList.begin(),tripletList.end());
    SparseMatrix<double> Willmore = Laplace*Inv_M*Laplace;
    

    resulting_op.setFromTriplets(tripletList2.begin(),tripletList2.end());


    resulting_op=(M+h*( (sigma+2*kappa*H0*H0)* Laplace +  4*kappa*H0*Gauss  + 2*kappa*Willmore ));
    // resulting_op=resulting_op;
    // TODO
    return resulting_op; // placeholder
}




/*
 * Performs membrane flow.
 *
 * Input: The timestep <h>, and membrane parameters.
 * Returns:
 */
void MembraneFlow::integrate(double h, double sigma, double kappa, double H0, double P,double V0) {

    // TODO
    // Note: Geometry Central has linear solvers: https://geometry-central.net/numerical/linear_solvers/
    
    SparseMatrix<double> M=geometry->massMatrix();
    SparseMatrix<double> Op_flow=buildFlowOperator(M,sigma,kappa,H0,h);
    // Eigen::SparseLU<SparseMatrix<double> > solver;
    Eigen::SimplicialLDLT<SparseMatrix<double>> solver;   

    Vector<double> x_pos(mesh->nVertices());
    Vector<double> y_pos(mesh->nVertices());
    Vector<double> z_pos(mesh->nVertices());
    // // rhs=other_m*( );
    for(Vertex v: mesh-> vertices())
        {
        Vector3 Position =geometry->inputVertexPositions[v.getIndex()];
        x_pos.coeffRef(v.getIndex())=Position.x;
        y_pos.coeffRef(v.getIndex())=Position.y;
        z_pos.coeffRef(v.getIndex())=Position.z;
        }

    // Vector<double> result = solver.solve(rhs);
    solver.analyzePattern(Op_flow);
    solver.factorize(Op_flow);
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
        
        // if(index_actual==0){
        //     std::cout << "The position was "<< ","<<Position[0]<< ","<<Position[1]<< ","<<Position[2] << ","<<"and now the position will be "<< new_x[index_actual]<< ","<< new_y[index_actual] <<  ","<<new_z[index_actual] <<    "\n \n";

        // }

        Update= { new_x[index_actual],new_y[index_actual],new_z[index_actual]  };
        geometry->inputVertexPositions[v] = Update ; // placeholder
    }

    // Now that we have updated the positions we need to rescale

    if(P>1.0){

    double k = std::pow(V0/geometry->totalVolume(),1/3.0);
    for(Vertex v: mesh->vertices()) geometry->inputVertexPositions[v]= geometry->inputVertexPositions[v]*k;
    }
    



}