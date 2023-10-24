// Implement member functions for ModifiedMeanCurvatureFlow class.
#include "modified-mean-curvature-flow.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
ModifiedMeanCurvatureFlow::ModifiedMeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry, and A (Laplace matrix)
    mesh = inputMesh;
    geometry = inputGeo;
      

    this->A = geometry->laplaceMatrix(); // placeholder
}

/*
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
SparseMatrix<double> ModifiedMeanCurvatureFlow::buildFlowOperator(const SparseMatrix<double>& M, double h) const {
    // TODO

    size_t nvert=mesh->nVertices();
    SparseMatrix<double> resulting_op(nvert,nvert);
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    double eps = 1e-5;
    for(size_t index; index<mesh->nVertices();index++)
    {
        tripletList.push_back(T(index,index,1.0) );
    }
    resulting_op.setFromTriplets(tripletList.begin(),tripletList.end());


    resulting_op=(M+h*this->A);
    // resulting_op=M*resulting_op;


    return resulting_op; // placeholder
}