// Implement member functions for NOrmalFlow class.
#include "normal-flow.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
NormalFlow::NormalFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {
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
SparseMatrix<double> NormalFlow::buildFlowOperator(const SparseMatrix<double>& M, double h) const {

    SparseMatrix<double> FlowMat=geometry->NormalFlowMat();
    size_t nvert=mesh->nVertices();
    SparseMatrix<double> resulting_op(3*nvert,3*nvert);
    size_t j;
    size_t k;
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    Vector3 fj;
    Vector3 fk;
    for(size_t index; index<3*nvert;index++)
    {
        // Vertex vi= mesh->vertex[index];
        // for(Halfedge he : vi.outgoingHalfedges()){
        //     // Here i can have  j and k 
        //     j=he.tipVertex().getIndex();
        //     k=he.next().tipVertex().getIndex();
        //     fj=geometry->inputVertexPositions[j];
        //     fk=geometry->inputVertexPositions[k]
        //     // This is the first component
        //     tripletList.push_back(T(3*index,3*j+1,fk.z/12))
        //     tripletList.push_back(T(3*index,3*k+2,fj.y/12))
            
        //     tripletList.push_back(T(3*index,3*j+2,-1*fk.y/12))
        //     tripletList.push_back(T(3*index,3*k+1,-1*fj.z/12))
            
        //     //Now is the second component

        //     tripletList.push_back(T(3*index+1,3*j+2,fk.x/12))
        //     tripletList.push_back(T(3*index+1,3*k,fj.z/12))

        //     tripletList.push_back(T(3*index+1,3*j,-1*fk.z/12))
        //     tripletList.push_back(T(3*index+1,3*k+2,-1*fj.x/12))

        //     // Now the third component

        //     tripletList.push_back(T(3*index+1,3*j,fk.y/12))
        //     tripletList.push_back(T(3*index+1,3*k+1,fj.x/12))

        //     tripletList.push_back(T(3*index+1,3*j+1,-1*fk.x/12))
        //     tripletList.push_back(T(3*index+1,3*k,-1*fj.y/12))
            

        // }
        tripletList.push_back(T(index,index,1.0) );
    }
    resulting_op.setFromTriplets(tripletList.begin(),tripletList.end());
    
    resulting_op=(resulting_op+h*FlowMat);
    // TODO
    return resulting_op; // placeholder
}


SparseMatrix<double> NormalFlow::buildFlowOperatorx(const SparseMatrix<double>& M, double h) const {

    SparseMatrix<double> FlowMat=geometry->NormalFlowMatx();
    size_t nvert=mesh->nVertices();
    SparseMatrix<double> resulting_op(nvert,nvert);
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

    for(size_t index; index<mesh->nVertices();index++)
    {
        tripletList.push_back(T(index,index,1.0) );
    }
    resulting_op.setFromTriplets(tripletList.begin(),tripletList.end());
    
    resulting_op=(M+h*FlowMat);
    // TODO
    return resulting_op; // placeholder
}

SparseMatrix<double> NormalFlow::buildFlowOperatory(const SparseMatrix<double>& M, double h) const {

    SparseMatrix<double> FlowMat=geometry->NormalFlowMaty();
    size_t nvert=mesh->nVertices();
    SparseMatrix<double> resulting_op(nvert,nvert);
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for(size_t index; index<mesh->nVertices();index++)
    {
        tripletList.push_back(T(index,index,1.0) );
    }
    resulting_op.setFromTriplets(tripletList.begin(),tripletList.end());


    resulting_op=(M+h*FlowMat);
    // resulting_op=resulting_op;
    
    return resulting_op; // placeholder
}

SparseMatrix<double> NormalFlow::buildFlowOperatorz(const SparseMatrix<double>& M, double h) const {

    SparseMatrix<double> FlowMat=geometry->NormalFlowMatz();
    size_t nvert=mesh->nVertices();
    SparseMatrix<double> resulting_op(nvert,nvert);
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for(size_t index; index<mesh->nVertices();index++)
    {
        tripletList.push_back(T(index,index,1.0) );
    }
    resulting_op.setFromTriplets(tripletList.begin(),tripletList.end());


    resulting_op=(M+(h/6)*FlowMat);
    // resulting_op=resulting_op;// TODO
    return resulting_op; // placeholder
}




/*
 * Performs mean curvature flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void NormalFlow::integrate(double h) {

    
    
    SparseMatrix<double> M=geometry->massMatrix();
    SparseMatrix<double> Op_flow=buildFlowOperator(M,h);
    
    // SparseMatrix<double> Op_flowx=buildFlowOperatorx(M,h);
    // SparseMatrix<double> Op_flowy=buildFlowOperatory(M,h);
    // SparseMatrix<double> Op_flowz=buildFlowOperatorz(M,h);
    
    size_t nvert=mesh->nVertices();
    
    
    Eigen::SparseLU<SparseMatrix<double> > solver;
    
    // Eigen::SparseLU<SparseMatrix<double> > solver1;
    // Eigen::SparseLU<SparseMatrix<double> > solver2;
    // Eigen::SparseLU<SparseMatrix<double> > solver3;
    

    // Vector<double> x_pos(mesh->nVertices());
    // Vector<double> y_pos(mesh->nVertices());
    // Vector<double> z_pos(mesh->nVertices());
    Vector<double> X_vect(3*mesh->nVertices());
    
    for(Vertex v: mesh-> vertices())
        {
        Vector3 Position =geometry->inputVertexPositions[v.getIndex()];
        X_vect.coeffRef(v.getIndex())=Position.x;
        X_vect.coeffRef(v.getIndex()+nvert)=Position.y;
        X_vect.coeffRef(v.getIndex()+2*nvert)=Position.z;

        // x_pos.coeffRef(v.getIndex())=Position.x;
        // y_pos.coeffRef(v.getIndex())=Position.y;
        // z_pos.coeffRef(v.getIndex())=Position.z;
        }

    
    // Vector<double> result = solver.solve(rhs);
    // solver1.compute(Op_flowx);
    // solver2.compute(Op_flowy);
    // solver3.compute(Op_flowz);
    solver.compute(Op_flow);
    
  
    // Lets try subtracting the mean x position

    // Vector<double> rhs=(x_pos);
   
    // Vector<double> new_x= solver1.solve(rhs);

    // rhs=y_pos;
    // Vector<double> new_y= solver2.solve(rhs);


    // rhs=z_pos;
    // Vector<double> new_z= solver3.solve(rhs);

    Vector<double> rhs=(X_vect);
   
    Vector<double> new_X=solver.solve(rhs);




    // Note: Update positions via geometry->inputVertexPositions
    Vector3 Update={0,0,0};
    Vector3 Update2={0,0,0};
    size_t index_actual;
    for (Vertex v : mesh->vertices()) {
        // Vector3 Update={0,0,0};
        // Update={0,0,0};
        // Vector3 Position =geometry->inputVertexPositions[v.getIndex()];
        
        // for(Halfedge he : v.incomingHalfedges()){
        //     Update+= (1/6.0)*cross(geometry->inputVertexPositions[he.tailVertex()] ,geometry->inputVertexPositions[he.next().tipVertex()] );

        // }   
        index_actual=v.getIndex();
        Update2={new_X[index_actual],new_X[index_actual+nvert],new_X[index_actual+2*nvert]};

        // if(index_actual==2){
        //     // std::cout << "The position was "<< ","<<Position[0]<< ","<<Position[1]<< ","<<Position[2] << ","<<"and now the position will be "<< new_x[index_actual]<< ","<< new_y[index_actual] <<  ","<<new_z[index_actual] <<    "\n \n";
        //     // std::cout<<"THe difference between the exact product and the matrix is:" << (geometry->inputVertexPositions[v]+h*Update-Update2).norm()<<"\n" ;
        //     // std::cout<<"The vectors are "<< geometry->inputVertexPositions[v]+h*Update <<"And the other vector "<< Update2<<"\n";
        //     // std::cout<<"THe difference between the actual and the new vector is "<< geometry->inputVertexPositions[v]-Update2<<"\n";
        // }

        // Update= { new_x[index_actual],new_y[index_actual],new_z[index_actual]  };
        geometry->inputVertexPositions[v] = Update2 ; // placeholder
        geometry->refreshQuantities();
    }

}