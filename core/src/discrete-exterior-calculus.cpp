// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {





// 
// 
// 
// 
// 



/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {
    
    // TODO
    SparseMatrix<double> Hodge_0(mesh.nVertices(),mesh.nVertices());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    Vertex End_condition;
    // double edgelength=0;
    double diag=0.0;

    for(size_t index=0; index< mesh.nVertices();index++){
        Vertex vert=mesh.vertex(index);
        diag=0.0;
        for(Halfedge he : vert.outgoingHalfedges()) {
            Edge Current_edge=he.edge();
            double edgelength2 =(inputVertexPositions[he.tipVertex().getIndex()]-inputVertexPositions[he.tailVertex().getIndex()]).norm2(); 
            // double edgelength=edgeLength(Current_edge);
            // std::cout<<"\t"<<edgelength<<"\t";
            diag= diag+cotan(he)*edgelength2;
            Halfedge he_2=he.twin();
            Current_edge=he_2.edge();
            edgelength2 =(inputVertexPositions[he.tipVertex().getIndex()]-inputVertexPositions[he.tailVertex().getIndex()]).norm2();
            diag= diag+cotan(he_2)*edgelength2;
         }

        
        tripletList.push_back(T(index,index,diag/8.0) );  
    }
    Hodge_0.setFromTriplets(tripletList.begin(),tripletList.end());
    return Hodge_0;// placeholder
}





/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {
    SparseMatrix<double> Hodge_1(mesh.nEdges(),mesh.nEdges());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    Halfedge he;
    double diag=0.0;
    for(size_t index=0;index<mesh.nEdges();index++)
    {

        // For each side i want to calculate this
        diag=0.0;
        Edge Edge_actual=mesh.edge(index);
        he=Edge_actual.halfedge();
        diag+=( cotan(he)+cotan(he.twin()) )/2.0;
        tripletList.push_back(T(index,index,diag) );  
        
    }
    
    Hodge_1.setFromTriplets(tripletList.begin(),tripletList.end());
    // TODO
    return Hodge_1; // placeholder
}


/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {
    SparseMatrix<double> Hodge_2(mesh.nFaces(),mesh.nFaces());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    Halfedge he;
    double diag=0.0;
    double l_ij=0.0;
    double l_jk=0.0;
    double l_ki=0.0;
    double s=0.0;
    for(size_t index=0;index<mesh.nFaces();index++)
    {

        // For each side i want to calculate this
        diag=0.0;
        Face Face_actual=mesh.face(index);
        he=Face_actual.halfedge();
        l_ij=edgeLength(he.edge());
        l_jk=edgeLength(he.next().edge());
        l_ki=edgeLength(he.next().next().edge());
        s=(l_ij+l_jk+l_ki)/2;


        diag= 1/sqrt( s*(s-l_ij)*(s-l_jk)*(s-l_ki))  ;
        tripletList.push_back(T(index,index,diag) );  
    }
    
    Hodge_2.setFromTriplets(tripletList.begin(),tripletList.end());
    


    // TODO
    return Hodge_2; // placeholder
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    // TODO
    // I want a signed adjecency matrix right?
   
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

   
    SparseMatrix<double> d0(mesh.nEdges(),mesh.nVertices());
    size_t n_rows=mesh.nEdges();
    size_t n_cols=mesh.nVertices();
    
    for(Edge E : mesh.edges()) {

        size_t E_index= E.getIndex();
        size_t V_index_tail= E.firstVertex().getIndex();
        size_t V_index_tip=E.secondVertex().getIndex();
        tripletList.push_back(T(E_index,V_index_tail,-1.0) );
        tripletList.push_back(T(E_index,V_index_tip,1.0) );
    }
    d0.setFromTriplets(tripletList.begin(),tripletList.end());
    
    
    
    return d0;
}


/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form3N() const {

    // TODO
    // I want a signed adjecency matrix right?
   
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

    SparseMatrix<double> d0(mesh.nEdges(),3*mesh.nVertices());
    // size_t n_rows=mesh.nEdges();
    size_t n_vert=mesh.nVertices();
    
    for(Edge E : mesh.edges()) {

        size_t E_index= E.getIndex();
        size_t V_index_tail= E.firstVertex().getIndex();
        size_t V_index_tip=E.secondVertex().getIndex();
        tripletList.push_back(T(E_index,V_index_tail,-1.0) );
        tripletList.push_back(T(E_index,V_index_tip,1.0) );

        tripletList.push_back(T(E_index,V_index_tail+n_vert,-1.0) );
        tripletList.push_back(T(E_index,V_index_tip+n_vert,1.0) );

        tripletList.push_back(T(E_index,V_index_tail+2*n_vert,-1.0) );
        tripletList.push_back(T(E_index,V_index_tip+2*n_vert,1.0) );


    }
    d0.setFromTriplets(tripletList.begin(),tripletList.end());
    
    
    
    return d0;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    // TODO

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

   
    SparseMatrix<double> d1(mesh.nFaces(),mesh.nEdges());
    size_t n_rows=mesh.nFaces();
    size_t n_cols=mesh.nEdges();
    SparseMatrix<double> d_0=buildExteriorDerivative0Form();
    
    for(Face f : mesh.faces()){
            // I will start with a halfedge
        size_t face_index= f.getIndex();
        Halfedge he=f.halfedge();
        size_t v1_index=he.tailVertex().getIndex();
        size_t edge_index=he.edge().getIndex();
        for(int iter=0;iter<3;iter++){
        
        if(d_0.coeffRef(edge_index,v1_index)==-1){
            tripletList.push_back(T(face_index,edge_index,1.0) );
            // edge_index,face_index,1
        }
        else{
            tripletList.push_back(T(face_index,edge_index,-1.0) );
        }
        he=he.next();
        v1_index=he.tailVertex().getIndex();
        edge_index=he.edge().getIndex();
        }
        
        
        }
        d1.setFromTriplets(tripletList.begin(),tripletList.end());




    return d1; // placeholder
}

} // namespace surface
} // namespace geometrycentral