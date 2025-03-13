// PLEASE READ:
//
// This file implements additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because
// we are "inside" the class, we no longer have to call
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
// Functions in this file can be called from other projects simply by using geometry->cotan(he),
// geometry->barycentricDualArea(v), etc. where "geometry" is a pointer to a VertexPositionGeometry. This avoids having
// to declare a GeometryRoutines object in every project, and also mimics the way that geometry routines are normally
// called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.


#include "geometrycentral/surface/vertex_position_geometry.h"
#include <complex>
#include "geometrycentral/surface/remeshing.h"

namespace geometrycentral {

namespace surface {
/*
 * Compute the Euler characteristic of the mesh.
 */
int VertexPositionGeometry::eulerCharacteristic() const {
    return (int)mesh.nVertices() - (int)mesh.nEdges() + (int)mesh.nFaces();
}

/*
 * Compute the mean length of all the edges in the mesh.
 *
 * Input:
 * Returns: The mean edge length.
 */
double VertexPositionGeometry::meanEdgeLength() const {

    double total = 0.0;
    for (Edge e : mesh.edges()) {
        total += edgeLength(e);
    }
    return total / mesh.nEdges();



}






/*
 * Compute the total surface area of the mesh.
 *
 * Input:
 * Returns: The surface area of the mesh.
 */
double VertexPositionGeometry::totalArea() const {

    double total = 0.0;
    // std::cout<<"The number of faces is " << mesh.nFaces() << " \n";
    int counter = 0;
    for (Face f : mesh.faces()) {
        if(std::isnan(faceArea(f))) continue;
        
        // if(f.isBoundaryLoop()) std::cout<<"The face is a boundary loop?\n";
        counter+=1;
        total += faceArea(f);
        // std::cout<< faceArea(f) << " ";
    }
    // std::cout<<"The number of considered faces is " << counter <<" \n";
    // std::cout<<" \n";
    // std::cout<<"The total area is " << total <<" \n";
    return total;
}

double VertexPositionGeometry::totalVolume() const{
    double total=0.0;
    for(Face f : mesh.faces()){
        Halfedge he = f.halfedge();
        const Vector3 fi = inputVertexPositions[he.tailVertex()] ;
        const Vector3 fj = inputVertexPositions[he.tipVertex()] ;
        const Vector3 fk = inputVertexPositions[he.next().tipVertex()] ;
        total+= (1/6.0)*dot(fi,cross(fj,fk));

    }
    return total;
}



/*
 * Computes the cotangent of the angle opposite to a halfedge. (Do NOT use built-in function for this)
 *
 * Input: The halfedge whose cotan weight is to be computed.
 * Returns: The cotan of the angle opposite the given halfedge.
 */
double VertexPositionGeometry::cotan(Halfedge he) const {
    //  I have a halfedge 
    
    Vector3 Position_1=inputVertexPositions[he.tipVertex()];
    Vector3 Position_2=inputVertexPositions[he.next().tipVertex()];
    Vector3 Position_3=inputVertexPositions[he.tailVertex()];
    
    Vector3 u=Position_1-Position_2;
    Vector3 v=Position_3-Position_2;

    double cotan = dot(v,u)/norm(cross(v,u));
    
    
    return cotan/2; // placeholder
}

/*
 * Computes the barycentric dual area of a vertex.
 *
 * Input: The vertex whose barycentric dual area is to be computed.
 * Returns: The barycentric dual area of the given vertex.
 */
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {

    // TODO
    // Ok i want to calculate the barycentric dual area which for a vertex is 1/3 or the total of the areas around
    
    double Dual_area=0.0;  
    for(Halfedge he : v.outgoingHalfedges()){
        if(!(he.isInterior())){
            continue;
        }
        Dual_area+=faceArea(he.face());
    }
    
    // for(Face f : v.adjacentFaces()) {
    //     Dual_area+=faceArea(f);

    // }
    
    return Dual_area/3.0; // placeholder
}

/*
 * Computes the angle (in radians) at a given corner. (Do NOT use built-in function for this)
 *
 *
 * Input: The corner at which the angle needs to be computed.
 * Returns: The angle clamped between 0 and π.
 */
double VertexPositionGeometry::angle(Corner c) const {
    Halfedge he_1= c.halfedge();
    Halfedge he_2=he_1.next().next();

    // std::cout <<he_2.tipVertex().getIndex() <<"\t";
    // std::cout <<he_1.tailVertex().getIndex();
    
    Vector3 Point_1=inputVertexPositions[he_1.tailVertex().getIndex()];
    Vector3 Point_2=inputVertexPositions[he_1.tipVertex().getIndex()];
    Vector3 Point_3=inputVertexPositions[he_2.tailVertex().getIndex()];
    // acos()
    // Halfedge 1 is pointing outside the corner and he_2 inside the corner
    Vector3 Edge_1=Point_3-Point_1;
    Vector3 Edge_2=Point_2-Point_1;

    double result_angle= acos(dot(Edge_1,Edge_2)/(norm(Edge_1)*norm(Edge_2)));



    // TODO
    return result_angle; // placeholder
}

/*
 * Computes the signed angle (in radians) between two adjacent faces. (Do NOT use built-in function for this)
 *
 * Input: The halfedge (shared by the two adjacent faces) on which the dihedral angle is computed.
 * Returns: The dihedral angle.
 */
double VertexPositionGeometry::dihedralAngle(Halfedge he) const {
    // Ok we will start with
    Vector3 Normal_1 = faceNormal(he.face());
    Vector3 Normal_2 = faceNormal(he.twin().face());

    // i HAVE THE TWO NORMALS 
    Vector3 e_ij= inputVertexPositions[he.tipVertex().getIndex()]-inputVertexPositions[he.tailVertex().getIndex()];
    e_ij=e_ij/norm(e_ij);
    double theta_ij=atan2( dot(e_ij, cross( Normal_1,Normal_2 )),dot(Normal_1,Normal_2)  );


    // TODO
    return theta_ij; // placeholder
}

/*
 * Computes the normal at a vertex using the "equally weighted" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "equally weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalEquallyWeighted(Vertex v) const {
    Vector3 Normal={0,0,0};
    for(Face f : v.adjacentFaces()) {
  // do science here
    Normal+=faceNormal(f); 
    }
    Normal=Normal/norm(Normal);
    // TODO
    return Normal; // placeholder
}

/*
 * Computes the normal at a vertex using the "tip angle weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "tip angle weights" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAngleWeighted(Vertex v) const {
    Vector3 Normal={0,0,0};
    // 
    Corner corner;
    // how can i find the corner 

    for(Halfedge he : v.outgoingHalfedges()){
        if(!he.isInterior()){
            continue;
        }
        corner=he.corner();
        Normal+=faceNormal(he.face())*angle(corner);
        // Normal+=faceNormal(he.face());
    }



//     for(Face f : v.adjacentFaces()) {
//   // do science here
//     for(Halfedge he : f.adjacentHalfedges()) 
//     {
//         if(he.tailVertex()==v){
//             corner=he.corner();
//         }
    
//     }
    
//     Normal+=faceNormal(f)*angle(corner); 
    
    
//     }
    // Normal=Normal;
    // /norm(Normal);
    // TODO
    // TODO
    return Normal; // placeholder
}

/*
 * Computes the normal at a vertex using the "inscribed sphere" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "inscribed sphere" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalSphereInscribed(Vertex v) const {

    Vector3 Normal={0,0,0};
    // 
    Vector3 Point_1=inputVertexPositions[v];
    Vector3 Point_2;
    Vector3 Point_3;
    Vector3 e_ij;
    Vector3 e_ik;

    // how can i find the corner 
    for(Face f : v.adjacentFaces()) {
  // do science here
    for(Halfedge he : f.adjacentHalfedges()) 
    {
        if(he.tailVertex()==v){
        Point_2=inputVertexPositions[he.tipVertex().getIndex()];        
        }
        if(he.tipVertex()==v){
        Point_3=inputVertexPositions[he.tailVertex().getIndex()];
        }
    
    }
    
    e_ij=Point_2-Point_1;
    e_ik=Point_3-Point_1;

    Normal+=cross(e_ij,e_ik)/(norm(e_ij)*norm(e_ik)*norm(e_ij)*norm(e_ik)); 
    
    
    }
    Normal=Normal/norm(Normal);

    // TODO
    return Normal; // placeholder
}

/*
 * Computes the normal at a vertex using the "face area weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "face area weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAreaWeighted(Vertex v) const {
    Vector3 Normal={0,0,0};
    for(Face f : v.adjacentFaces()) {
  // do science here
    Normal+=faceNormal(f)*faceArea(f); 
    }
    Normal=Normal/norm(Normal);
    // TODO
    return Normal; // placeholder
}

/*
 * Computes the normal at a vertex using the "Gauss curvature" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "Gauss curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalGaussianCurvature(Vertex v) const {

    Vector3 Position_1=inputVertexPositions[v];
    Vector3 Position_2;
    Vector3 Normal={0,0,0};
    Vector3 e_ij;
    for(Halfedge he: v.incomingHalfedges()){
        Position_2=inputVertexPositions[he.tailVertex().getIndex()];
        e_ij=Position_2-Position_1;
        Normal+=dihedralAngle(he)*e_ij/norm(e_ij);

    }
    // for(Halfedge he: v.outgoingHalfedges()){
    //     Position_2=inputVertexPositions[he.tipVertex().getIndex()];
    //     e_ij=Position_2-Position_1;
    //     Normal+=dihedralAngle(he)*e_ij/norm(e_ij);

    // }
    Normal=Normal/norm(Normal);
    return Normal;

}

/*
 * Computes the normal at a vertex using the "mean curvature" method (equivalent to the "area gradient" method).
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "mean curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalMeanCurvature(Vertex v) const {
    Vector3 Position_1=inputVertexPositions[v];
    Vector3 Position_2;
    Vector3 Normal={0,0,0};
    Vector3 e_ji;
    for(Halfedge he: v.incomingHalfedges()){
        Position_2=inputVertexPositions[he.tailVertex().getIndex()];
        e_ji=Position_1-Position_2;
        
        Normal+=0.5*(cotan(he.twin())+cotan(he))*e_ji;

    }

    // Normal=Normal/norm(Normal);
    return Normal;


    return {0, 0, 0}; // placeholder
}





/*
 * Computes the angle defect at a vertex.
 *
 * Input: The vertex whose angle defect is to be computed.
 * Returns: The angle defect of the given vertex.
 */
double VertexPositionGeometry::angleDefect(Vertex v) const {
    double d =2*PI;
    Corner corner;
    for(Face f : v.adjacentFaces()) {
    for(Halfedge he : f.adjacentHalfedges()) 
    {
        if(he.tailVertex()==v){
            corner=he.corner();
        }
    }
    
    d-=angle(corner); 
    
    
    }

    // TODO
    return d; // placeholder
}

/*
 * Computes the total angle defect of the mesh.
 *
 * Input:
 * Returns: The total angle defect
 */
double VertexPositionGeometry::totalAngleDefect() const {
    double total_d=0;
    for(Vertex v : mesh.vertices()) {
    total_d+=angleDefect(v);
    }
    

    // TODO
    return total_d; // placeholder
}


/*
Compute the gradient of an angle at vertex i with respect to vertex j 

*/
Vector3 VertexPositionGeometry::Angle_grad(Corner c, Vertex j, Vector3 Normal){

    Vector3 u;
    Vector3 v;
    Vector3 grad;
    if(c.vertex().getIndex()==j.getIndex()){
        u = inputVertexPositions[c.halfedge().next().vertex()] -inputVertexPositions[c.halfedge().vertex()];
        v = inputVertexPositions[c.halfedge().next().next().vertex()] -inputVertexPositions[c.halfedge().vertex()];
        grad =  ( cross(Normal,u)/u.norm2() -cross(Normal,v)/v.norm2());
        return -1*grad;
    }

    if(c.halfedge().next().vertex().getIndex()==j.getIndex()){
        u=inputVertexPositions[c.halfedge().next().vertex()] -inputVertexPositions[c.halfedge().vertex()];
        grad = -1*(cross(Normal,u)/u.norm2()); 
        return -1*grad;
    }
    
    if(c.halfedge().next().next().vertex().getIndex()==j.getIndex()){
        v = inputVertexPositions[c.halfedge().next().next().vertex()] -inputVertexPositions[c.halfedge().vertex()];
        grad=cross(Normal,v)/v.norm2();
        return -1*grad;
    }

    std::cout<< "The vertex is not in the triangle where the angle is \n";

return Vector3({0,0,0});
}





/*
 * Computes the (integrated) scalar mean curvature at a vertex.
 *
 * Input: The vertex whose mean curvature is to be computed.
 * Returns: The mean curvature at the given vertex.
 */
double VertexPositionGeometry::scalarMeanCurvature(Vertex v) const {

    Vector3 Position_1=inputVertexPositions[v];
    Vector3 Position_2;
    double scalar_mean_curve=0;
    Vector3 e_ij;
    for(Halfedge he: v.incomingHalfedges()){
        Position_2=inputVertexPositions[he.tailVertex().getIndex()];
        e_ij=Position_2-Position_1;
        scalar_mean_curve+=dihedralAngle(he)*norm(e_ij);

    }
    

    return 0.25*scalar_mean_curve; // placeholder
}

/*
 * Computes the circumcentric dual area of a vertex.
 *
 * Input: The vertex whose circumcentric dual area is to be computed.
 * Returns: The circumcentric dual area of the given vertex.
 */
double VertexPositionGeometry::circumcentricDualArea(Vertex v) const {
    double Area=0;

    Vector3 Normal={0,0,0};
    Vector3 Point_1=inputVertexPositions[v];
    Vector3 Point_2;
    Vector3 e_ij;
    Corner corner;
    for(Face f : v.adjacentFaces()) {
    for(Halfedge he : f.adjacentHalfedges()) 
    {
        if(he.tailVertex()==v){
            Point_2=inputVertexPositions[he.tipVertex()];
            e_ij=Point_2-Point_1;
            Area+=cotan(he)*norm(e_ij)*norm(e_ij);
           
        }
        if(he.tipVertex()==v){
            Point_2=inputVertexPositions[he.tailVertex()];
            e_ij=Point_2-Point_1;
            Area+=cotan(he)*norm(e_ij)*norm(e_ij);
           
        }
        


    }
    
    
    }
    

    // TODO
    return Area/4; // placeholder
}

/*
 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
 *
 * Input: The vertex on which the principal curvatures need to be computed.
 * Returns: A std::pair containing the minimum and maximum principal curvature values at a vertex.
 */
std::pair<double, double> VertexPositionGeometry::principalCurvatures(Vertex v) const {

    // TODO
    double H;
    double K;

    H=scalarMeanCurvature(v)/circumcentricDualArea(v);
    K=angleDefect(v)/circumcentricDualArea(v);


    return std::make_pair(H-sqrt(H*H-K), H+sqrt(H*H-K)); // placeholder
}


/*
 * Builds the sparse POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace matrix,
 * multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse positive definite Laplace matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrix() const {
    SparseMatrix<double> d0=this->buildExteriorDerivative0Form();
    // SparseMatrix<double> H1=this->buildHodgeStar1Form() ;
    // SparseMatrix<double> d1=this->buildExteriorDerivative1Form();
    // SparseMatrix<double> H2=this->buildHodgeStar2Form();
   
    SparseMatrix<double> laplace=(d0.transpose()*(buildHodgeStar1Form()*d0));
    SparseMatrix<double> diagonal(mesh.nVertices(),mesh.nVertices());

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    double eps = 1e-8;
    for(size_t index = 0; index<mesh.nVertices();index++)
    {
        tripletList.push_back(T(index,index,eps) );
    }
    diagonal.setFromTriplets(tripletList.begin(),tripletList.end());

    // std::cout << "we are here right?";
    return laplace+diagonal; // placeholder
}

SparseMatrix<double> VertexPositionGeometry::uniformlaplacianMatrix() const{
    SparseMatrix<double> L(mesh.nVertices(),mesh.nVertices());

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

    int id1;
    int  id2;
    double neigh;
    for(Vertex v : mesh.vertices()){
        id1 = v.getIndex();
        neigh = 1/v.degree();
        for( Halfedge he : v.outgoingHalfedges()){
        id2 = he.next().vertex().getIndex();
        tripletList.push_back(T(id1,id1,neigh) );
        tripletList.push_back(T(id1,id2,-1*neigh) );
        
        // tripletList.push_back(id1,id1,neigh);
        // tripletList.push_back(id1,id2,-1*neigh);

        }
    L.setFromTriplets(tripletList.begin(),tripletList.end());
    return L;
    }


}



SparseMatrix<double> VertexPositionGeometry::laplaceMatrix3d() const {
    size_t nvert=mesh.nVertices();
    
    // SparseMatrix<double> d0=this->buildExteriorDerivative0Form3N();
    SparseMatrix<double> laplace= laplaceMatrix();
    // SparseMatrix<double> diagonal(3*mesh.nVertices(),3*mesh.nVertices());

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    double eps = 1e-8;
    SparseMatrix<double> laplace3d(3*nvert,3*nvert);
    int row;
    int col;
    double value;
    for (int k=0; k<laplace.outerSize();++k){
        for( SparseMatrix<double>::InnerIterator it(laplace,k); it; ++it){
            value=it.value();
            row=it.row();
            col=it.col();
            // it.index();
            tripletList.push_back(T(row,col,value));
            tripletList.push_back(T(row+nvert,col+nvert,value));
            tripletList.push_back(T(row+2*nvert,col+2*nvert,value));
                        
        }

    }
    laplace3d.setFromTriplets(tripletList.begin(),tripletList.end());


    // for(size_t index = 0; index<3*nvert;index++)
    // {
    //     tripletList.push_back(T(index,index,eps) );
    // }
    // diagonal.setFromTriplets(tripletList.begin(),tripletList.end());

    // std::cout << "we are here right?";
    return laplace3d; // placeholder
}




SparseMatrix<double> VertexPositionGeometry::laplaceMatrix2() const {
    size_t N_vert=mesh.nVertices();
    SparseMatrix<double> laplace(N_vert,N_vert);


    // The idea would be to create a sparse matrix and fill it with triplets


    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for(size_t index = 0; index<N_vert;index++)
    {
        Vertex v_i=mesh.vertex(index);
        double val_off_diag;
        double val_diag=0;  
        for(Halfedge he : v_i.outgoingHalfedges()) {
            size_t v_f_index=he.tipVertex().getIndex();
            val_off_diag= -1*0.5*(cotan(he)+cotan(he.twin()) );
            tripletList.push_back(T(index,v_f_index,val_off_diag) );    
            val_diag+=-1*val_off_diag;
        }
        val_diag+=1e-8;
        tripletList.push_back(T(index,index,val_diag) );
    }
    laplace.setFromTriplets(tripletList.begin(),tripletList.end());

    // std::cout << "we are here right?";
    return laplace; // placeholder
}


SparseMatrix<double> VertexPositionGeometry::NormalFlowMat() const {
    size_t N_vert=mesh.nVertices();
    SparseMatrix<double> Flow_op(3*N_vert,3*N_vert);

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for(size_t index = 0; index<N_vert;index++)
    {
        Vertex v_i=mesh.vertex(index);
    
        for(Halfedge he : v_i.incomingHalfedges()) {
            size_t p_index=he.next().tipVertex().getIndex();
            if(p_index+2*N_vert==3*N_vert){
                std::cout <<"This is clearly wrong\n";
            }
            // I need to decide which vertex am i saving in th ematrix
            // Lets use this first one ok...
            Vector3 fj=inputVertexPositions[he.tailVertex()];
            // THis is the x component of vertex index 
            tripletList.push_back(T(index,p_index+2*N_vert,fj.y/6) );
            tripletList.push_back(T(index,p_index+N_vert,-1*fj.z/6) );

            tripletList.push_back(T(index+N_vert,p_index,fj.z/6) );
            tripletList.push_back(T(index+N_vert,p_index+2*N_vert,-1*fj.x/6) );

            tripletList.push_back(T(index+2*N_vert,p_index+N_vert,fj.x/6) );
            tripletList.push_back(T(index+2*N_vert,p_index,-1*fj.y/6) );
            

        }
    }
    Flow_op.setFromTriplets(tripletList.begin(),tripletList.end());

    // std::cout << "we are here right?";
    return Flow_op; // placeholder
}


SparseMatrix<double> VertexPositionGeometry::NormalFlowMatx() const {
    size_t N_vert=mesh.nVertices();
    SparseMatrix<double> Flow_op(N_vert,N_vert);

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for(size_t index = 0; index<N_vert;index++)
    {
        Vertex v_i=mesh.vertex(index);
        
        tripletList.push_back(T(index,index,vertexNormalSphereInscribed(v_i).x ) );
    }
    Flow_op.setFromTriplets(tripletList.begin(),tripletList.end());

    return Flow_op; 
}

SparseMatrix<double> VertexPositionGeometry::NormalFlowMaty() const {
    size_t N_vert=mesh.nVertices();
    SparseMatrix<double> Flow_op(N_vert,N_vert);

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for(size_t index = 0; index<N_vert;index++)
    {
        Vertex v_i=mesh.vertex(index);
        tripletList.push_back(T(index,index,vertexNormalSphereInscribed(v_i).y) );
    }
    Flow_op.setFromTriplets(tripletList.begin(),tripletList.end());

    return Flow_op; 
}

SparseMatrix<double> VertexPositionGeometry::NormalFlowMatz() const {
    size_t N_vert=mesh.nVertices();
    SparseMatrix<double> Flow_op(N_vert,N_vert);

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for(size_t index = 0; index<N_vert;index++)
    {
        Vertex v_i=mesh.vertex(index);
        tripletList.push_back(T(index,index,vertexNormalSphereInscribed(v_i).z ) );
    }
    Flow_op.setFromTriplets(tripletList.begin(),tripletList.end());

    return Flow_op; 
}


SparseMatrix<double> VertexPositionGeometry::GaussFlowMat() const {
    size_t N_vert=mesh.nVertices();
    SparseMatrix<double> Gauss(N_vert,N_vert);
    
    // The idea would be to create a sparse matrix and fill it with triplets

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for(size_t index = 0; index<N_vert;index++)
    {
        Vertex v_i=mesh.vertex(index);
        double val_off_diag;
        double val_diag=0;  
        for(Halfedge he : v_i.outgoingHalfedges()) {
            Vertex v_f=he.tipVertex();
            size_t v_f_index=v_f.getIndex();
            // val_off_diag=-0.5*
            val_off_diag= dihedralAngle(he)/(inputVertexPositions[v_i]-inputVertexPositions[v_f]).norm();
            
            tripletList.push_back(T(index,v_f_index,val_off_diag) );    
            val_diag+=-1*val_off_diag;
        }
        val_diag+=1e-8;
        tripletList.push_back(T(index,index,val_diag) );
    }
    Gauss.setFromTriplets(tripletList.begin(),tripletList.end());

    // std::cout << "we are here right?";
    return Gauss; // placeholder
}

SparseMatrix<double> VertexPositionGeometry::GaussFlowMat3d(VertexData<double> Scalar_MC) const {
    size_t N_vert=mesh.nVertices();
    SparseMatrix<double> Gauss(3*N_vert,3*N_vert);

    // The idea would be to create a sparse matrix and fill it with triplets

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for(size_t index = 0; index<N_vert;index++)
    {
        Vertex v_i=mesh.vertex(index);
        double val_off_diag;
        double val_diag=0;  
        double prefactor;
        for(Halfedge he : v_i.outgoingHalfedges()) {
            Vertex v_f=he.tipVertex();
            size_t v_f_index=v_f.getIndex();
            // val_off_diag=-0.5*
            prefactor=-1*(Scalar_MC[index] +Scalar_MC[v_f_index]);
            val_off_diag= prefactor*-1*0.5*dihedralAngle(he)/(inputVertexPositions[v_i]-inputVertexPositions[v_f]).norm();
            
            tripletList.push_back(T(index,v_f_index,val_off_diag) );
            tripletList.push_back(T(index+N_vert,v_f_index+N_vert,val_off_diag) );
            tripletList.push_back(T(index+2*N_vert,v_f_index+2*N_vert,val_off_diag) );
                
            val_diag+=-1*val_off_diag;
        }
        val_diag+=1e-8;
        tripletList.push_back(T(index,index,val_diag) );
        tripletList.push_back(T(index+N_vert,index+N_vert,val_diag) );
        tripletList.push_back(T(index+2*N_vert,index+2*N_vert,val_diag) );
    
    }
    Gauss.setFromTriplets(tripletList.begin(),tripletList.end());

    return Gauss; // placeholder
}





SparseMatrix<double> VertexPositionGeometry::MeanHFlowMat3d(VertexData<double> Scalar_MC,Vector<Vector3> face_Normals) const {
    size_t N_vert=mesh.nVertices();
    SparseMatrix<double> MeanH(3*N_vert,3*N_vert);

    // The idea would be to create a sparse matrix and fill it with triplets
    Vector3 Normal1;
    Vector3 Normal2;
    Vertex v_1;
    Vertex v_2;
    Vertex v_3;
    size_t v_1_index;
    size_t v_2_index;
    size_t v_3_index;

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for(size_t index = 0; index<N_vert;index++)
    {
        Vertex v_i=mesh.vertex(index);
        double val_off_diag;
        double val_diag=0;  
        double prefactor;
        for(Halfedge he : v_i.outgoingHalfedges()) {
            v_1=he.tipVertex();
            v_2=he.next().tipVertex();
            v_3=he.twin().next().tipVertex();
            v_1_index=v_1.getIndex();
            v_2_index=v_2.getIndex();
            v_3_index=v_3.getIndex();
            Normal1=face_Normals[he.face().getIndex()];
            Normal2=face_Normals[he.twin().face().getIndex()];
            // Normal1={0,0,0};
            // Normal2={0,0,0};
            
            // I need to do the prefactor again
            prefactor=(1/3.0)*(Scalar_MC[index]*Scalar_MC[index])+(2.0/3.0)*(Scalar_MC[v_1_index]*Scalar_MC[v_1_index]);
            // prefactor=1.0;
            Normal1=prefactor*Normal1/4; //The factor of 4 appear in the cross product
            Normal2=prefactor*Normal2/4;
            // I added it to the normals since it will be used in all the next multiplications and the cross product is a linear operator


            // This is N1XV2
            tripletList.push_back(T(index,v_2_index+2*N_vert,Normal1.y) ); //n1y v2z
            tripletList.push_back(T(index,v_2_index+N_vert,-1*Normal1.z) ); // - n1z v2y
            
            tripletList.push_back(T(index+N_vert,v_2_index,Normal1.z) ); //n1z v2x
            tripletList.push_back(T(index+N_vert,v_2_index+2*N_vert,-1*Normal1.x) ); // - n1x v2z
            
            tripletList.push_back(T(index+2*N_vert,v_2_index+N_vert,Normal1.x) ); //n1x v2y
            tripletList.push_back(T(index+2*N_vert,v_2_index,-1*Normal1.y) ); // - n1y v2x
            
            // tHIS IS N2xV3
            tripletList.push_back(T(index,v_3_index+2*N_vert,-1*Normal2.y) ); //n2y v3z
            tripletList.push_back(T(index,v_3_index+N_vert,Normal2.z) ); // - n2z v3y
            
            tripletList.push_back(T(index+N_vert,v_3_index,-1*Normal2.z) ); //n2z v3x
            tripletList.push_back(T(index+N_vert,v_3_index+2*N_vert,Normal2.x) ); // - n2x v3z
            
            tripletList.push_back(T(index+2*N_vert,v_3_index+N_vert,-1*Normal2.x) ); //n2x v3y
            tripletList.push_back(T(index+2*N_vert,v_3_index,Normal2.y) ); // - n2y v3x
            

            // This is (N2-N1)xV1
            tripletList.push_back(T(index,v_1_index+2*N_vert,(Normal2.y-Normal1.y)) ); //n1y v2z
            tripletList.push_back(T(index,v_1_index+N_vert,-1*(Normal2.z-Normal1.z)) ); // - n1z v2y
            
            tripletList.push_back(T(index+N_vert,v_1_index,(Normal2.z-Normal1.z)) ); //n1z v2x
            tripletList.push_back(T(index+N_vert,v_1_index+2*N_vert,-1*(Normal2.x-Normal1.x)) ); // - n1x v2z
            
            tripletList.push_back(T(index+2*N_vert,v_1_index+N_vert,(Normal2.x-Normal1.x)) ); //n1x v2y
            tripletList.push_back(T(index+2*N_vert,v_1_index,-1*(Normal2.y-Normal1.y)) ); // - n1y v2x
            
        }
    
    }
    MeanH.setFromTriplets(tripletList.begin(),tripletList.end());

    return MeanH; // placeholder
}


Vector3 VertexPositionGeometry::dihedralAngleGradient(Halfedge he, Vertex v,Vector<Vector3> face_Normals) const {
    // std::cout<< he.edge().isBoundary();

    double l = edgeLength(he.edge());



    if (he.edge().isBoundary()) {
    return Vector3{0, 0, 0};
  } else if (he.vertex() == v) { //This is only used for the SIJ_1
    return (cotan(he.next().next()) *
                face_Normals[he.face().getIndex()] +
            cotan(he.twin().next()) *
                face_Normals[he.twin().face().getIndex()]) /l;
  } else if (he.next().vertex() == v) { //This is for the firt s term
    return (cotan(he.twin().next().next()) *
                face_Normals[he.twin().face().getIndex()] +
            cotan(he.next()) *
                face_Normals[he.face().getIndex()]) /
           l;
  } else if (he.next().next().vertex() == v) { //Este ocurre para el segundo termino
    return (-(cotan(he.next().next()) +
              cotan(he.next())) *
            face_Normals[he.face().getIndex()]) /
           l;
  } else {
    // mem3dg_runtime_error("Unexpected combination of halfedge and vertex!");
    std::cout<< "THe dihedral angle gradient is not working\n";
    return Vector3{0, 0, 0};
  }
std::cout<< "THis is impossible to print\n";
return Vector3{0, 0, 0};
}





SparseMatrix<double> VertexPositionGeometry::Schlafi3d(VertexData<double> Scalar_MC,Vector<Vector3> face_Normals) const {
    size_t N_vert=mesh.nVertices();
    SparseMatrix<double> Schlafi(3*N_vert,3*N_vert);

    // The idea would be to create a sparse matrix and fill it with triplets

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    Vector3 Pos_i;
    Vector3 Sij_1;
    Vector3 Sij_2;
    for(size_t index = 0; index<N_vert;index++)
    {
        Vertex v_i=mesh.vertex(index);
        double val_off_diag;
        double val_diag=0;  
        double prefactor;
        Sij_1={0,0,0};
        Sij_2={0,0,0};
        for(Halfedge he : v_i.outgoingHalfedges()) {
            Vertex v_f=he.tipVertex();
            size_t v_f_index=v_f.getIndex();
            // val_off_diag=-0.5*
            prefactor=-1*(Scalar_MC[index]);
            Sij_1+=prefactor*edgeLength(he.edge())*dihedralAngleGradient(he,he.vertex(),face_Normals);
            prefactor=-1*(Scalar_MC[v_f_index]);
            Sij_2+=prefactor*( edgeLength(he.twin().edge())*dihedralAngleGradient(he.twin(),he.vertex(),face_Normals) +  edgeLength(he.next().edge())*dihedralAngleGradient(he.next(),he.vertex(),face_Normals) +  edgeLength(he.twin().next().next().edge())*dihedralAngleGradient(he.twin().next().next(), he.vertex(),face_Normals)); 


                
            // val_diag+=-1*val_off_diag;
        }
        Pos_i=inputVertexPositions[v_i];

        tripletList.push_back(T(index,index,(Sij_1.x+Sij_2.x)/Pos_i.x) );
        tripletList.push_back(T(index+N_vert,index+N_vert,(Sij_1.y+Sij_2.y)/Pos_i.y ) );
        tripletList.push_back(T(index+2*N_vert,index+2*N_vert,(Sij_1.z+Sij_2.z)/Pos_i.z) );
    
    }
    Schlafi.setFromTriplets(tripletList.begin(),tripletList.end());

    return Schlafi; // placeholder
}



/*
 * Builds the sparse  matrix containing the local Q for an interior edge.
 *
 * Input:
 * Returns: Local Q for willmore flow
 */

//         x2
//         /\
//        /  \
//     e1/    \e3
//      /  t0  \
//     /        \
//    /    e0    \
//  x0------------x1
//    \          /
//     \   t1   /
//      \      /
//     e2\    /e4
//        \  /
//         \/
//         x3
//
// Edge orientation: e0,e1,e2 point away from x0
//                      e3,e4 point away from x1


SparseMatrix<double> VertexPositionGeometry::LocalQWillmore(const Halfedge he, const Vector3 x0, const Vector3 x1, const Vector3 x2, const Vector3 x3) const{
    // Ok the thing we want to build is 
    SparseMatrix<double> Q(4,4);
    

    
    const Vector3 e0 = x1-x0;
    const Vector3 e1 = x2-x0;
    const Vector3 e2 = x3-x0;
    const Vector3 e3 = x2-x1;
    const Vector3 e4 = x3-x1;

    const double c01 = 2.0*cotan(he.next());
    const double c02 = 2.0*cotan(he.twin().next().next());
    const double c03 = 2.0*cotan(he.next().next());
    const double c04 = 2.0*cotan(he.twin().next());

    const double K0[] = {c03+c04, c01+c02, -c01-c03, -c02-c04};

    const double A0 = 0.5 * cross(e0,e1).norm();
    const double A1 = 0.5 * cross(e0,e2).norm();

    const double coef = -3. / ((A0+A1));

    assert(finite(coef));
    assert(finite(c01));
    assert(finite(c02));
    assert(finite(c03));
    assert(finite(c04));


    for (int i=0; i<4; ++i) {
    for (int j=0; j<i; ++j) {
        // tripletList.push_back(T(i,j,coef*K0[i] * K0[j]));
        // tripletList.push_back(T(j,i,coef*K0[i] * K0[j]));
        Q.insert(i,j) = coef*K0[i] * K0[j];
        Q.insert(j,i) = coef*K0[i] * K0[j];

    //   Q[i,j] = Q[j,i] = coef * K0[i] * K0[j];
    }
    Q.insert(i,i)= coef * K0[i] * K0[i];
    }



    return Q;
}

// I may change this so the function is a mutator instead of a 


SparseMatrix<double> VertexPositionGeometry::FullQWillmore() const{

    SparseMatrix<double> Q(mesh.nVertices(),mesh.nVertices());
    // I need the indices of the vertics
    SparseMatrix<double> LocalQ(4,4);
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for(size_t ei=0; ei< mesh.nEdges(); ++ei){
        Edge e = mesh.edge(ei);
        // Here we have a problem you see because when i call this halfedge the other function may call the other
        if( e.isBoundary() ) continue;
        Halfedge he = e.halfedge();
        int fromVertexIdx= he.tailVertex().getIndex();
        int toVertexIdx = he.tipVertex().getIndex();
        int leftApexIdx = he.next().tipVertex().getIndex();
        int rightApexIdx = he.twin().next().tipVertex().getIndex();
        int idx[] = {fromVertexIdx, // get index to tail of edge
                 toVertexIdx,   // get index to head of edge
                 leftApexIdx,   // opposite vertex on left-side  triangle
                 rightApexIdx}; // opposite vertex on right-side triangle

        LocalQ = LocalQWillmore(he,inputVertexPositions[idx[0]],inputVertexPositions[idx[1]],inputVertexPositions[idx[2]],inputVertexPositions[idx[3]] );
        for( int i=0; i<4;++i){
            for(int j=0;j<4;++j){
                tripletList.push_back(T(idx[i],idx[j],LocalQ.coeffRef(i,j)));
       
            }
        }
        // wE WILL USE TRIPLETS BECAUSE THIS MATRIX IS SPARSE.
    }



    Q.setFromTriplets(tripletList.begin(),tripletList.end());

    return Q;
}



SparseMatrix<double> VertexPositionGeometry::FullKWillmore() const {

    // Great now i need to do this
    SparseMatrix<double> K(mesh.nVertices(),mesh.nVertices());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    // I need to iterate over every edge.
    Halfedge he;
    Vertex vi;
    Vertex vj;
    Vertex vk;
    Vertex vl;
    Vector3 a;
    Vector3 b;
    Vector3 c;
    Vector3 d;
    double cos_beta;
    double a_norm;
    double b_norm;
    double c_norm;
    double d_norm;
    size_t V_index;


    
    for(Vertex V : mesh.vertices()){
        // For each vertex i need to iterate
        V_index=V.getIndex();
        for(Halfedge he : V.incomingHalfedges()) {
        
        vi=he.tipVertex();
        vj=he.tailVertex();
        vk=he.next().tipVertex();
        vl=he.twin().next().tipVertex();
        
        a=inputVertexPositions[vj]-inputVertexPositions[vk];
        b=inputVertexPositions[vl]-inputVertexPositions[vj];
        c=inputVertexPositions[vi]-inputVertexPositions[vl];
        d=inputVertexPositions[vk]-inputVertexPositions[vi];
        a_norm=a.norm();
        b_norm=b.norm();
        c_norm=c.norm();
        d_norm=d.norm();
        cos_beta= (dot(a,c)*dot(b,d) -dot(a,b)*dot(c,d)-dot(b,c)-dot(d,a))/(a_norm*b_norm*c_norm*d_norm);

        tripletList.push_back(T(  V_index , vj.getIndex()  , cos_beta/(a_norm*a_norm) -dot(b,c)/(a_norm,b_norm,c_norm,d_norm)  ));
        tripletList.push_back(T(  V_index , vk.getIndex()  , -1.0*cos_beta/(a_norm*a_norm) +dot(b,c)/(a_norm,b_norm,c_norm,d_norm)  ));
        
        tripletList.push_back(T(  V_index , vl.getIndex()  , dot(a,c)/(a_norm*c_norm*b_norm*d_norm) +dot(c,d)/(a_norm,b_norm,c_norm,d_norm) ));
        tripletList.push_back(T(  V_index , vj.getIndex()  , -1.0*dot(a,c)/(a_norm*c_norm*b_norm*d_norm) -dot(c,d)/(a_norm,b_norm,c_norm,d_norm) ));
        
        tripletList.push_back(T(  V_index , vi.getIndex()  , -1.0*dot(a,c)/(a_norm*c_norm*b_norm*d_norm) -dot(c,d)/(a_norm,b_norm,c_norm,d_norm) ));
        tripletList.push_back(T(  V_index , vl.getIndex()  , dot(a,c)/(a_norm*c_norm*b_norm*d_norm) +dot(c,d)/(a_norm,b_norm,c_norm,d_norm) ));
        
        tripletList.push_back(T(  V_index , vk.getIndex()  , -1.0*cos_beta/(d_norm*d_norm) +dot(b,c)/(a_norm,b_norm,c_norm,d_norm)  ));
        tripletList.push_back(T(  V_index , vi.getIndex()  , cos_beta/(d_norm*d_norm) -dot(b,c)/(a_norm,b_norm,c_norm,d_norm)  ));
        

        }
  



    }
    K.setFromTriplets(tripletList.begin(),tripletList.end());




    return K;    
}




/*
 * Builds the sparse diagonal mass matrix containing the barycentric dual area of each vertex.
 *
 * Input:
 * Returns: Sparse mass matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {
    SparseMatrix<double> Areas(mesh.nVertices(),mesh.nVertices());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    // TODO
    for(size_t index = 0; index<mesh.nVertices();index++){
    tripletList.push_back(T(index,index,barycentricDualArea(mesh.vertex(index))));    
    }
    Areas.setFromTriplets(tripletList.begin(),tripletList.end());
    
    
    return Areas; // placeholder
}



SparseMatrix<double> VertexPositionGeometry::massMatrix3d() const {
    SparseMatrix<double> M=massMatrix();
    size_t N_vert=mesh.nVertices();
    SparseMatrix<double> M3(3*N_vert,3*N_vert);
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    // TODO
    for(size_t index = 0; index<mesh.nVertices();index++){
    tripletList.push_back(T(index,index, M.coeffRef(index,index) ));    
    tripletList.push_back(T(index+N_vert,index+N_vert,M.coeffRef(index,index)));    
    tripletList.push_back(T(index+2*N_vert,index+2*N_vert,M.coeffRef(index,index) ));    
    
    }
    M3.setFromTriplets(tripletList.begin(),tripletList.end());
    
    
    return M3; // placeholder
}





/*
 * Builds the sparse complex POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace
 * matrix, multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse complex positive definite Laplace matrix for the mesh.
 */
SparseMatrix<std::complex<double>> VertexPositionGeometry::complexLaplaceMatrix() const {

    // TODO
    return identityMatrix<std::complex<double>>(1); // placeholder
}







/*
 * Compute the center of mass of a mesh.
 */
Vector3 VertexPositionGeometry::centerOfMass() const {

    // // Compute center of mass.
    // Vector3 center = {0.0, 0.0, 0.0};
    // for (Vertex v : mesh.vertices()) {
    //     center += inputVertexPositions[v];
    // }
    // center /= mesh.nVertices();


    double dual_area;
    double tot_area=0;
    Vector3 center = {0.0, 0.0, 0.0};
    
    for(Vertex v : mesh.vertices()){
        dual_area=barycentricDualArea(v);
        center += inputVertexPositions[v]*dual_area;
        tot_area+=dual_area;
    }
    center /= tot_area;





    return center;
}

/*
 * Centers a mesh about the origin.
 * Also rescales the mesh to unit radius if <rescale> == true.
 */
void VertexPositionGeometry::normalize(const Vector3& origin, bool rescale) {

    // Compute center of mass.
    Vector3 center = centerOfMass();

    // Translate to origin [of original mesh].
    double radius = 0;
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] -= center;
        radius = std::max(radius, inputVertexPositions[v].norm());
    }

    // Rescale.
    if (rescale) {
        for (Vertex v : mesh.vertices()) {
            inputVertexPositions[v] /= radius;
        }
    }

    // Translate to origin [of original mesh].
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] += origin;
    }
}








} // namespace surface
} // namespace geometrycentral