// Implement member functions for MeanCurvatureFlow class.
#include "Mem-3dg.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
Mem3DG::Mem3DG(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

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
Vector<double> Mem3DG::buildFlowOperator(double h, double V_bar, double nu, double c0,double P0,double KA,double KB) const {

    // Lets get our target area and curvature
    double V=geometry->totalVolume();
    double A_bar=4*PI*pow(3*V/(4*PI*nu),2.0/3.0);
    double H_bar=sqrt(4*PI/A_bar)*c0/2.0; //Coment this with another comment
    
    double A=geometry->totalArea();
    double D_P=P0*(1/V-1/V_bar);
    double lambda=KA*(A-A_bar )/A_bar;
    // lambda=1.0*KA;
    

    // I will start with pressure.

    double E_pressure= E_Pressure(P0,V_bar);
    double E_area=E_Surface(KA,A_bar);

    // std::cout<< "\n The surface tension energy is "<<E_area <<"\n";    

    // std::cout<< "The osmotic pressure energy is "<<E_pressure <<"\n";    



    
    return h*(Bending(H_bar,KB)+OsmoticPressure(D_P)+SurfaceTension(lambda)); // placeholder
}


Vector<double> Mem3DG::OsmoticPressure(double D_P) const {

    // You have the face normals
    
    size_t index;
    Vector3 Normal;
    size_t N_vert=mesh->nVertices();
    Vector<double> Force(3*N_vert);
    
    for(Vertex v : mesh->vertices()) {
     // do science here
        index=v.getIndex();
        Normal={0,0,0};
        for(Face f : v.adjacentFaces()) {
            Normal+=geometry->faceArea(f)*geometry->faceNormal(f);

        }
        Force.coeffRef(v.getIndex())=Normal.x/3.0;
        Force.coeffRef(v.getIndex()+N_vert)=Normal.y/3.0;
        Force.coeffRef(v.getIndex()+2*N_vert)=Normal.z/3.0;

    }

    // std::cout<< "THe osmotic pressure force in magnitude is: "<< D_P*sqrt(Force.transpose()*Force) <<"\n";
    return D_P*Force;
}


Vector<double> Mem3DG::SurfaceTension(double lambda) const {

    size_t index;
    Vector3 Normal;
    size_t N_vert=mesh->nVertices();
    Vector<double> Force(3*N_vert);
    

        
    for(Vertex v : mesh->vertices()) {
        // do science here
            index=v.getIndex();
            Normal=geometry->vertexNormalMeanCurvature(v);
            Force.coeffRef(index)=Normal.x;
            Force.coeffRef(index+N_vert)=Normal.y;
            Force.coeffRef(index+2*N_vert)=Normal.z;

        }

    // std::cout<< "THe surface tension force in magnitude is: "<< -1*lambda*sqrt(Force.transpose()*Force) <<"\n";
    return -1*lambda*Force;
}

geometrycentral::Vector3 Mem3DG::dihedralAngleGradient(geometrycentral::surface::Halfedge he, geometrycentral::surface::Vertex v){
    double l = geometry->edgeLengths[he.edge()];
    if (he.edge().isBoundary()) {
    return geometrycentral::Vector3{0, 0, 0};
  } else if (he.vertex() == v) {
    return (geometry->halfedgeCotanWeights[he.next().next()] *
                geometry->faceNormals[he.face()] +
            geometry->halfedgeCotanWeights[he.twin().next()] *
                geometry->faceNormals[he.twin().face()]) /
           l;
  } else if (he.next().vertex() == v) {
    return (geometry->halfedgeCotanWeights[he.twin().next().next()] *
                geometry->faceNormals[he.twin().face()] +
            geometry->halfedgeCotanWeights[he.next()] *
                geometry->faceNormals[he.face()]) /
           l;
  } else if (he.next().next().vertex() == v) {
    return (-(geometry->halfedgeCotanWeights[he.next().next()] +
              geometry->halfedgeCotanWeights[he.next()]) *
            geometry->faceNormals[he.face()]) /
           l;
  } else {
    // mem3dg_runtime_error("Unexpected combination of halfedge and vertex!");
    return geometrycentral::Vector3{0, 0, 0};
  }



}


Vector<double> Mem3DG::Bending(double H0,double KB) const {

    size_t index;
    size_t neigh_index;
    size_t N_vert=mesh->nVertices();
    Vector3 e_ji;
    Vector3 Hij;
    Vector3 Kij;
    Vector3 Sij_1;
    Vector3 Sij_2;
    Vector3 F1={0,0,0};
    Vector3 F2={0,0,0};
    Vector3 F3={0,0,0};
    Vector3 F4={0,0,0};
    Vector3 Position_1;
    Vector3 Position_2;

    Vector<double> Force1(3*N_vert);

    Vector<double> Scalar_MC(N_vert);
    // Vector<double> Vertex_area(N_vert);
    double factor;

    for(Vertex v : mesh->vertices()) {
        index=v.getIndex();
        Scalar_MC.coeffRef(index)=geometry->scalarMeanCurvature(v)/geometry->barycentricDualArea(v);
        // Vertex_area.coeffRef(index)=geometry->barycentricDualArea(v);
        }   
    
    // std::cout<< Scalar_MC<< "\n";
    // std::cout<<"The spontaneous is"<< H0 <<"\n";
    for(Vertex v: mesh->vertices()){
        F1={0,0,0};
        F2={0,0,0};
        F3={0,0,0};
        F4={0,0,0};

        index=v.getIndex();
        Position_1=geometry->inputVertexPositions[v];
        for(Halfedge he: v.incomingHalfedges()){

            
            neigh_index=he.tailVertex().getIndex();
            Position_2=geometry->inputVertexPositions[neigh_index];
            e_ji=Position_1-Position_2;  
            
            Kij=(0.5*geometry->dihedralAngle(he)/geometry->edgeLength(he.edge()))*e_ji;
            factor=-1*(Scalar_MC[index]-H0)-1*(Scalar_MC[neigh_index]-H0);
            F1=F1+factor*Kij;


            Hij=0.5*(geometry->cotan(he.twin())+geometry->cotan(he))*e_ji;
            // Lets try the code one 
            // Hij=

            factor=(1/3.0)*(Scalar_MC[index] -H0)*(Scalar_MC[index]+H0)+(2.0/3.0)*(Scalar_MC[neigh_index]-H0)*(Scalar_MC[neigh_index]+H0);

            F2=F2+factor*Hij;

            Sij_1=  0.5*(geometry->cotan(he.twin().next().next())* geometry->faceNormal(he.twin().face()) + geometry->cotan(he.next()) * geometry->faceNormal(he.face()));
            
            // THis says the implementation, but uses outgoing edges
            // Sij_1=  0.5*(geometry->cotan(he.next().next())* geometry->faceNormal(he.face()) + geometry->cotan(he.twin()) * geometry->faceNormal(he.twin().face()));

            factor=-1*(Scalar_MC[index]-H0);
            F3=F3+factor*Sij_1;

            Sij_2=-0.5*(   geometry->cotan(he.twin())*geometry->faceNormal(he.twin().face()) + geometry->cotan(he)*geometry->faceNormal(he.face()));
            factor= -1*(Scalar_MC[neigh_index]-H0);
            F4=F4+factor*Sij_2;



    }
    

    // Force1.coeffRef(index)=F1.x;
    // Force1.coeffRef(index+N_vert)=F1.y;
    // Force1.coeffRef(index+2*N_vert)=F1.z;
    
    Force1.coeffRef(index)=F1.x+F2.x+F3.x+F4.x;
    Force1.coeffRef(index+N_vert)=F1.y+F2.y+F3.y+F4.y;
    Force1.coeffRef(index+2*N_vert)=F1.z+F2.z+F3.z+F4.z;
    

    }

    // std::cout<< "THe bending force in magnitude is: "<< KB*sqrt(Force1.transpose()*Force1) <<"\n";

    return KB*Force1;
}


// 

double Mem3DG::E_Pressure(double P0, double V_bar) const {


    double V = geometry->totalVolume();
    
    return 0.5*P0*(V-V_bar)*(V-V_bar)/(V_bar*V_bar);
}

double Mem3DG::E_Surface(double KA, double A_bar) const {

    double A= geometry->totalArea();

    return 0.5*KA*A*A;

    return 0.5*KA*(A-A_bar)*(A-A_bar)/A_bar;
}


double Mem3DG::E_Bending(double H0,double KB) const{

    double Eb=0;

    // for (Vertex v : mesh->vertices()) {
        

    // }

    return 0;
}





/*
 * Performs mean curvature flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void Mem3DG::integrate(double h, double V_bar, double nu, double c0,double P0,double KA,double KB ) {

    Vector<double> Total_force=buildFlowOperator(h,V_bar,nu,c0,P0,KA,KB);
    Vector3 Update;
    size_t vindex;
    size_t Nvert=mesh->nVertices();
    for (Vertex v : mesh->vertices()) {
        vindex=v.getIndex();

        Update= { Total_force[vindex],Total_force[vindex+Nvert],Total_force[vindex+2*Nvert]  };
        geometry->inputVertexPositions[v] =geometry->inputVertexPositions[v]+ Update ; // placeholder
    }
    
  

}