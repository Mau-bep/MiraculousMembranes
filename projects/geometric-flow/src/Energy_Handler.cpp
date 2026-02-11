

#include "Energy_Handler.h"
#include "Beads.h"
#include <fstream>
#include <omp.h>

#include <chrono>
#include <Eigen/Core>


using namespace std;
using namespace geometrycentral;
using namespace geometrycentral::surface;
 

// THings to note, i should make this structure in main and then give a pointer of it to Mem3DG
// Because if i dont then it will be carried all the time by mem3dg and i dont think thats very smart.


E_Handler::E_Handler(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo){
    mesh = inputMesh;
    geometry = inputGeo;

}

E_Handler::E_Handler(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, std::vector<std::string> inputEnergyNames, std::vector<std::vector<double>> inputEnergyConstants ){
    mesh = inputMesh;
    geometry = inputGeo;
    Energies = inputEnergyNames;
    Energy_constants = inputEnergyConstants;
}

void E_Handler::Add_Bead(Bead *bead){

    Beads.push_back(bead);
    // std::cout<<"THe original direction is " << bead->Bead_I <<" \n";
    // // Beads[Beads.size()-1]->Bead_I = bead->Bead_I;
    // std::cout<<"Moving pointers\n";
    // std::cout<<"Te direction of the I is" << Beads[Beads.size()-1]->Bead_I <<" \n";
    // std::cout<<"The energy constant is "<< Beads[Beads.size()-1]->Bead_I->Energy_constants[0] << "MEDG \n";

}

void E_Handler::Add_Energy(std::string Energy_name, std::vector<double> Constants){
    Energies.push_back(Energy_name);
    Energy_constants.push_back(Constants);
}



void E_Handler::update_face_reference(){
    // std::cout<<"Updating face reference, guess we remeshed (: \n";

    Face_reference = FaceData<double>(*mesh);
    // Damn i need the target area here
    double A_bar = 0.0;
    for(size_t i = 0; i < Energies.size(); i++){
        if(Energies[i] == "Area_constraint"){
            A_bar = Energy_constants[i][1];
            // std::cout<<"The target area is " << A_bar <<  "\n";
            break;
        }
    }
    
    // Now i need to calculate the difference between the current area and the target area
    double Total_A;
    double A_face;
    for(Face f : mesh->faces()){
        A_face = geometry->faceArea(f);
        Total_A += A_face;
        Face_reference[f] = A_face; //We need to remove this 2 later
    }
    if(A_bar>1e-5){
        double factor = A_bar/Total_A; // If we assume this difference is small for all the vertices then we can just add A_diff/mesh->nFaces() to each face reference
        // for(Face f : mesh->faces()){
        // std::cout<<"We are multiplying a factor of " << factor << " to each face reference to account for the area constraint\n";
        Face_reference *= factor;
        // }
    }


    return;
}


double E_Handler::E_Volume_constraint(std::vector<double> Constants) const {
double V = geometry->totalVolume();
double KV = Constants[0];
double V_bar = Constants[1];
return 0.5*KV*(V-V_bar)*(V-V_bar)/(V_bar*V_bar);

}

double E_Handler::E_Area_constraint(std::vector<double> Constants) const {
    double KA = Constants[0];
    double A_bar = Constants[1];
    double A = geometry->totalArea();
    // return 0.5*KA*A*A;
    return 0.5*KA*(A-A_bar)*(A-A_bar)/(A_bar*A_bar);
}

double E_Handler::E_SurfaceTension(std::vector<double> Constants) const{
    return geometry->totalArea()*Constants[0];
}

double E_Handler::E_Bending(std::vector<double> Constants) const {
    double KB = Constants[0];
    double H0 = Constants[1];
    size_t index;
    double Eb=0;
    double H;
    double r_eff2;
    Vector3 Pos;
    double A;
    for(Vertex v : mesh->vertices()) {
        //boundary_fix
        // std::cout<<"boundary \n";
        if(v.isBoundary()) continue;
        // std::cout<<" indxe\n";

        index=v.getIndex();
        // std::cout<<" pos\n";
        Pos = geometry->inputVertexPositions[v];
        // std::cout<<"reff \n";
        r_eff2 = Pos.z*Pos.z + Pos.y*Pos.y ;      
        // std::cout<<"boundary? \n";
        if(r_eff2 > 1.6 && boundary ) continue;
  

        A = geometry->barycentricDualArea(v);
        H=(geometry->scalarMeanCurvature(v)/A-H0);
        // std::cout<<" isnan\n";
        if(std::isnan(H)) continue;
           
  
        Eb+=KB*H*H*A;
        
        }   
    // std::cout<<" done with ?\n";
    // std::cout<<Eb<<"\n";
    return Eb;

}

double E_Handler::E_Bending_2(std::vector<double> Constants) const {
    double KB = Constants[0];
    double H0 = Constants[1];
    size_t index;
    double Eb=0;
    double H;
    double r_eff2;
    Vector3 Pos;
    for(Vertex v : mesh->vertices()) {
        //boundary_fix
        // std::cout<<"boundary \n";
        if(v.isBoundary()) continue;
        // std::cout<<" indxe\n";

        index=v.getIndex();
        // std::cout<<" pos\n";
        Pos = geometry->inputVertexPositions[v];
        // std::cout<<"reff \n";
        r_eff2 = Pos.z*Pos.z + Pos.y*Pos.y ;      
        // std::cout<<"boundary? \n";
        if(r_eff2 > 1.6 && boundary ) continue;
        // std::cout<<" MC and dual area\n";
        H=(geometry->scalarMeanCurvature(v)-H0);
        // std::cout<<" isnan\n";
        if(std::isnan(H)){
          
        //   std::cout<<"Dual area: "<< geometry->barycentricDualArea(v);
        //   std::cout<<"Scalar mean Curv"<< geometry->scalarMeanCurvature(v);
        //   std::cout<<"One of the H is not a number\n";
        continue;
        }        
        // std::cout<<" adding\n";
        Eb+=KB*H*H/geometry->barycentricDualArea(v);
        
        }   
    // std::cout<<" done with ?\n";
    // std::cout<<Eb<<"\n";
    return Eb;

}



double E_Handler::E_Laplace(std::vector<double> Constants) const{

    double E_L = 0.0;

    double Lambda = Constants[0];
    double H0 = Constants[1];
    size_t index;
    Vector3 Pos;
    Vector3 Pos2;
    Eigen::Vector<double,12> Positions;

    double dual_Area;
    double Q1;
    double Q2;
    double Q3;
    double Wij;
    double H;
    Vector3 dx;
    Vector3 Hn;
    for(Vertex v : mesh->vertices()) {
        //boundary_fix
        if(v.isBoundary()) continue;
        
        index=v.getIndex();
        Pos = geometry->inputVertexPositions[v];
        dual_Area = geometry->barycentricDualArea(v);
        
        // I need to iterate over all the edges
        H = geometry->vertexNormalMeanCurvature(v).norm();
        // std::cout<<"One option is " << geometry->vertexNormalMeanCurvature(v).norm()/dual_Area << "\n";
        // std::cout<<"THe other  is " << geometry->scalarMeanCurvature(v)/dual_Area << "\n";
        // std::cout<<"THe curvature could be "<< H/dual_Area << " or " << geometry->scalarMeanCurvature(v)/dual_Area <<"\n";
        E_L += (H/dual_Area-H0  )*(H/dual_Area -H0) *dual_Area;


    }


    return Lambda * E_L;

}

double E_Handler::E_Edge_reg(std::vector<double> Constants) const{

    double E_edge = 0.0;
    double KE = Constants[0];

    Halfedge he; 
    Eigen::Vector<double,9> Positions;
    Eigen::Vector<double,3> Edge_lengths_prev;
    EdgeData<double> Edge_lengths(*mesh);
    for(Edge e: mesh->edges()){
        Edge_lengths[e] = geometry->edgeLength(e);
    }

    for(Face f: mesh->faces()){
        he = f.halfedge();
        if(he.edge().isBoundary() || he.next().edge().isBoundary() || he.next().next().edge().isBoundary()) continue;
        Positions << geometry->inputVertexPositions[he.vertex()].x, geometry->inputVertexPositions[he.vertex()].y, geometry->inputVertexPositions[he.vertex()].z,
                    geometry->inputVertexPositions[he.next().vertex()].x, geometry->inputVertexPositions[he.next().vertex()].y, geometry->inputVertexPositions[he.next().vertex()].z,
                    geometry->inputVertexPositions[he.next().next().vertex()].x, geometry->inputVertexPositions[he.next().next().vertex()].y, geometry->inputVertexPositions[he.next().next().vertex()].z;
        Edge_lengths_prev << Edge_lengths[he.edge()], 
                            Edge_lengths[he.next().next().edge()], 
                            Edge_lengths[he.next().edge()];
        E_edge += KE*geometry->Ej_edge_regular(Positions, Edge_lengths_prev);
    }

    return E_edge;

}
double E_Handler::E_Edge_reg_2(std::vector<double> Constants) const{

    double E_edge = 0.0;
    double KE = Constants[0];

    // Ok perfect so now we add the energy 
    Halfedge he;
    Eigen::Vector<double,12> Positions;
    Eigen::Vector<double,5> Edge_lengths_prev;

    EdgeData<double> Edge_lengths(*mesh);
    for(Edge e: mesh->edges()){
        Edge_lengths[e] = geometry->edgeLength(e);
    }
    
    for(Edge e: mesh->edges()) {
        if(e.isBoundary()) continue;
        // I need to get the edge length
        he = e.halfedge();
        
        Positions << geometry->inputVertexPositions[he.vertex()].x, geometry->inputVertexPositions[he.vertex()].y, geometry->inputVertexPositions[he.vertex()].z,
                    geometry->inputVertexPositions[he.next().vertex()].x, geometry->inputVertexPositions[he.next().vertex()].y, geometry->inputVertexPositions[he.next().vertex()].z,
                    geometry->inputVertexPositions[he.next().next().vertex()].x, geometry->inputVertexPositions[he.next().next().vertex()].y, geometry->inputVertexPositions[he.next().next().vertex()].z,
                    geometry->inputVertexPositions[he.twin().next().next().vertex()].x, geometry->inputVertexPositions[he.twin().next().next().vertex()].y, geometry->inputVertexPositions[he.twin().next().next().vertex()].z;
        Edge_lengths_prev << Edge_lengths[he.edge()], 
                            Edge_lengths[he.next().next().edge()], 
                            Edge_lengths[he.twin().next().edge()],
                            Edge_lengths[he.next().edge()],
                            Edge_lengths[he.twin().next().next().edge()];
        
        E_edge += KE*geometry->Ej_edge_regular(Positions, Edge_lengths_prev);
    }


    return E_edge;
}


double E_Handler::E_Face_reg(std::vector<double> Constants ) const{

    double E_face = 0.0;
    double KF = Constants[0];

    double A;
    for(Face f: mesh->faces()){
        A = geometry->faceArea(f);
        E_face += 0.5*KF*(A-Face_reference[f])*(A-Face_reference[f])/(Face_reference[f]*Face_reference[f]);
    }

    return E_face;

}





VertexData<Vector3> E_Handler::F_Volume_constraint(std::vector<double> Constants) const{
    double KV = Constants[0];
    double V_bar = Constants[1];
    double V = geometry->totalVolume();
    
    size_t index;
    Vector3 Normal;
    size_t N_vert=mesh->nVertices();

    VertexData<Vector3> Force(*mesh);
    for(Vertex v : mesh->vertices()) {
     // do science here
        index=v.getIndex();
        Normal={0,0,0};
        for(Face f : v.adjacentFaces()) {
            Normal+=geometry->faceArea(f)*geometry->faceNormal(f);

        }
        // Force[v.getIndex()]=D_P*Normal/3.0;
        Force[v.getIndex()]=Normal/3.0;
    }




    return -1*KV*((V-V_bar)/(V_bar*V_bar))*Force;
}

VertexData<Vector3> E_Handler::F_Volume_constraint_2(std::vector<double> Constants) const{
    
    double KV = Constants[0];
    double V_bar = Constants[1];
    double V = geometry->totalVolume();

    VertexData<Vector3> Force(*mesh);

    Eigen::Vector<double,9> Positions;
    Eigen::Vector<double,9> Grad;
    std::array<Vertex,3> Vertices;

    Halfedge he;
    Vector3 Force_vector;

    for(Face f: mesh->faces()){
        he = f.halfedge();
        Vertices[0] = he.vertex();
        Vertices[1] = he.next().vertex();
        Vertices[2] = he.next().next().vertex();

        Positions << geometry->inputVertexPositions[Vertices[0]].x , geometry->inputVertexPositions[Vertices[0]].y , geometry->inputVertexPositions[Vertices[0]].z,
                    geometry->inputVertexPositions[Vertices[1]].x, geometry->inputVertexPositions[Vertices[1]].y, geometry->inputVertexPositions[Vertices[1]].z,
                    geometry->inputVertexPositions[Vertices[2]].x, geometry->inputVertexPositions[Vertices[2]].y, geometry->inputVertexPositions[Vertices[2]].z;

        Grad = geometry->gradient_volume(Positions);
        for( size_t i = 0; i < 3; i++){
            Force_vector = Vector3{Grad[3*i], Grad[3*i+1], Grad[3*i+2]};
            Force[Vertices[i]] += Force_vector;
        }
    }

    return -1*KV*((V-V_bar)/(V_bar*V_bar))*Force;


}


VertexData<Vector3> E_Handler::F_Volume(std::vector<double> Constants) const{
    
    double KV = Constants[0];
    // double V_bar = Constants[1];
    // double V = geometry->totalVolume();

    VertexData<Vector3> Force(*mesh);

    Eigen::Vector<double,9> Positions;
    Eigen::Vector<double,9> Grad;
    std::array<Vertex,3> Vertices;

    Halfedge he;
    Vector3 Force_vector;

    for(Face f: mesh->faces()){
        he = f.halfedge();
        Vertices[0] = he.vertex();
        Vertices[1] = he.next().vertex();
        Vertices[2] = he.next().next().vertex();

        Positions << geometry->inputVertexPositions[Vertices[0]].x , geometry->inputVertexPositions[Vertices[0]].y , geometry->inputVertexPositions[Vertices[0]].z,
                    geometry->inputVertexPositions[Vertices[1]].x, geometry->inputVertexPositions[Vertices[1]].y, geometry->inputVertexPositions[Vertices[1]].z,
                    geometry->inputVertexPositions[Vertices[2]].x, geometry->inputVertexPositions[Vertices[2]].y, geometry->inputVertexPositions[Vertices[2]].z;

        Grad = geometry->gradient_volume(Positions);
        // std::cout<<"Grad is "<< Grad.transpose() << "\n";
        for( size_t i = 0; i < 3; i++){
            Force_vector = Vector3{Grad[3*i], Grad[3*i+1], Grad[3*i+2]};
            Force[Vertices[i]] += Force_vector;
            // std::cout<<"The force vector is" << Force_vector << "\n";
        }
    }

    return -1*KV*Force;


}

VertexData<Vector3> E_Handler::F_SurfaceTension(std::vector<double> Constants) const{

    double sigma = Constants[0]; 

    size_t index;
    Vector3 Normal;
    size_t N_vert=mesh->nVertices();
    VertexData<Vector3> Force(*mesh);
    Vector3 u;
    Halfedge he_grad;
        
    for(Vertex v : mesh->vertices()) {
            Normal={0,0,0};
            Force[v]={0,0,0};
            Force[v]=-1*2*geometry->vertexNormalMeanCurvature(v);
            
        }

    // std::cout<< "THe surface tension force in magnitude is: "<< -1*lambda*sqrt(Force.transpose()*Force) <<"\n";
    return sigma*Force;

}

VertexData<Vector3> E_Handler::F_SurfaceTension_2(std::vector<double> Constants) const{

    double sigma = Constants[0];

    VertexData<Vector3> Force(*mesh);
    
    Eigen::Vector<double,9> Positions;
    Eigen::Vector<double,9> Grad;
    std::array<Vertex,3> Vertices;
    
    Halfedge he;
    Vector3 Force_vector;
    for(Face f:  mesh->faces()){
        he = f.halfedge();
        Vertices[0] = he.vertex();
        Vertices[1] = he.next().vertex();
        Vertices[2] = he.next().next().vertex();

        Positions << geometry->inputVertexPositions[Vertices[0]].x , geometry->inputVertexPositions[Vertices[0]].y , geometry->inputVertexPositions[Vertices[0]].z,
                    geometry->inputVertexPositions[Vertices[1]].x, geometry->inputVertexPositions[Vertices[1]].y, geometry->inputVertexPositions[Vertices[1]].z,
                    geometry->inputVertexPositions[Vertices[2]].x, geometry->inputVertexPositions[Vertices[2]].y, geometry->inputVertexPositions[Vertices[2]].z;
        Grad = geometry->gradient_triangle_area(Positions);
        for( size_t i = 0; i < 3; i++){
            Force_vector = Vector3{Grad[3*i], Grad[3*i+1], Grad[3*i+2]};
            if( Vertices[i].isBoundary() ) Force[Vertices[i]] = {0,0,0};
            else Force[Vertices[i]] += -1*sigma*Force_vector;
        }

    }

    return Force;


}


VertexData<Vector3> E_Handler::F_Area_constraint(std::vector<double> Constants) const{
    double KA = Constants[0];
    double A_bar = Constants[1];
    double A = geometry->totalArea();
    return ((A-A_bar)/(A_bar*A_bar))*F_SurfaceTension_2({KA});


}
//This function correctly implements the spontaneous curvature
VertexData<Vector3> E_Handler::F_Bending(std::vector<double> Constants) const{


    double KB = Constants[0];
    double H0 = Constants[1];


    VertexData<Vector3> Force(*mesh);

    Eigen::Vector<double,6> Positions_edge;
    Eigen::Vector<double,9> Positions_face;
    Eigen::Vector<double,12> Positions_dihedral;

    Eigen::Vector<double,6> Grad_edge;
    Eigen::Vector<double,9> Grad_face;
    Eigen::Vector<double,12> Grad_dihedral;

    std::array<Vertex,2> Vertices_edge;
    std::array<Vertex,3> Vertices_face;
    std::array<Vertex,4> Vertices_dihedral;

    Halfedge he;
    Vector3 Force_vector;
    
    double constant;

    // So there are two terms, one that is the sum on the faces and one that is the sum on the edges, lets do the faces first
    // I think its best if i store the dual areas and the scalar mean curvatures in vectors.
    VertexData<double> Dual_areas(*mesh,0.0);
    VertexData<double> Scalar_MC(*mesh,0.0);
    for(Vertex v : mesh->vertices()){
        Dual_areas[v] = geometry->barycentricDualArea(v); //Ai
        Scalar_MC[v] = geometry->scalarMeanCurvature(v); //Hi
    }


    for(Face f: mesh->faces()){
        he = f.halfedge();

        Vertices_face[0] = he.vertex();
        Vertices_face[1] = he.next().vertex();
        Vertices_face[2] = he.next().next().vertex();

        if(Vertices_face[0].isBoundary() || Vertices_face[1].isBoundary() || Vertices_face[2].isBoundary() ) continue;

        Positions_face << geometry->inputVertexPositions[Vertices_face[0]].x , geometry->inputVertexPositions[Vertices_face[0]].y , geometry->inputVertexPositions[Vertices_face[0]].z,
                        geometry->inputVertexPositions[Vertices_face[1]].x, geometry->inputVertexPositions[Vertices_face[1]].y, geometry->inputVertexPositions[Vertices_face[1]].z,
                        geometry->inputVertexPositions[Vertices_face[2]].x, geometry->inputVertexPositions[Vertices_face[2]].y, geometry->inputVertexPositions[Vertices_face[2]].z;
        
        Grad_face = geometry->gradient_triangle_area(Positions_face);

        constant =  - 1*(Scalar_MC[Vertices_face[0]]*Scalar_MC[Vertices_face[0]]/(Dual_areas[Vertices_face[0]]*Dual_areas[Vertices_face[0]] ) -H0*H0)/3.0 
                    - 1*(Scalar_MC[Vertices_face[1]]*Scalar_MC[Vertices_face[1]]/(Dual_areas[Vertices_face[1]]*Dual_areas[Vertices_face[1]] ) -H0*H0)/3.0 
                    - 1*(Scalar_MC[Vertices_face[2]]*Scalar_MC[Vertices_face[2]]/(Dual_areas[Vertices_face[2]]*Dual_areas[Vertices_face[2]] ) -H0*H0)/3.0;

        Grad_face *= constant;
        for(size_t i = 0; i < 3; i++){
            Force_vector = Vector3{Grad_face[3*i], Grad_face[3*i+1], Grad_face[3*i+2]};
            if( Vertices_face[i].isBoundary() ) Force[Vertices_face[i]] = {0,0,0};
            else Force[Vertices_face[i]] += -1*KB*Force_vector;
        }

    }

    // Now we need the quantities for the second term which is a summation on the edges
    double dih;
    double lij;

    Vector3 Vectorsum1 = {0, 0, 0};
    Vector3 Vectorsum2 = {0, 0, 0};

    for(Edge e: mesh->edges()){
        
        if(e.isBoundary()) continue;

        Vertices_edge[0] = e.halfedge().vertex();
        Vertices_edge[1] = e.halfedge().twin().vertex();
        

        Positions_edge << geometry->inputVertexPositions[Vertices_edge[0]].x , geometry->inputVertexPositions[Vertices_edge[0]].y , geometry->inputVertexPositions[Vertices_edge[0]].z,
                        geometry->inputVertexPositions[Vertices_edge[1]].x, geometry->inputVertexPositions[Vertices_edge[1]].y, geometry->inputVertexPositions[Vertices_edge[1]].z;
        
        Vertices_dihedral[0] = e.halfedge().vertex();
        Vertices_dihedral[1] = e.halfedge().next().vertex();
        Vertices_dihedral[3] = e.halfedge().next().next().vertex();
        Vertices_dihedral[2] = e.halfedge().twin().next().next().vertex();
        
        Positions_dihedral << geometry->inputVertexPositions[Vertices_dihedral[0]].x , geometry->inputVertexPositions[Vertices_dihedral[0]].y , geometry->inputVertexPositions[Vertices_dihedral[0]].z,
                            geometry->inputVertexPositions[Vertices_dihedral[1]].x, geometry->inputVertexPositions[Vertices_dihedral[1]].y, geometry->inputVertexPositions[Vertices_dihedral[1]].z,
                            geometry->inputVertexPositions[Vertices_dihedral[2]].x, geometry->inputVertexPositions[Vertices_dihedral[2]].y, geometry->inputVertexPositions[Vertices_dihedral[2]].z,
                            geometry->inputVertexPositions[Vertices_dihedral[3]].x, geometry->inputVertexPositions[Vertices_dihedral[3]].y, geometry->inputVertexPositions[Vertices_dihedral[3]].z;
        
        
        
        Grad_edge = geometry->gradient_edge_length(Positions_edge);
        Grad_dihedral = geometry->gradient_dihedral_angle(Positions_dihedral);


        dih = geometry->dihedralAngle(e.halfedge());
        lij = geometry->edgeLength(e);

        // std::cout<<"THis is one dihedral" << dih <<" and the other one is " << geometry->Dihedral_angle(Positions_dihedral) <<" \n";
        
        constant =  (0.5)*(Scalar_MC[Vertices_edge[0]]/(Dual_areas[Vertices_edge[0]])-H0) +
                    (0.5)*(Scalar_MC[Vertices_edge[1]]/(Dual_areas[Vertices_edge[1]])-H0);
                    
        Grad_edge *= constant*dih;
        Grad_dihedral *= constant*lij;

        // Vectorsum1={0, 0, 0};
        for(size_t i = 0; i < 2; i++){
            Force_vector = Vector3{Grad_edge[3*i], Grad_edge[3*i+1], Grad_edge[3*i+2]};
            if( Vertices_edge[i].isBoundary() ) Force[Vertices_edge[i]] = {0,0,0};
            else Force[Vertices_edge[i]] += -1*KB*Force_vector;
        }
        
        for(size_t i = 0; i < 4; i++){
            Force_vector = Vector3{Grad_dihedral[3*i], Grad_dihedral[3*i+1], Grad_dihedral[3*i+2]};

            if( Vertices_dihedral[i].isBoundary() ) Force[Vertices_dihedral[i]] = {0,0,0};
            else Force[Vertices_dihedral[i]] += -1*KB*Force_vector;
        }

    }

    return Force;



}

VertexData<Vector3> E_Handler::F_Bending_2(std::vector<double> Constants) const{

    double KB = Constants[0];
    double H0 = Constants[1];


    VertexData<Vector3> Force(*mesh);

    Eigen::Vector<double,6> Positions_edge;
    Eigen::Vector<double,9> Positions_face;
    Eigen::Vector<double,12> Positions_dihedral;

    Eigen::Vector<double,6> Grad_edge;
    Eigen::Vector<double,9> Grad_face;
    Eigen::Vector<double,12> Grad_dihedral;

    std::array<Vertex,2> Vertices_edge;
    std::array<Vertex,3> Vertices_face;
    std::array<Vertex,4> Vertices_dihedral;

    Halfedge he;
    Vector3 Force_vector;
    
    double constant;

    // So there are two terms, one that is the sum on the faces and one that is the sum on the edges, lets do the faces first
    // I think its best if i store the dual areas and the scalar mean curvatures in vectors.
    VertexData<double> Dual_areas(*mesh,0.0);
    VertexData<double> Scalar_MC(*mesh,0.0);
    for(Vertex v : mesh->vertices()){
        Dual_areas[v] = geometry->barycentricDualArea(v);
        Scalar_MC[v] = geometry->scalarMeanCurvature(v);
    }


    for(Face f: mesh->faces()){
        he = f.halfedge();

        Vertices_face[0] = he.vertex();
        Vertices_face[1] = he.next().vertex();
        Vertices_face[2] = he.next().next().vertex();

        if(Vertices_face[0].isBoundary() || Vertices_face[1].isBoundary() || Vertices_face[2].isBoundary() ) continue;

        Positions_face << geometry->inputVertexPositions[Vertices_face[0]].x , geometry->inputVertexPositions[Vertices_face[0]].y , geometry->inputVertexPositions[Vertices_face[0]].z,
                        geometry->inputVertexPositions[Vertices_face[1]].x, geometry->inputVertexPositions[Vertices_face[1]].y, geometry->inputVertexPositions[Vertices_face[1]].z,
                        geometry->inputVertexPositions[Vertices_face[2]].x, geometry->inputVertexPositions[Vertices_face[2]].y, geometry->inputVertexPositions[Vertices_face[2]].z;
        
        Grad_face = geometry->gradient_triangle_area(Positions_face);

        constant =  - 1*(Scalar_MC[Vertices_face[0]] - H0)*(Scalar_MC[Vertices_face[0]] - H0)/(3.0*Dual_areas[Vertices_face[0]]*Dual_areas[Vertices_face[0]] ) 
                    - 1*(Scalar_MC[Vertices_face[1]] - H0)*(Scalar_MC[Vertices_face[1]] - H0)/(3.0*Dual_areas[Vertices_face[1]]*Dual_areas[Vertices_face[1]] ) 
                    - 1*(Scalar_MC[Vertices_face[2]] - H0)*(Scalar_MC[Vertices_face[2]] - H0)/(3.0*Dual_areas[Vertices_face[2]]*Dual_areas[Vertices_face[2]] );

        Grad_face *= constant;
        for(size_t i = 0; i < 3; i++){
            Force_vector = Vector3{Grad_face[3*i], Grad_face[3*i+1], Grad_face[3*i+2]};
            if( Vertices_face[i].isBoundary() ) Force[Vertices_face[i]] = {0,0,0};
            else Force[Vertices_face[i]] += -1*KB*Force_vector;
        }

    }

    // Now we need the quantities for the second term which is a summation on the edges
    double dih;
    double lij;

    Vector3 Vectorsum1 = {0, 0, 0};
    Vector3 Vectorsum2 = {0, 0, 0};

    for(Edge e: mesh->edges()){
        
        if(e.isBoundary()) continue;

        Vertices_edge[0] = e.halfedge().vertex();
        Vertices_edge[1] = e.halfedge().twin().vertex();
        

        Positions_edge << geometry->inputVertexPositions[Vertices_edge[0]].x , geometry->inputVertexPositions[Vertices_edge[0]].y , geometry->inputVertexPositions[Vertices_edge[0]].z,
                        geometry->inputVertexPositions[Vertices_edge[1]].x, geometry->inputVertexPositions[Vertices_edge[1]].y, geometry->inputVertexPositions[Vertices_edge[1]].z;
        
        Vertices_dihedral[0] = e.halfedge().vertex();
        Vertices_dihedral[1] = e.halfedge().next().vertex();
        Vertices_dihedral[3] = e.halfedge().next().next().vertex();
        Vertices_dihedral[2] = e.halfedge().twin().next().next().vertex();
        
        Positions_dihedral << geometry->inputVertexPositions[Vertices_dihedral[0]].x , geometry->inputVertexPositions[Vertices_dihedral[0]].y , geometry->inputVertexPositions[Vertices_dihedral[0]].z,
                            geometry->inputVertexPositions[Vertices_dihedral[1]].x, geometry->inputVertexPositions[Vertices_dihedral[1]].y, geometry->inputVertexPositions[Vertices_dihedral[1]].z,
                            geometry->inputVertexPositions[Vertices_dihedral[2]].x, geometry->inputVertexPositions[Vertices_dihedral[2]].y, geometry->inputVertexPositions[Vertices_dihedral[2]].z,
                            geometry->inputVertexPositions[Vertices_dihedral[3]].x, geometry->inputVertexPositions[Vertices_dihedral[3]].y, geometry->inputVertexPositions[Vertices_dihedral[3]].z;
        
        
        
        Grad_edge = geometry->gradient_edge_length(Positions_edge);
        Grad_dihedral = geometry->gradient_dihedral_angle(Positions_dihedral);


        dih = geometry->dihedralAngle(e.halfedge());
        lij = geometry->edgeLength(e);

        // std::cout<<"THis is one dihedral" << dih <<" and the other one is " << geometry->Dihedral_angle(Positions_dihedral) <<" \n";
        
        constant =  (0.5)*(Scalar_MC[Vertices_edge[0]]-H0)/(Dual_areas[Vertices_edge[0]]) +
                    (0.5)*(Scalar_MC[Vertices_edge[1]]-H0)/(Dual_areas[Vertices_edge[1]]);
                    
        Grad_edge *= constant*dih;
        Grad_dihedral *= constant*lij;

        // Vectorsum1={0, 0, 0};
        for(size_t i = 0; i < 2; i++){
            Force_vector = Vector3{Grad_edge[3*i], Grad_edge[3*i+1], Grad_edge[3*i+2]};
            if( Vertices_edge[i].isBoundary() ) Force[Vertices_edge[i]] = {0,0,0};
            else Force[Vertices_edge[i]] += -1*KB*Force_vector;
        }
        for(size_t i = 0; i < 4; i++){
            Force_vector = Vector3{Grad_dihedral[3*i], Grad_dihedral[3*i+1], Grad_dihedral[3*i+2]};
            if( Vertices_dihedral[i].isBoundary() ) Force[Vertices_dihedral[i]] = {0,0,0};
            else Force[Vertices_dihedral[i]] += -1*KB*Force_vector;
        }

    }

    return Force;



}


VertexData<Vector3> E_Handler::F_Laplace(std::vector<double> Constants) const{

    double Lambda = Constants[0];
    double H0 = Constants[1];

    VertexData<Vector3> Force(*mesh);

    Eigen::Vector<double,12> Positions_quad;
    Eigen::Vector<double,9> Positions_face;

    Eigen::Vector<double,12> Grad_quad;
    Eigen::Vector<double,9> Grad_face;
    Eigen::Vector<double,6> Grad_edge;

    std::array<Vertex,3> Vertices_face;
    std::array<Vertex,4> Vertices_quad;
   
    EdgeData<double> Wij(*mesh,0.0);
    VertexData<double> Dual_areas(*mesh,0.0);
    VertexData<double> Hi(*mesh,0.0);

    Halfedge he;
    double constant;

    Vector3 Q;
    Vector3 Force_vector;
    Vector3 dx;

    for(Vertex v : mesh->vertices()){
        Dual_areas[v] = geometry->barycentricDualArea(v);
        Hi[v] = geometry->vertexNormalMeanCurvature(v).norm();
    }

    for(Edge e: mesh->edges()){
        he = e.halfedge();
        Positions_quad << geometry->inputVertexPositions[he.vertex()].x, geometry->inputVertexPositions[he.vertex()].y, geometry->inputVertexPositions[he.vertex()].z,
                        geometry->inputVertexPositions[he.next().vertex()].x, geometry->inputVertexPositions[he.next().vertex()].y, geometry->inputVertexPositions[he.next().vertex()].z,
                        geometry->inputVertexPositions[he.next().next().vertex()].x, geometry->inputVertexPositions[he.next().next().vertex()].y, geometry->inputVertexPositions[he.next().next().vertex()].z,
                        geometry->inputVertexPositions[he.twin().next().next().vertex()].x, geometry->inputVertexPositions[he.twin().next().next().vertex()].y, geometry->inputVertexPositions[he.twin().next().next().vertex()].z;
        Wij[e] = geometry->Cotan_weight(Positions_quad);
        
    }
    // Ok so we have the two terms and we need to start calculating stuff

    for(Vertex v: mesh->vertices()){

        Q = Vector3{0,0,0};

        Q = -1*geometry->vertexNormalMeanCurvature(v);
        
        constant = -(1.0/3.0)*(Hi[v]*Hi[v]/(Dual_areas[v]*Dual_areas[v])-H0*H0);
        // 
        for(Face f: v.adjacentFaces()){
            he = f.halfedge();
            Vertices_face[0] = he.vertex();
            Vertices_face[1] = he.next().vertex();
            Vertices_face[2] = he.next().next().vertex();

            Positions_face << geometry->inputVertexPositions[Vertices_face[0]].x , geometry->inputVertexPositions[Vertices_face[0]].y , geometry->inputVertexPositions[Vertices_face[0]].z,
                            geometry->inputVertexPositions[Vertices_face[1]].x, geometry->inputVertexPositions[Vertices_face[1]].y, geometry->inputVertexPositions[Vertices_face[1]].z,
                            geometry->inputVertexPositions[Vertices_face[2]].x, geometry->inputVertexPositions[Vertices_face[2]].y, geometry->inputVertexPositions[Vertices_face[2]].z;
            
            Grad_face = geometry->gradient_triangle_area(Positions_face);

            for(size_t i = 0; i < 3; i++){
                Force_vector = Vector3{Grad_face[i*3], Grad_face[i*3+1], Grad_face[i*3+2]};
                Force_vector *= constant;
                // Force_vector *= Lambda;
                Force[Vertices_face[i]] -= Force_vector;
            }
        }


        // THe second term is also tricky more so even.

        // I  have Q at this point, which is 

        constant = (Hi[v]/Dual_areas[v]-H0)/(2*Hi[v]);

        for(Halfedge he: v.outgoingHalfedges()){
            // 
            Vertices_quad[0] = he.vertex();
            Vertices_quad[1] = he.next().vertex();
            Vertices_quad[2] = he.next().next().vertex();
            Vertices_quad[3] = he.twin().next().next().vertex();


            Positions_quad << geometry->inputVertexPositions[Vertices_quad[0]].x, geometry->inputVertexPositions[Vertices_quad[0]].y, geometry->inputVertexPositions[Vertices_quad[0]].z,
                            geometry->inputVertexPositions[Vertices_quad[1]].x, geometry->inputVertexPositions[Vertices_quad[1]].y, geometry->inputVertexPositions[Vertices_quad[1]].z,
                            geometry->inputVertexPositions[Vertices_quad[2]].x, geometry->inputVertexPositions[Vertices_quad[2]].y, geometry->inputVertexPositions[Vertices_quad[2]].z,
                            geometry->inputVertexPositions[Vertices_quad[3]].x, geometry->inputVertexPositions[Vertices_quad[3]].y, geometry->inputVertexPositions[Vertices_quad[3]].z;

            Grad_quad = geometry->gradient_cotan_weight(Positions_quad);

            dx = geometry->inputVertexPositions[Vertices_quad[1]] - geometry->inputVertexPositions[v];
            for(size_t i = 0; i < 4; i++){
                Force_vector = Vector3{Grad_quad[i*3], Grad_quad[i*3+1], Grad_quad[i*3+2]};
                
                Force_vector *= constant*(dot(Q,dx));

                Force[Vertices_quad[i]] -= Force_vector;
            }

            // What about the other term now the other term is 
            
            Force_vector = -1*Q*Wij[he.edge()];
            Force_vector *= constant;
            Force[he.vertex()] -= Force_vector;
            Force[he.twin().vertex()] += Force_vector;





        }




    }

    // That is the first term
    return Force*Lambda;



}

VertexData<Vector3> E_Handler::F_Edge_reg(std::vector<double> Constants) const{

    VertexData<Vector3> Force(*mesh);
    double KE = Constants[0];

    // Ok perfect so now we add the energy
    Halfedge he;

    Eigen::Vector<double,9> Positions;
    Eigen::Vector<double,3> Edge_lengths_prev;
    Eigen::Vector<double,9> Grad_E;

    std::array<Vertex,3> Vertices;
    Vector3 Force_vector;
    EdgeData<double> Edge_lengths(*mesh);
    for(Edge e: mesh->edges()){
        Edge_lengths[e] = geometry->edgeLength(e);
    }
    // Till here everything is the same
    for(Face f: mesh->faces()) {

        
        // I need to get the edge length
        he = f.halfedge();

        Vertices[0] = he.vertex();
        Vertices[1] = he.next().vertex();
        Vertices[2] = he.next().next().vertex();
        if(Vertices[0].isBoundary() || Vertices[1].isBoundary() || Vertices[2].isBoundary() ) continue;
        
        Positions << geometry->inputVertexPositions[he.vertex()].x, geometry->inputVertexPositions[he.vertex()].y, geometry->inputVertexPositions[he.vertex()].z,
                    geometry->inputVertexPositions[he.next().vertex()].x, geometry->inputVertexPositions[he.next().vertex()].y, geometry->inputVertexPositions[he.next().vertex()].z,
                    geometry->inputVertexPositions[he.next().next().vertex()].x, geometry->inputVertexPositions[he.next().next().vertex()].y, geometry->inputVertexPositions[he.next().next().vertex()].z;
        Edge_lengths_prev << Edge_lengths[he.edge()], 
                            Edge_lengths[he.next().next().edge()], 
                            Edge_lengths[he.next().edge()];
        
        Grad_E = geometry->gradient_edge_regular(Positions, Edge_lengths_prev);

        // Now we have the gradient of the edge regularization energy, we need to add it to the force
        for(size_t i = 0; i < 3; i++){
            Force_vector = Vector3{Grad_E[3*i], Grad_E[3*i+1], Grad_E[3*i+2]};
            Force_vector *= KE;
            if( Vertices[i].isBoundary() ) Force[Vertices[i]] = {0,0,0};
            else Force[Vertices[i]] -= Force_vector;
        }
    
    }

    return Force;   
}

VertexData<Vector3> E_Handler::F_Edge_reg_2(std::vector<double> Constants) const{

    VertexData<Vector3> Force(*mesh);
    double KE = Constants[0];

    // Ok perfect so now we add the energy 
    Halfedge he;

    // std::vector<int> 
    Eigen::Vector<double,12> Positions;
    Eigen::Vector<double,5> Edge_lengths_prev;
    Eigen::Vector<double,12> Grad_E;

    std::array<Vertex,4> Vertices;
    Vector3 Force_vector;


    EdgeData<double> Edge_lengths(*mesh);
    for(Edge e: mesh->edges()){
        Edge_lengths[e] = geometry->edgeLength(e);
    }
    
    // Till here everything is the same


    for(Edge e: mesh->edges()) {
        // I need to get the edge length
        he = e.halfedge();

        Vertices[0] = he.vertex();
        Vertices[1] = he.next().vertex();
        Vertices[2] = he.next().next().vertex();
        Vertices[3] = he.twin().next().next().vertex();
        
        Positions << geometry->inputVertexPositions[he.vertex()].x, geometry->inputVertexPositions[he.vertex()].y, geometry->inputVertexPositions[he.vertex()].z,
                    geometry->inputVertexPositions[he.next().vertex()].x, geometry->inputVertexPositions[he.next().vertex()].y, geometry->inputVertexPositions[he.next().vertex()].z,
                    geometry->inputVertexPositions[he.next().next().vertex()].x, geometry->inputVertexPositions[he.next().next().vertex()].y, geometry->inputVertexPositions[he.next().next().vertex()].z,
                    geometry->inputVertexPositions[he.twin().next().next().vertex()].x, geometry->inputVertexPositions[he.twin().next().next().vertex()].y, geometry->inputVertexPositions[he.twin().next().next().vertex()].z;
        Edge_lengths_prev << Edge_lengths[he.edge()], 
                            Edge_lengths[he.next().next().edge()], 
                            Edge_lengths[he.twin().next().edge()],
                            Edge_lengths[he.next().edge()],
                            Edge_lengths[he.twin().next().next().edge()];
        
        
        Grad_E = geometry->gradient_edge_regular(Positions, Edge_lengths_prev);


        // Now we have the gradient of the edge regularization energy, we need to add it to the force
        for(size_t i = 0; i < 4; i++){
            Force_vector = Vector3{Grad_E[3*i], Grad_E[3*i+1], Grad_E[3*i+2]};
            Force_vector *= KE;
            if( Vertices[i].isBoundary() ) Force[Vertices[i]] = {0,0,0};
            else Force[Vertices[i]] -= Force_vector;
          
        }

    }


    return Force;
}


VertexData<Vector3> E_Handler::F_Face_reg(std::vector<double> Constants) const{

    double KF = Constants[0];

    VertexData<Vector3> Force(*mesh);
    
    Eigen::Vector<double,9> Positions;
    Eigen::Vector<double,9> Grad;
    std::array<Vertex,3> Vertices;
    
    Halfedge he;
    Vector3 Force_vector;
    double A;
    for(Face f:  mesh->faces()){
        A = geometry->faceArea(f);
        he = f.halfedge();
        Vertices[0] = he.vertex();
        Vertices[1] = he.next().vertex();
        Vertices[2] = he.next().next().vertex();

        Positions << geometry->inputVertexPositions[Vertices[0]].x , geometry->inputVertexPositions[Vertices[0]].y , geometry->inputVertexPositions[Vertices[0]].z,
                    geometry->inputVertexPositions[Vertices[1]].x, geometry->inputVertexPositions[Vertices[1]].y, geometry->inputVertexPositions[Vertices[1]].z,
                    geometry->inputVertexPositions[Vertices[2]].x, geometry->inputVertexPositions[Vertices[2]].y, geometry->inputVertexPositions[Vertices[2]].z;
        Grad = geometry->gradient_triangle_area(Positions);
        for( size_t i = 0; i < 3; i++){
            Force_vector = Vector3{Grad[3*i], Grad[3*i+1], Grad[3*i+2]};
            if( Vertices[i].isBoundary() ) {
                std::cout<<"THere is a boundary in this sim\n";
                Force[Vertices[i]] = {0,0,0};}

            else Force[Vertices[i]] += -1*KF*Force_vector*(A- Face_reference[f])/(Face_reference[f]*Face_reference[f]);
        }

    }

    return Force;


}



SparseMatrix<double> E_Handler::H_SurfaceTension(std::vector<double> Constants){

    // Ok so this functino will assemble the Hessi an for the surface tension energy 
    int nVerts = mesh->nVertices();
    int nBeads = Beads.size();
    double KA = Constants[0];
    // Eigen::MatrixXd Hessian = Eigen::MatrixXd::Zero(nVerts,nVerts);
    SparseMatrix<double> Hessian(3*(nVerts+nBeads),3*(nVerts+nBeads) );
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    // Eigen::Matrix<double, nVerts, nVerts> Hessian;
    Vertex v1;
    Vertex v2;
    Vertex v3;
    Halfedge he;
    size_t index1;
    size_t index2;
    size_t index3;
    array<int,3> indices;
    // 
    for(Face f: mesh->faces()){
        // 
        he = f.halfedge();
        v1 = he.vertex();
        v2 = he.next().vertex();
        v3 = he.next().next().vertex();
        indices[0] = v1.getIndex();
        indices[1] = v2.getIndex();
        indices[2] = v3.getIndex();
        
        // So we will add a modification for boundaries then check if it works :p
        if(v1.isBoundary() || v2.isBoundary() || v3.isBoundary() ) continue;

        Eigen::Vector<double,9> Positions;

        Positions << geometry->inputVertexPositions[v1].x,geometry->inputVertexPositions[v1].y,geometry->inputVertexPositions[v1].z,
                    geometry->inputVertexPositions[v2].x,geometry->inputVertexPositions[v2].y,geometry->inputVertexPositions[v2].z,
                    geometry->inputVertexPositions[v3].x,geometry->inputVertexPositions[v3].y,geometry->inputVertexPositions[v3].z;
         
        Eigen::Matrix<double,9,9> Hessian_block = geometry->hessian_triangle_area(Positions);
        // Now i need to load this quantities onto a bigger matrix ...
        for(int i = 0; i < 9; i++){
            for(int j = 0; j <9; j++){

                // So i am at vertex ij  now the i and j correspond to a vertex  i = 0 1 2 (vertex 1 )  3 4 5 (vertex 2 ) 6 7 8 (vertex 3)
                if( Hessian_block(i,j) > 1e-12 || Hessian_block(i,j) < -1e-12) tripletList.push_back(T(indices[i/3]*3+i%3,indices[j/3]*3+j%3,Hessian_block(i,j)) );  
            }
        } 

    }
    Hessian.setFromTriplets(tripletList.begin(),tripletList.end());
    // That gives me the hessian of the surface tension :O.


    return KA*Hessian;
}

SparseMatrix<double> E_Handler::H_Bending(std::vector<double> Constants) {

    double KB = Constants[0];
    double H0 = Constants[1];

    int N_verts = mesh->nVertices();
    int N_beads = Beads.size();
    // std::cout<<"The number of beads is "<< N_beads <<"\n";
    SparseMatrix<double> Hessian(3*(N_verts+N_beads),3*(N_verts+N_beads));
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

    Eigen::Vector<double,6> Positions_edge;
    Eigen::Vector<double,9> Positions_face;
    Eigen::Vector<double,12> Positions_dihedral;

    Eigen::Vector<double,6> Grad_edge;
    Eigen::Vector<double,9> Grad_face;
    Eigen::Vector<double,12> Grad_dihedral;

    std::array<Vertex,2> Vertices_edge;
    std::array<Vertex,3> Vertices_face;
    std::array<Vertex,4> Vertices_dihedral;

    std::array<Vertex,3 > Vertices_face_2;
    std::array<Vertex,2> Vertices_edge_2;
    std::array<Vertex,4> Vertices_dihedral_2;


    Eigen::Matrix<double,6,6> Hessian_block_edge;
    Eigen::Matrix<double,9,9> Hessian_block_face;
    Eigen::Matrix<double,12,12> Hessian_block_dihedral;

    Halfedge he;
    double constant;

    VertexData<double> Dual_areas(*mesh,0.0);
    VertexData<double> Scalar_MC(*mesh,0.0);
    EdgeData<double> Edge_lengths(*mesh,0.0);
    EdgeData<double> Dihedral_angles(*mesh,0.0);


    for(Vertex v : mesh->vertices()){
        Dual_areas[v] = geometry->barycentricDualArea(v);
        Scalar_MC[v] = geometry->scalarMeanCurvature(v);

    }


    for(Edge e: mesh->edges()){
        Edge_lengths[e] = geometry->edgeLength(e);
        Dihedral_angles[e] = geometry->dihedralAngle(e.halfedge());
    }

    std::vector<Eigen::Vector<double,9>> Gradients_areas;
    std::vector<Eigen::Vector<double,6>> Gradients_edges;
    std::vector<Eigen::Vector<double,12>> Gradients_dihedrals;


    int id= 0;
    for(Face f:mesh -> faces()){
    
        he = f.halfedge();

        Vertices_face[0] = he.vertex();
        Vertices_face[1] = he.next().vertex();
        Vertices_face[2] = he.next().next().vertex();
        Positions_face <<   geometry->inputVertexPositions[Vertices_face[0]].x, geometry->inputVertexPositions[Vertices_face[0]].y, geometry->inputVertexPositions[Vertices_face[0]].z,
                            geometry->inputVertexPositions[Vertices_face[1]].x, geometry->inputVertexPositions[Vertices_face[1]].y, geometry->inputVertexPositions[Vertices_face[1]].z,
                            geometry->inputVertexPositions[Vertices_face[2]].x, geometry->inputVertexPositions[Vertices_face[2]].y, geometry->inputVertexPositions[Vertices_face[2]].z;
        
        Grad_face = geometry->gradient_triangle_area(Positions_face);
        Gradients_areas.push_back(Grad_face);
        
    }

    double dih1 = 0.0;
    double dih2 = 0.0;
    for(Edge e: mesh->edges()){
        Vertices_edge[0] = e.halfedge().vertex();
        Vertices_edge[1] = e.halfedge().twin().vertex();

        Vertices_dihedral[0] = e.halfedge().vertex();
        Vertices_dihedral[1] = e.halfedge().next().vertex();
        Vertices_dihedral[3] = e.halfedge().next().next().vertex();
        Vertices_dihedral[2] = e.halfedge().twin().next().next().vertex();

        Positions_edge << geometry->inputVertexPositions[Vertices_edge[0]].x , geometry->inputVertexPositions[Vertices_edge[0]].y , geometry->inputVertexPositions[Vertices_edge[0]].z,
                        geometry->inputVertexPositions[Vertices_edge[1]].x, geometry->inputVertexPositions[Vertices_edge[1]].y, geometry->inputVertexPositions[Vertices_edge[1]].z;
        Positions_dihedral << geometry->inputVertexPositions[Vertices_dihedral[0]].x , geometry->inputVertexPositions[Vertices_dihedral[0]].y , geometry->inputVertexPositions[Vertices_dihedral[0]].z,
                            geometry->inputVertexPositions[Vertices_dihedral[1]].x, geometry->inputVertexPositions[Vertices_dihedral[1]].y, geometry->inputVertexPositions[Vertices_dihedral[1]].z,
                            geometry->inputVertexPositions[Vertices_dihedral[2]].x, geometry->inputVertexPositions[Vertices_dihedral[2]].y, geometry->inputVertexPositions[Vertices_dihedral[2]].z,
                            geometry->inputVertexPositions[Vertices_dihedral[3]].x, geometry->inputVertexPositions[Vertices_dihedral[3]].y, geometry->inputVertexPositions[Vertices_dihedral[3]].z;

        Grad_edge = geometry->gradient_edge_length(Positions_edge);
        Gradients_edges.push_back(Grad_edge);

        Grad_dihedral = geometry->gradient_dihedral_angle(Positions_dihedral);
        Gradients_dihedrals.push_back(Grad_dihedral);
        // std::cout<<"The edge " << e.getIndex() << "has a edge grad = " << Gradients_edges[e.getIndex()].transpose() <<"\n and dihedral = " << Gradients_dihedrals[e.getIndex()].transpose()<< "\n";
    }

    // I have calculated all the quantities that i will be needing, this takes a lot of memory but lets hope its worth it ahah

    // Lets iterate over the faces   
    Eigen::Matrix<double,9,12> M_9_12;
    Eigen::Matrix<double,9,6> M_9_6;
    Eigen::Matrix<double,9,9> M_9_9;
    Eigen::Matrix<double,6,6> M_6_6;
    Eigen::Matrix<double,6,12> M_6_12;
    Eigen::Matrix<double,12,6> M_12_6;
    Eigen::Matrix<double,12,12> M_12_12;
    Eigen::Matrix<double,6,9> M_6_9;
    Eigen::Matrix<double,12,9> M_12_9;
    size_t max_id = 0;
    for(Face f : mesh->faces()){

        he = f.halfedge();
        Grad_face = Gradients_areas[f.getIndex()];
            
        Vertices_face[0] = he.vertex();
        Vertices_face[1] = he.next().vertex();
        Vertices_face[2] = he.next().next().vertex();

        if(Vertices_face[0].isBoundary() || Vertices_face[1].isBoundary() || Vertices_face[2].isBoundary() ) continue;

        for(Vertex v: f.adjacentVertices()){
            
            constant = -1*(1.0/6.0)*(Scalar_MC[v.getIndex()])/(Dual_areas[v.getIndex()]*Dual_areas[v.getIndex()]);
            for(Edge e: v.adjacentEdges()){
                he = e.halfedge();
                Vertices_dihedral[0] = he.vertex();
                Vertices_dihedral[1] = he.next().vertex();
                Vertices_dihedral[3] = he.next().next().vertex();
                Vertices_dihedral[2] = he.twin().next().next().vertex();

                Vertices_edge[0] = he.vertex();
                Vertices_edge[1] = he.next().vertex();

                M_9_12 = constant*Edge_lengths[e.getIndex()]*Grad_face * Gradients_dihedrals[e.getIndex()].transpose();

                for(size_t row = 0; row < 9; row ++){
                    for(size_t col = 0; col < 12; col++){
                        tripletList.push_back(T(Vertices_face[row/3].getIndex()*3+row%3, Vertices_dihedral[col/3].getIndex()*3+col%3, M_9_12(row,col)));
              
                        
                }
                }

                M_9_6 = constant * Dihedral_angles[e]* Grad_face * Gradients_edges[e.getIndex()].transpose();

                for(size_t row = 0; row < 9; row ++){
                    for(size_t col = 0; col < 6; col++){
                        tripletList.push_back(T(Vertices_face[row/3].getIndex()*3+row%3, Vertices_edge[col/3].getIndex()*3+col%3, M_9_6(row,col)));
                    
                    }
                }

                M_6_9 = constant * Dihedral_angles[e]*Gradients_edges[e.getIndex()]*Grad_face.transpose();

                    for(size_t row = 0; row < 6; row++){
                        for(size_t col = 0; col < 9; col++){
                        tripletList.push_back(T(Vertices_edge[row/3].getIndex()*3+row%3,Vertices_face[col/3].getIndex()*3+col%3, M_6_9(row,col)));
                     
                        }
                    }

                M_12_9 = constant * Edge_lengths[e.getIndex()]*Gradients_dihedrals[e.getIndex()] * Grad_face.transpose();

                    for(size_t row = 0; row < 12; row++){
                        for( size_t col = 0; col < 9; col++){
                            tripletList.push_back(T(Vertices_dihedral[row/3].getIndex()*3+row%3,Vertices_face[col/3].getIndex()*3+col%3,M_12_9(row,col)));
                           
                        }
                    }

            
            
            }

            constant = (2.0/9.0) * (Scalar_MC[v.getIndex()]) * (Scalar_MC[v.getIndex()])/(Dual_areas[v.getIndex()] * Dual_areas[v.getIndex()] * Dual_areas[v.getIndex()]);
            
            
            for(Face f2: v.adjacentFaces()){
                he = f2.halfedge();
                
                Vertices_face_2[0] = he.vertex();
                Vertices_face_2[1] = he.next().vertex();
                Vertices_face_2[2] = he.next().next().vertex();

                M_9_9 = constant * Grad_face * Gradients_areas[f2.getIndex()].transpose();
                    for(size_t row = 0; row < 9; row++){
                        for(size_t col = 0; col < 9; col++){

                            tripletList.push_back(T(Vertices_face[row/3].getIndex()*3+row%3,Vertices_face_2[col/3].getIndex()*3+col%3, M_9_9(row,col) ));
                            // if(isnan(M_9_9(row,col))) std::cout<<" nan flag 3 \n";
                        }
                    }

            }

            // And last but not least we have the hessian term

            constant = -1*(1.0/3.0)*(Scalar_MC[v.getIndex()]*Scalar_MC[v.getIndex()] /(Dual_areas[v.getIndex()]*Dual_areas[v.getIndex()]) - H0*H0 );
            
            Positions_face << geometry->inputVertexPositions[Vertices_face[0]].x , geometry->inputVertexPositions[Vertices_face[0]].y , geometry->inputVertexPositions[Vertices_face[0]].z,
                        geometry->inputVertexPositions[Vertices_face[1]].x, geometry->inputVertexPositions[Vertices_face[1]].y, geometry->inputVertexPositions[Vertices_face[1]].z,
                        geometry->inputVertexPositions[Vertices_face[2]].x, geometry->inputVertexPositions[Vertices_face[2]].y, geometry->inputVertexPositions[Vertices_face[2]].z;
        
            M_9_9 = constant * geometry->hessian_triangle_area(Positions_face);

            for(size_t row = 0; row < 9; row++){
                for(size_t col = 0; col < 9; col++){
                    tripletList.push_back(T(Vertices_face[row/3].getIndex()*3+row%3, Vertices_face[col/3].getIndex()*3+col%3, M_9_9(row,col)));
                    // if(isnan(M_9_9(row,col))) std::cout<<" nan flag 4 \n";
                }
            }


        }
    }
        // And those are all the terms at the sum over all the vertices (On paper is term I)

        for(Edge e: mesh->edges()){
            he = e.halfedge();

            Vertices_edge[0] = he.vertex();
            Vertices_edge[1] = he.twin().vertex();
            
            Vertices_dihedral[0] = he.vertex();
            Vertices_dihedral[1] = he.next().vertex();
            Vertices_dihedral[3] = he.next().next().vertex();
            Vertices_dihedral[2] = he.twin().next().next().vertex();

            if(Vertices_edge[0].isBoundary() || Vertices_edge[1].isBoundary() || Vertices_dihedral[2].isBoundary() || Vertices_dihedral[3].isBoundary() ) continue;

            for(Vertex v: e.adjacentVertices()){
                constant = 1.0/(8.0*Dual_areas[v.getIndex()]);
                for( Edge e2 : v.adjacentEdges()){
                    Vertices_edge_2[0] = e2.halfedge().vertex();
                    Vertices_edge_2[1] = e2.halfedge().twin().vertex();

                    M_6_6 = constant*Dihedral_angles[e.getIndex()]*Dihedral_angles[e2.getIndex()]*Gradients_edges[e.getIndex()]*Gradients_edges[e2.getIndex()].transpose();

                    for(size_t row = 0; row < 6; row++){
                        for(size_t col = 0; col < 6; col++){
                            tripletList.push_back(T(Vertices_edge[row/3].getIndex()*3+row%3, Vertices_edge_2[col/3].getIndex()*3+col%3,M_6_6(row,col)));
                            // if(isnan(M_6_6(row,col))) std::cout<<" nan flag 5 \n";
                        }
                    }                    


                    he = e2.halfedge();
                    Vertices_dihedral_2[0] = he.vertex();
                    Vertices_dihedral_2[1] = he.next().vertex();
                    Vertices_dihedral_2[3] = he.next().next().vertex();
                    Vertices_dihedral_2[2] = he.twin().next().next().vertex(); 

                    M_6_12 = constant* Dihedral_angles[e.getIndex()]*Edge_lengths[e2.getIndex()]*Gradients_edges[e.getIndex()]*Gradients_dihedrals[e2.getIndex()].transpose();

                    for(size_t row = 0; row < 6; row++){
                        for(size_t col =0; col < 12; col++){
                            tripletList.push_back(T(Vertices_edge[row/3].getIndex()*3+row%3, Vertices_dihedral_2[col/3].getIndex()*3+col%3,M_6_12(row,col)));
                            // if(isnan(M_6_12(row,col))) std::cout<<" nan flag 6 \n";
                        }
                    }

                     

                    M_12_6 = constant * Edge_lengths[e.getIndex()]*Dihedral_angles[e2.getIndex()]*Gradients_dihedrals[e.getIndex()]*Gradients_edges[e2.getIndex()].transpose();

                    for(size_t row = 0; row < 12; row++){
                        for(size_t col = 0; col < 6; col++){
                            tripletList.push_back(T(Vertices_dihedral[row/3].getIndex()*3+row%3, Vertices_edge_2[col/3].getIndex()*3+col%3, M_12_6(row,col)));
                            // if(isnan(M_12_6(row,col))) std::cout<<" nan flag 7 \n";
                        }
                    }

                    M_12_12 = constant* Edge_lengths[e.getIndex()] * Edge_lengths[e2.getIndex()]* Gradients_dihedrals[e.getIndex()]*Gradients_dihedrals[e2.getIndex()].transpose();

                    for(size_t row = 0; row < 12; row++){
                        for(size_t col = 0; col < 12; col++){
                            tripletList.push_back(T(Vertices_dihedral[row/3].getIndex()*3+row%3,Vertices_dihedral_2[col/3].getIndex()*3+col%3,M_12_12(row,col)));
                            // if(isnan(M_12_12(row,col))) std::cout<<"  8 \n";
                        }
                    }




                }


                // Those are the first terms
                // I MOVED THIS TERMS TO THE OTHER SUMMATION






            constant = (1.0/2.0)*(Scalar_MC[v.getIndex()]/(Dual_areas[v.getIndex()]) -H0) ;

            M_6_12 = constant * Gradients_edges[e.getIndex()]*Gradients_dihedrals[e.getIndex()].transpose();

            for( size_t row = 0; row < 6; row++){
                for(size_t col = 0; col < 12; col++){
                    tripletList.push_back(T(Vertices_edge[row/3].getIndex()*3+row%3,Vertices_dihedral[col/3].getIndex()*3+col%3,M_6_12(row,col)));
            
                }
            }


            M_12_6 = constant * Gradients_dihedrals[e.getIndex()]* Gradients_edges[e.getIndex()].transpose();

            for( size_t row = 0; row < 12; row++){
                for(size_t col = 0; col < 6; col++){
                    tripletList.push_back(T(Vertices_dihedral[row/3].getIndex()*3+row%3,Vertices_edge[col/3].getIndex()*3+col%3,M_12_6(row,col)));
                }
            }

            Positions_dihedral <<   geometry->inputVertexPositions[Vertices_dihedral[0]].x , geometry->inputVertexPositions[Vertices_dihedral[0]].y , geometry->inputVertexPositions[Vertices_dihedral[0]].z,
                                    geometry->inputVertexPositions[Vertices_dihedral[1]].x, geometry->inputVertexPositions[Vertices_dihedral[1]].y, geometry->inputVertexPositions[Vertices_dihedral[1]].z,
                                    geometry->inputVertexPositions[Vertices_dihedral[2]].x, geometry->inputVertexPositions[Vertices_dihedral[2]].y, geometry->inputVertexPositions[Vertices_dihedral[2]].z,
                                    geometry->inputVertexPositions[Vertices_dihedral[3]].x, geometry->inputVertexPositions[Vertices_dihedral[3]].y, geometry->inputVertexPositions[Vertices_dihedral[3]].z;

            


            M_12_12 = constant * Edge_lengths[e.getIndex()]* geometry->hessian_dihedral_angle(Positions_dihedral);

            for( size_t row = 0; row < 12; row++){
                for(size_t col = 0; col < 12; col++){
                    if(M_12_12(row,col) > 1e-12 || M_12_12(row,col) < -1e-12) tripletList.push_back(T(Vertices_dihedral[row/3].getIndex()*3+row%3,Vertices_dihedral[col/3].getIndex()*3+col%3,M_12_12(row,col)));
                }
            }

            
            Positions_edge << geometry->inputVertexPositions[Vertices_edge[0]].x , geometry->inputVertexPositions[Vertices_edge[0]].y , geometry->inputVertexPositions[Vertices_edge[0]].z,
                        geometry->inputVertexPositions[Vertices_edge[1]].x, geometry->inputVertexPositions[Vertices_edge[1]].y, geometry->inputVertexPositions[Vertices_edge[1]].z;

            M_6_6 = constant * Dihedral_angles[e.getIndex()]* geometry->hessian_edge_length(Positions_edge);

            for( size_t row = 0; row < 6; row++){
                for(size_t col = 0; col < 6; col++){
                    tripletList.push_back(T(Vertices_edge[row/3].getIndex()*3+row%3,Vertices_edge[col/3].getIndex()*3+col%3,M_6_6(row,col)));
                }
            }

            
            } //Sumation over adjacent vertices


        // Here we do the last terms


        } //Summation over edges


        Hessian.setFromTriplets(tripletList.begin(),tripletList.end());


        return KB*Hessian;
  

    }

    SparseMatrix<double> E_Handler::H_Bending_2(std::vector<double> Constants) {

    double KB = Constants[0];
    double H0 = Constants[1];

    int N_verts = mesh->nVertices();
    int N_beads = Beads.size();
    // std::cout<<"The number of beads is "<< N_beads <<"\n";
    SparseMatrix<double> Hessian(3*(N_verts+N_beads),3*(N_verts+N_beads));
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

    Eigen::Vector<double,6> Positions_edge;
    Eigen::Vector<double,9> Positions_face;
    Eigen::Vector<double,12> Positions_dihedral;

    Eigen::Vector<double,6> Grad_edge;
    Eigen::Vector<double,9> Grad_face;
    Eigen::Vector<double,12> Grad_dihedral;

    std::array<Vertex,2> Vertices_edge;
    std::array<Vertex,3> Vertices_face;
    std::array<Vertex,4> Vertices_dihedral;

    std::array<Vertex,3 > Vertices_face_2;
    std::array<Vertex,2> Vertices_edge_2;
    std::array<Vertex,4> Vertices_dihedral_2;


    Eigen::Matrix<double,6,6> Hessian_block_edge;
    Eigen::Matrix<double,9,9> Hessian_block_face;
    Eigen::Matrix<double,12,12> Hessian_block_dihedral;

    Halfedge he;
    double constant;

    VertexData<double> Dual_areas(*mesh,0.0);
    VertexData<double> Scalar_MC(*mesh,0.0);
    EdgeData<double> Edge_lengths(*mesh,0.0);
    EdgeData<double> Dihedral_angles(*mesh,0.0);


    for(Vertex v : mesh->vertices()){
        Dual_areas[v] = geometry->barycentricDualArea(v);
        Scalar_MC[v] = geometry->scalarMeanCurvature(v);

    }


    for(Edge e: mesh->edges()){
        Edge_lengths[e] = geometry->edgeLength(e);
        Dihedral_angles[e] = geometry->dihedralAngle(e.halfedge());
    }

    std::vector<Eigen::Vector<double,9>> Gradients_areas;
    std::vector<Eigen::Vector<double,6>> Gradients_edges;
    std::vector<Eigen::Vector<double,12>> Gradients_dihedrals;


    int id= 0;
    for(Face f:mesh -> faces()){
    
        he = f.halfedge();

        Vertices_face[0] = he.vertex();
        Vertices_face[1] = he.next().vertex();
        Vertices_face[2] = he.next().next().vertex();
        Positions_face <<   geometry->inputVertexPositions[Vertices_face[0]].x, geometry->inputVertexPositions[Vertices_face[0]].y, geometry->inputVertexPositions[Vertices_face[0]].z,
                            geometry->inputVertexPositions[Vertices_face[1]].x, geometry->inputVertexPositions[Vertices_face[1]].y, geometry->inputVertexPositions[Vertices_face[1]].z,
                            geometry->inputVertexPositions[Vertices_face[2]].x, geometry->inputVertexPositions[Vertices_face[2]].y, geometry->inputVertexPositions[Vertices_face[2]].z;
        
        Grad_face = geometry->gradient_triangle_area(Positions_face);
        Gradients_areas.push_back(Grad_face);
        
    }

    double dih1 = 0.0;
    double dih2 = 0.0;
    for(Edge e: mesh->edges()){
        Vertices_edge[0] = e.halfedge().vertex();
        Vertices_edge[1] = e.halfedge().twin().vertex();

        Vertices_dihedral[0] = e.halfedge().vertex();
        Vertices_dihedral[1] = e.halfedge().next().vertex();
        Vertices_dihedral[3] = e.halfedge().next().next().vertex();
        Vertices_dihedral[2] = e.halfedge().twin().next().next().vertex();

        Positions_edge << geometry->inputVertexPositions[Vertices_edge[0]].x , geometry->inputVertexPositions[Vertices_edge[0]].y , geometry->inputVertexPositions[Vertices_edge[0]].z,
                        geometry->inputVertexPositions[Vertices_edge[1]].x, geometry->inputVertexPositions[Vertices_edge[1]].y, geometry->inputVertexPositions[Vertices_edge[1]].z;
        Positions_dihedral << geometry->inputVertexPositions[Vertices_dihedral[0]].x , geometry->inputVertexPositions[Vertices_dihedral[0]].y , geometry->inputVertexPositions[Vertices_dihedral[0]].z,
                            geometry->inputVertexPositions[Vertices_dihedral[1]].x, geometry->inputVertexPositions[Vertices_dihedral[1]].y, geometry->inputVertexPositions[Vertices_dihedral[1]].z,
                            geometry->inputVertexPositions[Vertices_dihedral[2]].x, geometry->inputVertexPositions[Vertices_dihedral[2]].y, geometry->inputVertexPositions[Vertices_dihedral[2]].z,
                            geometry->inputVertexPositions[Vertices_dihedral[3]].x, geometry->inputVertexPositions[Vertices_dihedral[3]].y, geometry->inputVertexPositions[Vertices_dihedral[3]].z;

        Grad_edge = geometry->gradient_edge_length(Positions_edge);
        Gradients_edges.push_back(Grad_edge);

        Grad_dihedral = geometry->gradient_dihedral_angle(Positions_dihedral);
        Gradients_dihedrals.push_back(Grad_dihedral);
        // std::cout<<"The edge " << e.getIndex() << "has a edge grad = " << Gradients_edges[e.getIndex()].transpose() <<"\n and dihedral = " << Gradients_dihedrals[e.getIndex()].transpose()<< "\n";
    }

    // I have calculated all the quantities that i will be needing, this takes a lot of memory but lets hope its worth it ahah

    // Lets iterate over the faces   
    Eigen::Matrix<double,9,12> M_9_12;
    Eigen::Matrix<double,9,6> M_9_6;
    Eigen::Matrix<double,9,9> M_9_9;
    Eigen::Matrix<double,6,6> M_6_6;
    Eigen::Matrix<double,6,12> M_6_12;
    Eigen::Matrix<double,12,6> M_12_6;
    Eigen::Matrix<double,12,12> M_12_12;
    Eigen::Matrix<double,6,9> M_6_9;
    Eigen::Matrix<double,12,9> M_12_9;
    size_t max_id = 0;
    for(Face f : mesh->faces()){

        he = f.halfedge();
        Grad_face = Gradients_areas[f.getIndex()];
            
        Vertices_face[0] = he.vertex();
        Vertices_face[1] = he.next().vertex();
        Vertices_face[2] = he.next().next().vertex();

        if(Vertices_face[0].isBoundary() || Vertices_face[1].isBoundary() || Vertices_face[2].isBoundary() ) continue;

        for(Vertex v: f.adjacentVertices()){
            
            constant = -1*(1.0/6.0)*(Scalar_MC[v.getIndex()]-H0)/(Dual_areas[v.getIndex()]*Dual_areas[v.getIndex()]);
            for(Edge e: v.adjacentEdges()){
                he = e.halfedge();
                Vertices_dihedral[0] = he.vertex();
                Vertices_dihedral[1] = he.next().vertex();
                Vertices_dihedral[3] = he.next().next().vertex();
                Vertices_dihedral[2] = he.twin().next().next().vertex();

                Vertices_edge[0] = he.vertex();
                Vertices_edge[1] = he.next().vertex();

                M_9_12 = constant*Edge_lengths[e.getIndex()]*Grad_face * Gradients_dihedrals[e.getIndex()].transpose();

                for(size_t row = 0; row < 9; row ++){
                    for(size_t col = 0; col < 12; col++){
                        tripletList.push_back(T(Vertices_face[row/3].getIndex()*3+row%3, Vertices_dihedral[col/3].getIndex()*3+col%3, M_9_12(row,col)));
              
                        
                }
                }

                M_9_6 = constant * Dihedral_angles[e]* Grad_face * Gradients_edges[e.getIndex()].transpose();

                for(size_t row = 0; row < 9; row ++){
                    for(size_t col = 0; col < 6; col++){
                        tripletList.push_back(T(Vertices_face[row/3].getIndex()*3+row%3, Vertices_edge[col/3].getIndex()*3+col%3, M_9_6(row,col)));
                    
                    }
                }

                M_6_9 = constant * Dihedral_angles[e]*Gradients_edges[e.getIndex()]*Grad_face.transpose();

                    for(size_t row = 0; row < 6; row++){
                        for(size_t col = 0; col < 9; col++){
                        tripletList.push_back(T(Vertices_edge[row/3].getIndex()*3+row%3,Vertices_face[col/3].getIndex()*3+col%3, M_6_9(row,col)));
                     
                        }
                    }

                M_12_9 = constant * Edge_lengths[e.getIndex()]*Gradients_dihedrals[e.getIndex()] * Grad_face.transpose();

                    for(size_t row = 0; row < 12; row++){
                        for( size_t col = 0; col < 9; col++){
                            tripletList.push_back(T(Vertices_dihedral[row/3].getIndex()*3+row%3,Vertices_face[col/3].getIndex()*3+col%3,M_12_9(row,col)));
                           
                        }
                    }

                




            
            
            
            
            
            
            }

            constant = (2.0/9.0) * (Scalar_MC[v.getIndex()]-H0) * (Scalar_MC[v.getIndex()]-H0)/(Dual_areas[v.getIndex()] * Dual_areas[v.getIndex()] * Dual_areas[v.getIndex()]);
            
            
            for(Face f2: v.adjacentFaces()){
                he = f2.halfedge();
                
                Vertices_face_2[0] = he.vertex();
                Vertices_face_2[1] = he.next().vertex();
                Vertices_face_2[2] = he.next().next().vertex();

                M_9_9 = constant * Grad_face * Gradients_areas[f2.getIndex()].transpose();
                    for(size_t row = 0; row < 9; row++){
                        for(size_t col = 0; col < 9; col++){

                            tripletList.push_back(T(Vertices_face[row/3].getIndex()*3+row%3,Vertices_face_2[col/3].getIndex()*3+col%3, M_9_9(row,col) ));
                            // if(isnan(M_9_9(row,col))) std::cout<<" nan flag 3 \n";
                        }
                    }

            }

            // And last but not least we have the hessian term

            constant = -1*(1.0/3.0)*(Scalar_MC[v.getIndex()] - H0)*(Scalar_MC[v.getIndex()] - H0)/(Dual_areas[v.getIndex()]*Dual_areas[v.getIndex()]);
            
            Positions_face << geometry->inputVertexPositions[Vertices_face[0]].x , geometry->inputVertexPositions[Vertices_face[0]].y , geometry->inputVertexPositions[Vertices_face[0]].z,
                        geometry->inputVertexPositions[Vertices_face[1]].x, geometry->inputVertexPositions[Vertices_face[1]].y, geometry->inputVertexPositions[Vertices_face[1]].z,
                        geometry->inputVertexPositions[Vertices_face[2]].x, geometry->inputVertexPositions[Vertices_face[2]].y, geometry->inputVertexPositions[Vertices_face[2]].z;
        
            M_9_9 = constant * geometry->hessian_triangle_area(Positions_face);

            for(size_t row = 0; row < 9; row++){
                for(size_t col = 0; col < 9; col++){
                    tripletList.push_back(T(Vertices_face[row/3].getIndex()*3+row%3, Vertices_face[col/3].getIndex()*3+col%3, M_9_9(row,col)));
                    // if(isnan(M_9_9(row,col))) std::cout<<" nan flag 4 \n";
                }
            }


        }
    }
        // And those are all the terms at the sum over all the vertices (On paper is term I)

        for(Edge e: mesh->edges()){
            he = e.halfedge();

            Vertices_edge[0] = he.vertex();
            Vertices_edge[1] = he.twin().vertex();
            
            Vertices_dihedral[0] = he.vertex();
            Vertices_dihedral[1] = he.next().vertex();
            Vertices_dihedral[3] = he.next().next().vertex();
            Vertices_dihedral[2] = he.twin().next().next().vertex();

            if(Vertices_edge[0].isBoundary() || Vertices_edge[1].isBoundary() || Vertices_dihedral[2].isBoundary() || Vertices_dihedral[3].isBoundary() ) continue;

            for(Vertex v: e.adjacentVertices()){
                constant = 1.0/(8.0*Dual_areas[v.getIndex()]);
                for( Edge e2 : v.adjacentEdges()){
                    Vertices_edge_2[0] = e2.halfedge().vertex();
                    Vertices_edge_2[1] = e2.halfedge().twin().vertex();

                    M_6_6 = constant*Dihedral_angles[e.getIndex()]*Dihedral_angles[e2.getIndex()]*Gradients_edges[e.getIndex()]*Gradients_edges[e2.getIndex()].transpose();

                    for(size_t row = 0; row < 6; row++){
                        for(size_t col = 0; col < 6; col++){
                            tripletList.push_back(T(Vertices_edge[row/3].getIndex()*3+row%3, Vertices_edge_2[col/3].getIndex()*3+col%3,M_6_6(row,col)));
                            // if(isnan(M_6_6(row,col))) std::cout<<" nan flag 5 \n";
                        }
                    }                    


                    he = e2.halfedge();
                    Vertices_dihedral_2[0] = he.vertex();
                    Vertices_dihedral_2[1] = he.next().vertex();
                    Vertices_dihedral_2[3] = he.next().next().vertex();
                    Vertices_dihedral_2[2] = he.twin().next().next().vertex(); 

                    M_6_12 = constant* Dihedral_angles[e.getIndex()]*Edge_lengths[e2.getIndex()]*Gradients_edges[e.getIndex()]*Gradients_dihedrals[e2.getIndex()].transpose();

                    for(size_t row = 0; row < 6; row++){
                        for(size_t col =0; col < 12; col++){
                            tripletList.push_back(T(Vertices_edge[row/3].getIndex()*3+row%3, Vertices_dihedral_2[col/3].getIndex()*3+col%3,M_6_12(row,col)));
                            // if(isnan(M_6_12(row,col))) std::cout<<" nan flag 6 \n";
                        }
                    }

                     

                    M_12_6 = constant * Edge_lengths[e.getIndex()]*Dihedral_angles[e2.getIndex()]*Gradients_dihedrals[e.getIndex()]*Gradients_edges[e2.getIndex()].transpose();

                    for(size_t row = 0; row < 12; row++){
                        for(size_t col = 0; col < 6; col++){
                            tripletList.push_back(T(Vertices_dihedral[row/3].getIndex()*3+row%3, Vertices_edge_2[col/3].getIndex()*3+col%3, M_12_6(row,col)));
                            // if(isnan(M_12_6(row,col))) std::cout<<" nan flag 7 \n";
                        }
                    }

                    M_12_12 = constant* Edge_lengths[e.getIndex()] * Edge_lengths[e2.getIndex()]* Gradients_dihedrals[e.getIndex()]*Gradients_dihedrals[e2.getIndex()].transpose();

                    for(size_t row = 0; row < 12; row++){
                        for(size_t col = 0; col < 12; col++){
                            tripletList.push_back(T(Vertices_dihedral[row/3].getIndex()*3+row%3,Vertices_dihedral_2[col/3].getIndex()*3+col%3,M_12_12(row,col)));
                            // if(isnan(M_12_12(row,col))) std::cout<<"  8 \n";
                        }
                    }




                }


                // Those are the first terms
                // I MOVED THIS TERMS TO THE OTHER SUMMATION






            constant = (1.0/2.0)*(Scalar_MC[v.getIndex()]-H0)/(Dual_areas[v.getIndex()]) ;

            M_6_12 = constant * Gradients_edges[e.getIndex()]*Gradients_dihedrals[e.getIndex()].transpose();

            for( size_t row = 0; row < 6; row++){
                for(size_t col = 0; col < 12; col++){
                    tripletList.push_back(T(Vertices_edge[row/3].getIndex()*3+row%3,Vertices_dihedral[col/3].getIndex()*3+col%3,M_6_12(row,col)));
            
                }
            }


            M_12_6 = constant * Gradients_dihedrals[e.getIndex()]* Gradients_edges[e.getIndex()].transpose();

            for( size_t row = 0; row < 12; row++){
                for(size_t col = 0; col < 6; col++){
                    tripletList.push_back(T(Vertices_dihedral[row/3].getIndex()*3+row%3,Vertices_edge[col/3].getIndex()*3+col%3,M_12_6(row,col)));
                }
            }

            Positions_dihedral <<   geometry->inputVertexPositions[Vertices_dihedral[0]].x , geometry->inputVertexPositions[Vertices_dihedral[0]].y , geometry->inputVertexPositions[Vertices_dihedral[0]].z,
                                    geometry->inputVertexPositions[Vertices_dihedral[1]].x, geometry->inputVertexPositions[Vertices_dihedral[1]].y, geometry->inputVertexPositions[Vertices_dihedral[1]].z,
                                    geometry->inputVertexPositions[Vertices_dihedral[2]].x, geometry->inputVertexPositions[Vertices_dihedral[2]].y, geometry->inputVertexPositions[Vertices_dihedral[2]].z,
                                    geometry->inputVertexPositions[Vertices_dihedral[3]].x, geometry->inputVertexPositions[Vertices_dihedral[3]].y, geometry->inputVertexPositions[Vertices_dihedral[3]].z;

            


            M_12_12 = constant * Edge_lengths[e.getIndex()]* geometry->hessian_dihedral_angle(Positions_dihedral);

            for( size_t row = 0; row < 12; row++){
                for(size_t col = 0; col < 12; col++){
                    if(M_12_12(row,col) > 1e-12 || M_12_12(row,col) < -1e-12) tripletList.push_back(T(Vertices_dihedral[row/3].getIndex()*3+row%3,Vertices_dihedral[col/3].getIndex()*3+col%3,M_12_12(row,col)));
                }
            }

            
            Positions_edge << geometry->inputVertexPositions[Vertices_edge[0]].x , geometry->inputVertexPositions[Vertices_edge[0]].y , geometry->inputVertexPositions[Vertices_edge[0]].z,
                        geometry->inputVertexPositions[Vertices_edge[1]].x, geometry->inputVertexPositions[Vertices_edge[1]].y, geometry->inputVertexPositions[Vertices_edge[1]].z;

            M_6_6 = constant * Dihedral_angles[e.getIndex()]* geometry->hessian_edge_length(Positions_edge);

            for( size_t row = 0; row < 6; row++){
                for(size_t col = 0; col < 6; col++){
                    tripletList.push_back(T(Vertices_edge[row/3].getIndex()*3+row%3,Vertices_edge[col/3].getIndex()*3+col%3,M_6_6(row,col)));
                }
            }

            
            } //Sumation over adjacent vertices


        // Here we do the last terms


        } //Summation over edges


        Hessian.setFromTriplets(tripletList.begin(),tripletList.end());


        return KB*Hessian;
  

    }


SparseMatrix<double> E_Handler::H_Volume(std::vector<double> Constants){
    

    // std::cout<<"Hessian volume\n";
    double KV = Constants[0];

    int N_verts = mesh->nVertices();
    int N_beads = Beads.size();
    SparseMatrix<double> Hessian(3*(N_verts+N_beads),3*(N_verts+N_beads));
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

    Eigen::Matrix<double,9,9> Hessian_block_vol;
    Eigen::Vector<double, 9> Positions;
    std::vector<Vertex> Vertices(3);
    Halfedge he;
    Eigen::Matrix3d Zeros = Eigen::Matrix3d::Zero();



    for(Face f: mesh->faces()){
        he = f.halfedge();
        Vertices[0] = he.vertex();
        Vertices[1] = he.next().vertex();
        Vertices[2] = he.next().next().vertex();

         Positions << geometry->inputVertexPositions[Vertices[0]].x , geometry->inputVertexPositions[Vertices[0]].y , geometry->inputVertexPositions[Vertices[0]].z,
                    geometry->inputVertexPositions[Vertices[1]].x, geometry->inputVertexPositions[Vertices[1]].y, geometry->inputVertexPositions[Vertices[1]].z,
                    geometry->inputVertexPositions[Vertices[2]].x, geometry->inputVertexPositions[Vertices[2]].y, geometry->inputVertexPositions[Vertices[2]].z;

        // Now we do the hesian thingy 
        Hessian_block_vol = geometry->hessian_volume(Positions);
        
        for(size_t row = 0; row < 9; row++){
            for(size_t col = 0; col< 9; col++ ){
                if( Hessian_block_vol(row,col) > 1e-12 || Hessian_block_vol(row,col) < -1e-12) tripletList.push_back(T(Vertices[row/3].getIndex()*3+row%3,Vertices[col/3].getIndex()*3+col%3,Hessian_block_vol(row,col) ));

            }
        }


    }
    Hessian.setFromTriplets(tripletList.begin(),tripletList.end());



    // std::cout<<"Hessian volume done\n";
    return KV*Hessian;
}

SparseMatrix<double> E_Handler::H_Laplace(std::vector<double> Constants){
    double KB = Constants[0];
    double H0 = Constants[1];
    
    int N_verts = mesh->nVertices();
    int N_beads = Beads.size();
    SparseMatrix<double> Hessian(3*(N_verts+N_beads),3*(N_verts+N_beads));
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

    
    Eigen::Vector<double,3> dX;
    Eigen::Vector<double,3> dX2;
    
    // Eigen::Vector<double,6> Positions_edge;
    Eigen::Vector<double, 9> Positions_face;
    Eigen::Vector<double,12> Positions_cotan;

    Eigen::Vector<double,9> Grad_face;
    Eigen::Vector<double,12> Grad_cotan_sample;
    Eigen::Vector<double,12> Grad_cotan;
    Eigen::Vector<double,12> Grad_cotan2;



    std::array<Vertex,2> Vertices_edge;
    std::array<Vertex,3> Vertices_face;
    std::array<Vertex,4> Vertices_cotan;

    std::array<Vertex,3 > Vertices_face_2;
    std::array<Vertex,2> Vertices_edge_2;
    std::array<Vertex,4> Vertices_cotan_2;


    Halfedge he;
    double constant;

    VertexData<double> Dual_areas(*mesh,0.0);
    VertexData<double> Hi(*mesh,0.0);
    EdgeData<double> Cotan_W(*mesh,0.0);
    VertexData<Vector3> Q(*mesh,{0.0,0.0,0.0});
    std::vector<Eigen::Vector<double,3>> Qs;
    
    // std::cout<<"3\n";
    
    
    for(Vertex v: mesh->vertices()){

        Dual_areas[v] = geometry->barycentricDualArea(v);
        Q[v] = geometry->vertexNormalMeanCurvature(v);
        Qs.push_back( {Q[v].x,Q[v].y,Q[v].z});
        Hi[v] = Q[v].norm();
        
    }

    // std::cout<<"4\n";
    
    std::vector<Eigen::Vector<double,9>> Gradients_areas;
    std::vector<Eigen::Vector<double,12>> Gradients_cotan;

    for(Face f: mesh->faces()){
        he = f.halfedge();

        Vertices_face[0] = he.vertex();
        Vertices_face[1] = he.next().vertex();
        Vertices_face[2] = he.next().next().vertex();
        Positions_face <<   geometry->inputVertexPositions[Vertices_face[0]].x, geometry->inputVertexPositions[Vertices_face[0]].y, geometry->inputVertexPositions[Vertices_face[0]].z,
                            geometry->inputVertexPositions[Vertices_face[1]].x, geometry->inputVertexPositions[Vertices_face[1]].y, geometry->inputVertexPositions[Vertices_face[1]].z,
                            geometry->inputVertexPositions[Vertices_face[2]].x, geometry->inputVertexPositions[Vertices_face[2]].y, geometry->inputVertexPositions[Vertices_face[2]].z;
        
        Grad_face = geometry->gradient_triangle_area(Positions_face);
        Gradients_areas.push_back(Grad_face);
        
         
    }

    // std::cout<<"5\n";
    
    // std::cout<<"ITERATING OVER EDGES \n ";
    for(Edge e: mesh->edges()){
        // std::cout << e.getIndex()<< " ";

        he = e.halfedge();

        Vertices_cotan[0] = he.vertex();
        Vertices_cotan[1] = he.next().vertex();
        Vertices_cotan[2] = he.next().next().vertex();
        Vertices_cotan[3] = he.twin().next().next().vertex();

        Positions_cotan <<  geometry->inputVertexPositions[Vertices_cotan[0]].x, geometry->inputVertexPositions[Vertices_cotan[0]].y, geometry->inputVertexPositions[Vertices_cotan[0]].z,
                            geometry->inputVertexPositions[Vertices_cotan[1]].x, geometry->inputVertexPositions[Vertices_cotan[1]].y, geometry->inputVertexPositions[Vertices_cotan[1]].z,
                            geometry->inputVertexPositions[Vertices_cotan[2]].x, geometry->inputVertexPositions[Vertices_cotan[2]].y, geometry->inputVertexPositions[Vertices_cotan[2]].z,
                            geometry->inputVertexPositions[Vertices_cotan[3]].x, geometry->inputVertexPositions[Vertices_cotan[3]].y, geometry->inputVertexPositions[Vertices_cotan[3]].z;
        
        Grad_cotan_sample = geometry->gradient_cotan_weight(Positions_cotan);
        // std::cout<<"Grad cotan is" << Grad_cotan_sample.transpose() << "\n";
        Gradients_cotan.push_back(Grad_cotan_sample);
        Cotan_W[e] = geometry->Cotan_weight( Positions_cotan);

    }
    // std::cout<<"\n";
    // std::cout<<"Gradients cotan at edge 0 is " << Gradients_cotan[0] << "\n";
    // std::cout<<"6\n";
    
    Eigen::Matrix<double,9,9> M_9_9;
    Eigen::Matrix<double,9,12> M_9_12;
    // Eigen::Matrix<double,12,9> M_12_9;
    // Eigen::Matrix<double,9,6> M_9_6;
    // Eigen::Matrix<double,6,9> M_6_9;
    Eigen::Matrix<double,12,12> M_12_12;
    // Eigen::Matrix<double,12,6> M_12_6;
    // Eigen::Matrix<double,6,12> M_6_12;
    Eigen::Matrix<double,9,3> M_9_3;
    // Eigen::Matrix<double,3,9> M_3_9;
    Eigen::Matrix<double,3,12> M_3_12;
    Eigen::Matrix<double,12,3> M_12_3;
    
    Eigen::Matrix<double,3,3> M_3_3;
    Eigen::Matrix<double,3,3> Id = Eigen::Matrix3d::Identity();
    

    
    for(Vertex v: mesh->vertices()){
        // std::cout<<"v.getIndex() = " << v.getIndex() << "\n";
        for(Face f: v.adjacentFaces()){
            // std::cout<<"f.getIndex() = " << f.getIndex() << "\n";
            he = f.halfedge();
            Grad_face = Gradients_areas[f.getIndex()];
            Vertices_face[0] = he.vertex();
            Vertices_face[1] = he.next().vertex();
            Vertices_face[2] = he.next().next().vertex();
            Positions_face << geometry->inputVertexPositions[Vertices_face[0]].x , geometry->inputVertexPositions[Vertices_face[0]].y , geometry->inputVertexPositions[Vertices_face[0]].z,
                            geometry->inputVertexPositions[Vertices_face[1]].x, geometry->inputVertexPositions[Vertices_face[1]].y, geometry->inputVertexPositions[Vertices_face[1]].z,
                            geometry->inputVertexPositions[Vertices_face[2]].x, geometry->inputVertexPositions[Vertices_face[2]].y, geometry->inputVertexPositions[Vertices_face[2]].z;
        
            // std::cout<<"1\n";
    

            if(Vertices_face[0].isBoundary() || Vertices_face[1].isBoundary() || Vertices_face[2].isBoundary() ) continue;
            // std::cout<<"2\n";
    

            constant = -1*(1.0/3.0)*( Hi[v]*Hi[v]/(Dual_areas[v]*Dual_areas[v])-H0*H0 );

          
            M_9_9 = constant * geometry->hessian_triangle_area(Positions_face);

            for(size_t row = 0; row < 9; row++){
                for(size_t col = 0; col < 9; col++){
                    tripletList.push_back(T(Vertices_face[row/3].getIndex()*3+row%3, Vertices_face[col/3].getIndex()*3+col%3, M_9_9(row,col)));
                    // if(isnan(M_9_9(row,col))) std::cout<<" nan flag 4 \n";
                }
            }
            // std::cout<<"3\n";
    
            constant = (2.0/9.0)*Hi[v]*Hi[v]/(Dual_areas[v]*Dual_areas[v]*Dual_areas[v]);
            
            for(Face f2: v.adjacentFaces()){
                he = f2.halfedge();

                Vertices_face_2[0] = he.vertex();
                Vertices_face_2[1] = he.next().vertex();
                Vertices_face_2[2] = he.next().next().vertex();

                M_9_9 = constant * Grad_face * Gradients_areas[f2.getIndex()].transpose();
                for(size_t row = 0; row < 9; row++){
                    for(size_t col = 0; col < 9; col++){

                        tripletList.push_back(T(Vertices_face[row/3].getIndex()*3+row%3,Vertices_face_2[col/3].getIndex()*3+col%3, M_9_9(row,col) ));
                        // if(isnan(M_9_9(row,col))) std::cout<<" nan flag 3 \n";
                    }
                }

            }


             //Here i am in a face in a vertex
            // std::cout<<"4\n";
    
            constant = -(1.0/6.0)*(1/(Dual_areas[v]*Dual_areas[v]));
            for(Halfedge he1: v.outgoingHalfedges()){
                // he = e.halfedge();
                
                // if(he1.vertex().getIndex() != v.getIndex()) he1 = he.twin();
                Vertices_edge[0] = he1.vertex();
                Vertices_edge[1] = he1.twin().vertex();
                
                dX <<   geometry->inputVertexPositions[Vertices_edge[0]].x - geometry->inputVertexPositions[Vertices_edge[1]].x,
                        geometry->inputVertexPositions[Vertices_edge[0]].y - geometry->inputVertexPositions[Vertices_edge[1]].y,
                        geometry->inputVertexPositions[Vertices_edge[0]].z - geometry->inputVertexPositions[Vertices_edge[1]].z;
                 // dX = geometry->inputVertexPositions[]
    
                Vertices_cotan[0] = he1.edge().halfedge().vertex();
                Vertices_cotan[1] = he1.edge().halfedge().next().vertex();
                Vertices_cotan[2] = he1.edge().halfedge().next().next().vertex();
                Vertices_cotan[3] = he1.edge().halfedge().twin().next().next().vertex();

    
                 // THere are two terms, one is the gradient of the weight and the other gradient
                // std::cout<<"The he1 edge index is " << he1.edge().getIndex() << "\n";
                // std::cout<<"The size of Gradients cotan is " << Gradients_cotan.size() << "\n";
                // std::cout<<"Gradients cotan at 0 is" << Gradients_cotan[0] << "\n";
                Grad_cotan = Eigen::Vector<double,12>(Gradients_cotan[he1.edge().getIndex()]);

                // std::cout<<"5.5\n";
                M_9_12 = constant * Grad_face * Grad_cotan.transpose()*(Qs[v.getIndex()].transpose()*dX) ;

                for(size_t row = 0; row < 9; row++){
                    for(size_t col = 0; col < 12; col++){
                        tripletList.push_back(T(Vertices_face[row/3].getIndex()*3+row%3, Vertices_cotan[col/3].getIndex()*3+col%3, M_9_12(row,col)  ));
                        tripletList.push_back(T(Vertices_cotan[col/3].getIndex()*3+col%3, Vertices_face[row/3].getIndex()*3+row%3, M_9_12(row,col)  ));
                    }
                }
    

                M_9_3 = constant * Grad_face* Qs[v.getIndex()].transpose()*Cotan_W[he1.edge()];
                 // Now this term is positive por vertex 1 and negative for vertex 0
                for(size_t row = 0; row < 9; row++){
                    for(size_t col = 0; col < 3; col++){
                        tripletList.push_back(T(Vertices_face[row/3].getIndex()*3+row%3, Vertices_edge[0].getIndex()*3+col%3 ,    M_9_3(row,col) ));
                        tripletList.push_back(T(Vertices_face[row/3].getIndex()*3+row%3, Vertices_edge[1].getIndex()*3+col%3 ,  -1*M_9_3(row,col) ));
                        
                        tripletList.push_back(T(Vertices_edge[0].getIndex()*3+col%3, Vertices_face[row/3].getIndex()*3+row%3, M_9_3(row,col)  ));
                        tripletList.push_back(T(Vertices_edge[1].getIndex()*3+col%3, Vertices_face[row/3].getIndex()*3+row%3, -1*M_9_3(row,col)  ));

                    }
                }


            }

        

        }


    }
    

    // // AHora faltan los terminos que solo dependen de V Y E;
    for(Vertex v: mesh->vertices()){
        constant = 1/(8.0*Dual_areas[v]*Hi[v]*Hi[v]);

        for(Halfedge he1: v.outgoingHalfedges()){
            Vertices_edge[0] = he1.vertex();
            Vertices_edge[1] = he1.twin().vertex();

            
            dX <<   geometry->inputVertexPositions[Vertices_edge[0]].x - geometry->inputVertexPositions[Vertices_edge[1]].x,
                    geometry->inputVertexPositions[Vertices_edge[0]].y - geometry->inputVertexPositions[Vertices_edge[1]].y,
                    geometry->inputVertexPositions[Vertices_edge[0]].z - geometry->inputVertexPositions[Vertices_edge[1]].z;

            Vertices_cotan[0] = he1.edge().halfedge().vertex();
            Vertices_cotan[1] = he1.edge().halfedge().next().vertex();
            Vertices_cotan[2] = he1.edge().halfedge().next().next().vertex();
            Vertices_cotan[3] = he1.edge().halfedge().twin().next().next().vertex();
            Grad_cotan = Eigen::Vector<double,12>(Gradients_cotan[he1.edge().getIndex()]);
            

            for(Halfedge he2: v.outgoingHalfedges()){
            Vertices_edge_2[0] = he2.vertex();
            Vertices_edge_2[1] = he2.twin().vertex();
            dX2 <<  geometry->inputVertexPositions[Vertices_edge_2[0]].x - geometry->inputVertexPositions[Vertices_edge_2[1]].x,
                    geometry->inputVertexPositions[Vertices_edge_2[0]].y - geometry->inputVertexPositions[Vertices_edge_2[1]].y,
                    geometry->inputVertexPositions[Vertices_edge_2[0]].z - geometry->inputVertexPositions[Vertices_edge_2[1]].z;
                
            Vertices_cotan_2[0] = he2.edge().halfedge().vertex();
            Vertices_cotan_2[1] = he2.edge().halfedge().next().vertex();
            Vertices_cotan_2[2] = he2.edge().halfedge().next().next().vertex();
            Vertices_cotan_2[3] = he2.edge().halfedge().twin().next().next().vertex();
            Grad_cotan2 = Eigen::Vector<double,12>(Gradients_cotan[he2.edge().getIndex()]);

    //         // 
            M_12_12 = (constant)* Grad_cotan * Grad_cotan2.transpose()*( dX2.transpose()*Qs[v.getIndex()] )*(Qs[v.getIndex()].transpose()*dX) ;

            for(int row = 0; row < 12; row++){
                for(int col = 0; col < 12; col++){
                    tripletList.push_back(T(Vertices_cotan[row/3].getIndex()*3+row%3, Vertices_cotan_2[col/3].getIndex()*3 + col%3, M_12_12(row,col) ));    
                }
            }
    //         // 
            M_3_12 =  (constant*Cotan_W[he1.edge()])*Qs[v.getIndex()]*Grad_cotan2.transpose()*(Qs[v.getIndex()].transpose()*dX2 );

            for(int row = 0; row < 3; row++){
                for(int col = 0; col < 12; col++){
                    tripletList.push_back(T( Vertices_edge[0].getIndex()*3+row%3, Vertices_cotan_2[col/3].getIndex()*3+col%3, M_3_12(row,col)));
                    tripletList.push_back(T( Vertices_edge[1].getIndex()*3+row%3, Vertices_cotan_2[col/3].getIndex()*3+col%3, -1*M_3_12(row,col)));
                }

            }//Creo que esta correcto pero podemos revisitar

            // M_12_3 (constant*Cotan_W[he2.edge()] ) *Grad_cotan * Qs[v.getIndex()].transpose() ;// *  (Qs[v.getIndex()].transpose()*dX);
            M_12_3 = (constant*Cotan_W[he2.edge()])*Grad_cotan * Qs[v.getIndex()].transpose() * (Qs[v.getIndex()].transpose()*dX );

            for(int row = 0; row <12; row ++){
                for(int col = 0; col < 3; col++){
                    tripletList.push_back(T( Vertices_cotan[row/3].getIndex()*3+row%3, Vertices_edge_2[0].getIndex()*3+col%3, M_12_3(row,col)));
                    tripletList.push_back(T( Vertices_cotan[row/3].getIndex()*3+row%3, Vertices_edge_2[1].getIndex()*3+col%3, -1*M_12_3(row,col)));                    
                }
            }


            M_3_3 = (constant*Cotan_W[he1.edge()]*Cotan_W[he2.edge()])*Qs[v.getIndex()]*Qs[v.getIndex()].transpose();

            for(int row = 0; row < 3; row++){
                for(int col = 0; col < 3; col++){
                    tripletList.push_back(T( Vertices_edge[0].getIndex()*3+row%3, Vertices_edge_2[0].getIndex()*3+col%3,  M_3_3(row,col) ));
                    tripletList.push_back(T( Vertices_edge[1].getIndex()*3+row%3, Vertices_edge_2[1].getIndex()*3+col%3,  M_3_3(row,col) ));
                    tripletList.push_back(T( Vertices_edge[0].getIndex()*3+row%3, Vertices_edge_2[1].getIndex()*3+col%3,  -1*M_3_3(row,col) ));
                    tripletList.push_back(T( Vertices_edge[1].getIndex()*3+row%3, Vertices_edge_2[0].getIndex()*3+col%3,   -1*M_3_3(row,col) ));                
                }
            }
            
            }


        }


        // Here i am again at the vertex
        for(Halfedge he1 : v.outgoingHalfedges()){
            Vertices_edge[0] = he1.vertex();
            Vertices_edge[1] = he1.twin().vertex();

            
            dX <<   geometry->inputVertexPositions[Vertices_edge[0]].x-geometry->inputVertexPositions[Vertices_edge[1]].x,
                    geometry->inputVertexPositions[Vertices_edge[0]].y-geometry->inputVertexPositions[Vertices_edge[1]].y,
                    geometry->inputVertexPositions[Vertices_edge[0]].z-geometry->inputVertexPositions[Vertices_edge[1]].z;

            Vertices_cotan[0] = he1.edge().halfedge().vertex();
            Vertices_cotan[1] = he1.edge().halfedge().next().vertex();
            Vertices_cotan[2] = he1.edge().halfedge().next().next().vertex();
            Vertices_cotan[3] = he1.edge().halfedge().twin().next().next().vertex();
            Grad_cotan = Eigen::Vector<double,12>(Gradients_cotan[he1.edge().getIndex()]);
            
            constant = 2*(Hi[v]/Dual_areas[v] - H0)/(16.0*Hi[v]);
            for(Halfedge he2 : v.outgoingHalfedges()){

                Vertices_edge_2[0] = he2.vertex();
                Vertices_edge_2[1] = he2.twin().vertex();
                dX2 <<   geometry->inputVertexPositions[Vertices_edge_2[0]].x-geometry->inputVertexPositions[Vertices_edge_2[1]].x,
                                    geometry->inputVertexPositions[Vertices_edge_2[0]].y-geometry->inputVertexPositions[Vertices_edge_2[1]].y,
                                    geometry->inputVertexPositions[Vertices_edge_2[0]].z-geometry->inputVertexPositions[Vertices_edge_2[1]].z;
                    
                Vertices_cotan_2[0] = he2.edge().halfedge().vertex();
                Vertices_cotan_2[1] = he2.edge().halfedge().next().vertex();
                Vertices_cotan_2[2] = he2.edge().halfedge().next().next().vertex();
                Vertices_cotan_2[3] = he2.edge().halfedge().twin().next().next().vertex();
                Grad_cotan2 = Eigen::Vector<double,12>(Gradients_cotan[he2.edge().getIndex()]);

                M_12_12 = constant * Grad_cotan * Grad_cotan2.transpose() * ( dX.transpose()*(  Id - (Qs[v.getIndex()]*Qs[v.getIndex()].transpose())/(Hi[v]*Hi[v]) )*dX2);

                for(int row = 0; row < 12; row++){
                    for(int col = 0; col < 12; col++){
                        tripletList.push_back(T(Vertices_cotan[row/3].getIndex()*3+row%3, Vertices_cotan_2[col/3].getIndex()*3 + col%3, M_12_12(row,col) ));    
                    }
                }

                M_3_12 =  constant * Cotan_W[he1.edge()]*Id*( Id - (Qs[v.getIndex()]*Qs[v.getIndex()].transpose())/(Hi[v]*Hi[v]) )*dX2*Grad_cotan2.transpose();

                for(int row = 0; row < 3; row++){
                    for(int col = 0; col < 12; col++){
                        tripletList.push_back(T( Vertices_edge[0].getIndex()*3+row%3, Vertices_cotan_2[col/3].getIndex()*3+col%3, M_3_12(row,col)));
                        tripletList.push_back(T( Vertices_edge[1].getIndex()*3+row%3, Vertices_cotan_2[col/3].getIndex()*3+col%3, -1*M_3_12(row,col)));
                    }
                }

                M_12_3 = constant * Cotan_W[he2.edge()]* Grad_cotan* dX.transpose()*( Id - (Qs[v.getIndex()]*Qs[v.getIndex()].transpose())/(Hi[v]*Hi[v]) ); 
                
                for(int row = 0; row < 12; row++){
                    for(int col = 0; col < 3; col++){
                        tripletList.push_back(T( Vertices_cotan[row/3].getIndex()*3+row%3, Vertices_edge_2[0].getIndex()*3+col%3, M_12_3(row,col)));
                        tripletList.push_back(T( Vertices_cotan[row/3].getIndex()*3+row%3, Vertices_edge_2[1].getIndex()*3+col%3, -1*M_12_3(row,col)));
                    }
                }


                M_3_3 = constant* Cotan_W[he1.edge()]*Cotan_W[he2.edge()]* ( Id - (Qs[v.getIndex()]*Qs[v.getIndex()].transpose())/(Hi[v]*Hi[v]) );

                for(int row = 0; row < 3; row++){
                    for(int col = 0; col < 3; col++){
                        tripletList.push_back(T( Vertices_edge[0].getIndex()*3+row%3, Vertices_edge_2[0].getIndex()*3+col%3,    M_3_3(row,col) ));
                        tripletList.push_back(T( Vertices_edge[1].getIndex()*3+row%3, Vertices_edge_2[1].getIndex()*3+col%3,    M_3_3(row,col) ));
                        tripletList.push_back(T( Vertices_edge[0].getIndex()*3+row%3, Vertices_edge_2[1].getIndex()*3+col%3, -1*M_3_3(row,col) ));
                        tripletList.push_back(T( Vertices_edge[1].getIndex()*3+row%3, Vertices_edge_2[0].getIndex()*3+col%3, -1*M_3_3(row,col) ));                
                    }
                }


            }

            // And now finally the Hessian term;

            // 
            constant = 2*(Hi[v]/Dual_areas[v] - H0)/(4.0*Hi[v]);

            M_3_12 = constant * Qs[v.getIndex()]*Grad_cotan.transpose();
            
            for(int row = 0; row < 3; row++){
                    for(int col = 0; col < 12; col++){
                        tripletList.push_back(T( Vertices_edge[0].getIndex()*3+row%3, Vertices_cotan[col/3].getIndex()*3+col%3, M_3_12(row,col)));
                        tripletList.push_back(T( Vertices_cotan[col/3].getIndex()*3+col%3, Vertices_edge[0].getIndex()*3+row%3, M_3_12(row,col)));
                        
                        tripletList.push_back(T( Vertices_edge[1].getIndex()*3+row%3, Vertices_cotan[col/3].getIndex()*3+col%3,  -1*M_3_12(row,col)));
                        tripletList.push_back(T( Vertices_cotan[col/3].getIndex()*3+col%3, Vertices_edge[1].getIndex()*3+row%3,   -1*M_3_12(row,col)));                
                    }
                }

                        // Out last one is 
            Positions_cotan <<  geometry->inputVertexPositions[Vertices_cotan[0]].x, geometry->inputVertexPositions[Vertices_cotan[0]].y, geometry->inputVertexPositions[Vertices_cotan[0]].z,
                                geometry->inputVertexPositions[Vertices_cotan[1]].x, geometry->inputVertexPositions[Vertices_cotan[1]].y, geometry->inputVertexPositions[Vertices_cotan[1]].z,
                                geometry->inputVertexPositions[Vertices_cotan[2]].x, geometry->inputVertexPositions[Vertices_cotan[2]].y, geometry->inputVertexPositions[Vertices_cotan[2]].z,
                                geometry->inputVertexPositions[Vertices_cotan[3]].x, geometry->inputVertexPositions[Vertices_cotan[3]].y, geometry->inputVertexPositions[Vertices_cotan[3]].z;
        
            M_12_12 = constant*geometry->hessian_cotan_weight(Positions_cotan)*(Qs[v.getIndex()].transpose()*dX );

            for(int row = 0; row < 12; row++){
                for(int col = 0; col < 12; col++){
                    // if(fabs(M_12_12(row,col))>1e-10)
                    tripletList.push_back(T(Vertices_cotan[row/3].getIndex()*3+row%3, Vertices_cotan[col/3].getIndex()*3 + col%3, M_12_12(row,col)  ));    
                }
            }
            

        }




    }
    // std::cout<<"DOne with hessian\n";
    Hessian.setFromTriplets(tripletList.begin(),tripletList.end());
    
    
    return KB*Hessian;
}


SparseMatrix<double> E_Handler::H_Edge_reg(std::vector<double> Constants){
    double KE = Constants[0];
    // std::cout<<"Doing edge reg\n";
    int N_verts = mesh->nVertices();
    int N_beads = Beads.size();
    SparseMatrix<double> Hessian(3*(N_verts+N_beads),3*(N_verts+N_beads));
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    Eigen::Matrix<double,9,9> Hessian_block_edge;
    Eigen::Vector<double, 9> Positions;
    std::vector<Vertex> Vertices(3);
    Halfedge he;
    Eigen::Matrix3d Zeros = Eigen::Matrix3d::Zero();
    EdgeData<double> Edge_lengths(*mesh, 0.0);
    Eigen::Vector<double, 3> Edge_lengths_prev;

    for(Edge e: mesh->edges()){
        Edge_lengths[e] = geometry->edgeLength(e);
        // std::cout<<"" << e.getIndex() << " and its length is " << Edge_lengths[e] << "\n";
    }

    for(Face f: mesh->faces()){
        he = f.halfedge();

        Vertices[0] = he.vertex();
        Vertices[1] = he.next().vertex();
        Vertices[2] = he.next().next().vertex();

        if(Vertices[0].isBoundary() || Vertices[1].isBoundary() || Vertices[2].isBoundary() ) continue;

        Positions << geometry->inputVertexPositions[Vertices[0]].x , geometry->inputVertexPositions[Vertices[0]].y , geometry->inputVertexPositions[Vertices[0]].z,
                    geometry->inputVertexPositions[Vertices[1]].x, geometry->inputVertexPositions[Vertices[1]].y, geometry->inputVertexPositions[Vertices[1]].z,
                    geometry->inputVertexPositions[Vertices[2]].x, geometry->inputVertexPositions[Vertices[2]].y, geometry->inputVertexPositions[Vertices[2]].z;
        Edge_lengths_prev << Edge_lengths[he.edge()], 
                            Edge_lengths[he.next().next().edge()], 
                            Edge_lengths[he.twin().next().edge()];

        // Now we do the hesian thingy
        Hessian_block_edge = geometry->hessian_edge_regular(Positions, Edge_lengths_prev);

        // std::cout<<"block calculated\n";
        for(size_t row = 0; row < 9; row++){
            for(size_t col = 0; col< 9; col++ ){
                if( Hessian_block_edge(row,col) > 1e-12 || Hessian_block_edge(row,col) < -1e-12) tripletList.push_back(T(Vertices[row/3].getIndex()*3+row%3,Vertices[col/3].getIndex()*3+col%3,Hessian_block_edge(row,col) ));

            }
        }
    }

    // std::cout<<"Setting balues\n";
    Hessian.setFromTriplets(tripletList.begin(),tripletList.end());


    return KE*Hessian;



}   

SparseMatrix<double> E_Handler::H_Edge_reg_2(std::vector<double> Constants){

    // Lets get this hessian
    double KE = Constants[0];

    int N_verts = mesh->nVertices();
    int N_beads = Beads.size();
    SparseMatrix<double> Hessian(3*(N_verts+N_beads),3*(N_verts+N_beads));
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

    Eigen::Matrix<double,12,12> Hessian_block_edge;
    Eigen::Vector<double, 12> Positions;
    std::vector<Vertex> Vertices(4);
    Halfedge he;
    Eigen::Matrix3d Zeros = Eigen::Matrix3d::Zero();

    EdgeData<double> Edge_lengths(*mesh, 0.0);

    Eigen::Vector<double, 5> Edge_lengths_prev;

    for(Edge e: mesh->edges()){
        Edge_lengths[e] = geometry->edgeLength(e);
        // std::cout<<"" << e.getIndex() << " and its length is " << Edge_lengths[e] << "\n";
    }

    for(Edge e: mesh->edges()){
        he = e.halfedge();
        Vertices[0] = he.vertex();
        Vertices[1] = he.next().vertex();
        Vertices[2] = he.next().next().vertex();
        Vertices[3] = he.twin().next().next().vertex();

        Positions << geometry->inputVertexPositions[Vertices[0]].x , geometry->inputVertexPositions[Vertices[0]].y , geometry->inputVertexPositions[Vertices[0]].z,
                    geometry->inputVertexPositions[Vertices[1]].x, geometry->inputVertexPositions[Vertices[1]].y, geometry->inputVertexPositions[Vertices[1]].z,
                    geometry->inputVertexPositions[Vertices[2]].x, geometry->inputVertexPositions[Vertices[2]].y, geometry->inputVertexPositions[Vertices[2]].z,
                    geometry->inputVertexPositions[Vertices[3]].x, geometry->inputVertexPositions[Vertices[3]].y, geometry->inputVertexPositions[Vertices[3]].z;

        Edge_lengths_prev << Edge_lengths[he.edge()], 
                            Edge_lengths[he.next().next().edge()], 
                            Edge_lengths[he.twin().next().edge()],
                            Edge_lengths[he.next().edge()],
                            Edge_lengths[he.twin().next().next().edge()];
        // Now we do the hesian thingy 
        Hessian_block_edge = geometry->hessian_edge_regular(Positions, Edge_lengths_prev);
        // std::cout<<"block calculated\n";
        for(size_t row = 0; row < 12; row++){
            for(size_t col = 0; col< 12; col++ ){
                if( Hessian_block_edge(row,col) > 1e-12 || Hessian_block_edge(row,col) < -1e-12) tripletList.push_back(T(Vertices[row/3].getIndex()*3+row%3,Vertices[col/3].getIndex()*3+col%3,Hessian_block_edge(row,col) ));

            }
        }
    }
        
    Hessian.setFromTriplets(tripletList.begin(),tripletList.end());

    return KE*Hessian;

}






SparseMatrix<double> E_Handler::H_Face_reg(std::vector<double> Constants){

    // Ok so this functino will assemble the Hessi an for the surface tension energy 
    int nVerts = mesh->nVertices();
    int nBeads = Beads.size();
    double KA = Constants[0];
    // Eigen::MatrixXd Hessian = Eigen::MatrixXd::Zero(nVerts,nVerts);
    SparseMatrix<double> Hessian(3*(nVerts+nBeads),3*(nVerts+nBeads) );
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    // Eigen::Matrix<double, nVerts, nVerts> Hessian;
    Vertex v1;
    Vertex v2;
    Vertex v3;
    Halfedge he;
    size_t index1;
    size_t index2;
    size_t index3;
    array<int,3> indices;
    // 
    for(Face f: mesh->faces()){
        // 
        he = f.halfedge();
        v1 = he.vertex();
        v2 = he.next().vertex();
        v3 = he.next().next().vertex();
        indices[0] = v1.getIndex();
        indices[1] = v2.getIndex();
        indices[2] = v3.getIndex();
        
        // So we will add a modification for boundaries then check if it works :p
        if(v1.isBoundary() || v2.isBoundary() || v3.isBoundary() ) continue;

        Eigen::Vector<double,9> Positions;

        Positions << geometry->inputVertexPositions[v1].x,geometry->inputVertexPositions[v1].y,geometry->inputVertexPositions[v1].z,
                    geometry->inputVertexPositions[v2].x,geometry->inputVertexPositions[v2].y,geometry->inputVertexPositions[v2].z,
                    geometry->inputVertexPositions[v3].x,geometry->inputVertexPositions[v3].y,geometry->inputVertexPositions[v3].z;
         
        Eigen::Matrix<double,9,9> Hessian_block = KA*((geometry->faceArea(f)-Face_reference[f])/(Face_reference[f]*Face_reference[f])) *geometry->hessian_triangle_area(Positions);
        // I need to multiply the hessian block by the constant

        Eigen::Vector<double,9> Gradient_face = geometry->gradient_triangle_area(Positions);
        Eigen::Matrix<double,9,9> Hessian_block_2 = (KA/(Face_reference[f]*Face_reference[f])) * Gradient_face * Gradient_face.transpose();
        // Now i need to load this quantities onto a bigger matrix ...
        for(int i = 0; i < 9; i++){
            for(int j = 0; j <9; j++){

                // So i am at vertex ij  now the i and j correspond to a vertex  i = 0 1 2 (vertex 1 )  3 4 5 (vertex 2 ) 6 7 8 (vertex 3)
                if( Hessian_block(i,j) > 1e-12 || Hessian_block(i,j) < -1e-12) tripletList.push_back(T(indices[i/3]*3+i%3,indices[j/3]*3+j%3,Hessian_block(i,j)) );  
                if( Hessian_block_2(i,j) > 1e-12 || Hessian_block_2(i,j) < -1e-12) tripletList.push_back(T(indices[i/3]*3+i%3,indices[j/3]*3+j%3,Hessian_block_2(i,j)) );
            }
        } 
        


    }
    Hessian.setFromTriplets(tripletList.begin(),tripletList.end());
    // That gives me the hessian of the surface tension :O.


    return KA*Hessian;
}












void E_Handler::Calculate_energies(double* E){

    // So this is the calculation of the energies
    // std::cout<<"reassinginn \n";
    *E = 0;

    // std::cout<<"The value of energy is " << *E <<"  \n";
    Energy_values.resize(Energies.size());

    int bead_count = 0;
    // std::cout<<" iterating\n";
    for(size_t i = 0; i < Energies.size(); i++){
        // std::cout<<" The energy is "<< Energies[i]<<"\n";
        if(Energies[i] == "Volume_constraint"){
            if(Energy_constants[i][0]>1e-5) Energy_values[i] = E_Volume_constraint(Energy_constants[i]);
            else Energy_values[i] = 0.0;
            *E += Energy_values[i];
            continue;
        
        }

        if(Energies[i]=="Area_constraint"){
            if(Energy_constants[i][0]>1e-5) Energy_values[i] = E_Area_constraint(Energy_constants[i]);
            else Energy_values[i] = 0.0;

            // std::cout<<"THe target area is " << Energy_constants[i][1] <<" \n";
            *E += Energy_values[i];
            continue;
        }

        if(Energies[i]=="Surface_tension" || Energies[i] == "H1_Surface_tension" || Energies[i] == "H2_Surface_tension"){
            if(Energy_constants[i][0]<1e-5){
                Energy_values[i] = 0.0;
                continue;
            }
            Energy_values[i] = E_SurfaceTension(Energy_constants[i]);
            *E += Energy_values[i];
            // std::cout<<"The energy value is " << Energy_values[i]<<" \n";
            // std::cout<<"The value of E is" << *E <<" \n";
            continue;
        }

        if(Energies[i]=="Bending" || Energies[i] == "H1_Bending" || Energies[i]=="H2_Bending" ){
            if(Energy_constants[i][0]<1e-5){
                Energy_values[i] = 0.0;
                continue;
            }

            Energy_values[i] = E_Bending(Energy_constants[i]);
            // std::cout<<"succesfully bent \n";
            *E += Energy_values[i];
            
            continue;
        }
        if(Energies[i]=="Bead" || Energies[i]=="H1_Bead" || Energies[i]=="H2_Bead")
        {
            // Energy_values[i] = Beads[bead_count]->Energy();
            Energy_values[i] = Beads[bead_count]->Bead_I->Tot_Energy();
            // double E_bead = Beads[bead_count]->Energy();
            // if(Energy_values[i]-E_bead > 1e1) std::cout<<"The oldE is" << E_bead <<" and the new one is" << Energy_values[i] <<" \n";
            *E += Energy_values[i];
            bead_count++;
            // std::cout<<"The energy value is " << Energy_values[i]<<" \n";
            // std::cout<<"The value of E is" << *E <<" \n";
            continue;
        }

        if(Energies[i] =="Laplace"){
            Energy_values[i] = E_Laplace(Energy_constants[i]);
            *E += Energy_values[i];
            // std::cout<<"The energy value is " << Energy_values[i]<<" \n";
            // std::cout<<"The value of E is" << *E <<" \n";
            continue;
        }

        if(Energies[i]=="Edge_reg"){
            Energy_values[i] = E_Edge_reg(Energy_constants[i]);
            *E += Energy_values[i];
            // std::cout<<"The energy value is " << Energy_values[i]<<" \n";
            // std::cout<<"The value of E is" << *E <<" \n";
            continue;
        }
        if(Energies[i]=="Face_reg"){
            if(Face_reference.size() == 0){
               update_face_reference();
            }
            Energy_values[i] = E_Face_reg(Energy_constants[i]);
            *E += Energy_values[i];
            // std::cout<<"The energy value is " << Energy_values[i]<<" \n";
            // std::cout<<"The value of E is" << *E <<" \n";
            continue;
            
            }

    }



 
    return;
}

void E_Handler::Calculate_Lag_norm(double* Norm){

    // Now  i need to get this norm
    // std::cout<<""
    Calculate_gradient();
    Calculate_Jacobian();



    Vector3 Force;
    double val = 0;
    Eigen::VectorXd LambdaJ = Jacobian_constraints.transpose()*Lagrange_mult;
    
    double Lambdaval = 0.0;    
    for(size_t vi = 0; vi < mesh->nVertices(); vi++){
        Force = Current_grad[vi];

        val = LambdaJ(3*vi)+Force.x;
        *Norm += val*val;
        
        val = LambdaJ(3*vi+1)+Force.y;
        *Norm += val*val;
    
        val = LambdaJ(3*vi+2)+Force.z;
        *Norm += val*val;

    }
    // std::cout<<"The forces contributed " << *Norm << "\n";
    // std::cout<<"The Lagrange multipliers contributed " << Lambdaval << "\n";
    for(size_t bi = 0; bi < Beads.size(); bi ++){
        val = Beads[bi]->Total_force.norm2();
        *Norm +=val;

    }
  
    // I need to add the beads 
    for(int i = 0 ; i < N_constraints; i++){
        if(Constraints[i]=="Volume"){
        val = -1*(geometry->totalVolume()- Trgt_vol);
        *Norm += val*val;
        }
        if(Constraints[i]=="Area"){
            val = -1*(geometry->totalArea() - Trgt_area);
            *Norm +=val*val;
        }
    }
 



}


void E_Handler::Calculate_gradient(){
    // This function will calculate the gradient of the energy

    // std::cout<<"1 \n";
    Previous_grad = Current_grad;
    // std::cout<<"2 \n";
    Current_grad = VertexData<Vector3>(*mesh, Vector3{0.0, 0.0, 0.0});
    // std::cout<<"3 \n";
    VertexData<Vector3> Force_temp(*mesh, Vector3{0.0, 0.0, 0.0});

    int bead_count = 0;
    double grad_norm = 0;
    // std::cout<<"4 \n";
    // std::cout<<"THe size of Energies is " << Energies.size() << "\n";
    // std::cout<<"4 1 \n";
    for(size_t i = 0; i < Energies.size(); i++){
        // std::cout<<"Energy is " << Energies[i]<<" \n";
// 
        if(Energies[i] == "Volume_constraint"){
            // std::cout<<"DOing volume constraint\n";
            if(Energy_constants[i][0]<1e-5){
                if(Gradient_norms.size() == i){
                Gradient_norms.push_back(0.0);
            }
            else{
                Gradient_norms[i] = 0.0;
            }
                continue;
            }

            Force_temp = F_Volume_constraint(Energy_constants[i]);
            grad_norm = 0.0;
            for(Vertex v : mesh->vertices()){
                grad_norm += Force_temp[v].norm2();
            }
            if(Gradient_norms.size() == i){
                Gradient_norms.push_back(grad_norm);
            
            }
            else{
                Gradient_norms[i] = grad_norm;
            }

            Current_grad+=Force_temp;


            continue;
        }

        if(Energies[i]== "Area_constraint"){
            if(Energy_constants[i][0]<1e-5){
                if(Energy_constants[i][0]<1e-5){
                if(Gradient_norms.size() == i){
                Gradient_norms.push_back(0.0);
            }
            else{
                Gradient_norms[i] = 0.0;
            }
                continue;
            }
            continue;   
            
            }
            Force_temp = F_Area_constraint(Energy_constants[i]);
            grad_norm = 0.0;
            for(Vertex v : mesh->vertices()){
                grad_norm += Force_temp[v].norm2();
            }
            if(Gradient_norms.size() == i){
                Gradient_norms.push_back(grad_norm);
            
            }
            else{
                Gradient_norms[i] = grad_norm;
            }

            Current_grad+=Force_temp;


            continue;
        }

        if(Energies[i]=="Surface_tension"){
            // std::cout<<"Surface tension \n";
            Force_temp = F_SurfaceTension_2(Energy_constants[i]);
            
            grad_norm = 0.0;
            for(Vertex v : mesh->vertices()){
                grad_norm += Force_temp[v].norm2();
            }
            if(Gradient_norms.size() == i){
                Gradient_norms.push_back(grad_norm);
            
            }
            else{
                Gradient_norms[i] = grad_norm;
            }

            Current_grad+=Force_temp;

            continue;

           
        }

        if(Energies[i]=="Bending"){
            // std::cout<<"Bending\n";
            Force_temp = F_Bending(Energy_constants[i]);
            grad_norm = 0.0;

            for(Vertex v :mesh->vertices()){
                grad_norm += Force_temp[v].norm2();
            }
            if(Gradient_norms.size() == i){
                Gradient_norms.push_back(grad_norm);
            
            }
            else{
                Gradient_norms[i] = grad_norm;
            }

            Current_grad+=Force_temp;


            continue;
        }
        // I need to add the beads 

        if(Energies[i]=="Bead"){
            // std::cout<<"Bead E \n";
            // std::cout<<"The number of beads is " << Beads.size() << "\n";
            // std::cout<<"This is bead" << Beads[bead_count]->Bead_id << "\n";
            Force_temp = Beads[bead_count]->Bead_I->Gradient();
            // Force_tem2 = Beads[bead_count]->Gradient();
            // std::cout<<"The bead feels a force of " << Beads[bead_count]->Total_force.norm2() << "\n";
            // std::cout<<"Interactions \n";
            // Beads[bead_count]->Bead_interactions();
            // std::cout<<"? \n";
            grad_norm = 0;
            for(size_t j = 0; j < mesh->nVertices(); j++){
             grad_norm+= Force_temp[j].norm2();
            }
            // std::cout<<"not this \n";
            if(Gradient_norms.size() == i){
                Gradient_norms.push_back(grad_norm);
            
            }
            else{
                Gradient_norms[i] = grad_norm;
            }
            // std::cout<<"? ?\n";
            
            Current_grad+=Force_temp;
            bead_count +=1;
            // if(bead_count==2) std::cout<<"\t done for the iteration\n";

            continue;
        }

        if(Energies[i]=="Laplace"){

            // std::cout<<"Bending\n";
            Force_temp = F_Laplace(Energy_constants[i]);
            grad_norm = 0.0;

            for(Vertex v :mesh->vertices()){
                grad_norm += Force_temp[v].norm2();
            }
            if(Gradient_norms.size() == i){
                Gradient_norms.push_back(grad_norm);
            
            }
            else{
                Gradient_norms[i] = grad_norm;
            }
            Current_grad+=Force_temp;
            continue;
        }

        if(Energies[i]=="Edge_reg"){
            Force_temp = F_Edge_reg(Energy_constants[i]);
            grad_norm = 0.0;
            for(Vertex v : mesh->vertices()){
                grad_norm += Force_temp[v].norm2();
            }
            if(Gradient_norms.size() == i){
                Gradient_norms.push_back(grad_norm);
            
            }
            else{
                Gradient_norms[i] = grad_norm;
            }
            Current_grad+=Force_temp;
            continue;
        }

        if(Energies[i]=="Face_reg"){
            if(Face_reference.size()==0) update_face_reference();

            Force_temp = F_Face_reg(Energy_constants[i]);
            grad_norm = 0.0;
            for(Vertex v: mesh->vertices()){
                grad_norm += Force_temp[v].norm2();

            }
            if(Gradient_norms.size() == i){
                Gradient_norms.push_back(grad_norm);
            }
            else{
                Gradient_norms[i] = grad_norm;
            }
            Current_grad+=Force_temp;
            continue;
        }

        }
    // std::cout<<"5 \n";
    grad_norm = 0.0;
    for(size_t i = 0; i < mesh->nVertices(); i++){
        grad_norm+=Current_grad[i].norm2();
    }
    // std::cout<<"6 \n";
    if(Gradient_norms.size() == Energies.size()){
            Gradient_norms.push_back(grad_norm);
        }
        else{
            Gradient_norms[Energies.size()] = grad_norm;
            }

    return;
}

void E_Handler::Calculate_Jacobian(){

    int N_verts = mesh->nVertices();
    int N_constraints = 0;
    int N_beads = Beads.size();

    // std::vector<double> Constraint_values(0);
    
    VertexData<Vector3> grad_sur(*mesh);
    VertexData<Vector3> grad_vol(*mesh);
    std::vector<double> Energy_constants_val;
    for(std::string constraint : Constraints){
        if(constraint == "Volume"){
            // I need to find which are the energy constants
            grad_vol = F_Volume(std::vector<double>{-1.0});
            N_constraints += 1;
            // std::cout<<"Volume constraint added \n";
            // Constraint_values.push_back(geometry->totalVolume());

        }
        if(constraint == "Area"){
            grad_sur = F_SurfaceTension_2(std::vector<double>{-1.0});
            N_constraints += 1;
            // Constraint_values.push_back(geometry->totalArea());
        }
        if(constraint == "CM"){
            N_constraints += 3;        
        }
        if(constraint == "CMx"){
            N_constraints += 1;        
            // Constraint_values.push_back(0.0);
        }

        if(constraint == "CMy"){
            N_constraints += 1;        
            // Constraint_values.push_back(0.0);
        }
        if(constraint == "CMz"){
            N_constraints += 1;        
            // Constraint_values.push_back(0.0);
        }
        if(constraint == "Rx"){
            N_constraints += 1;        
            // Constraint_values.push_back(0.0);
        }

        if(constraint == "Ry"){
            N_constraints += 1;        
            // Constraint_values.push_back(0.0);
        }
        if(constraint == "Rz"){
            N_constraints += 1;        
            // Constraint_values.push_back(0.0);
        }
        

        
    }

    Jacobian_constraints.resize( N_constraints,3*(N_verts+N_beads));
    // std::cout<<"The size of the Jacobian is " << Jacobian_constraints.rows() << " " << Jacobian_constraints.cols() << "\n";

    // Now we need to fill the Jacobian matrix with the constraints

    Eigen::VectorXd Column = Eigen::VectorXd::Zero((N_beads+N_verts)*3);
    Vector3 Pos;
    for( size_t i =0; i < Constraints.size(); i++){
        Column = Eigen::VectorXd::Zero( (N_verts+N_beads)*3);
        if(Constraints[i] == "Volume"){
            // std::cout<<"Volume constraint \n";
            for(size_t vi = 0; vi < mesh->nVertices() ; vi++  ){
                Column[3*vi] = grad_vol[vi].x;
                Column[3*vi+1] = grad_vol[vi].y;
                Column[3*vi+2] = grad_vol[vi].z;
            }
            Jacobian_constraints.row(i) = Column;
            
            continue;
        }

        if(Constraints[i]== "Area"){

            for(size_t vi = 0; vi < mesh->nVertices(); vi++){
                Column[3*vi] = grad_sur[vi].x;
                Column[3*vi+1] = grad_sur[vi].y;
                Column[3*vi+2] = grad_sur[vi].z;
            }
            Jacobian_constraints.row(i) = Column;
        }
        if(Constraints[i]=="CMx"){
            // std::cout<<"CMz for the win \n";
            for(size_t i = 0; i < mesh->nVertices()+N_beads; i++){
                Column[3*i] = 1;
            }
            Jacobian_constraints.row(i) = Column;
            continue;
        }
        if(Constraints[i]=="CMy"){
            // std::cout<<"CMz for the win \n";
            for(size_t i = 0; i < mesh->nVertices()+N_beads; i++){
                Column[3*i+1] = 1;
            }
            
            Jacobian_constraints.row(i) = Column;
            continue;
        }
        if(Constraints[i]=="CMz"){
            // std::cout<<"CMz for the win \n";
            for(size_t i = 0; i < mesh->nVertices()+N_beads; i++){
                Column[3*i+2] = 1;
            }
            
            Jacobian_constraints.row(i) = Column;
            continue;
        }
        if(Constraints[i]=="Rx"){
            // std::cout<<"RX for the win \n";
            for(size_t vi = 0; vi < N_verts; vi++){
                Pos = geometry->inputVertexPositions[vi];
                Column[3*vi+1] = Pos.z;
                Column[3*vi+2] = -Pos.y;
            }
            for(size_t bi = 0; bi < N_beads; bi++){
                Pos = Beads[bi]->Pos;
                Column[3*(N_verts+bi)+1] = Pos.z;
                Column[3*(N_verts+bi)+2] = -Pos.y;
            }
            Jacobian_constraints.row(i) = Column;
            continue;
        }
        if(Constraints[i]=="Ry"){
            // std::cout<<"RY for the win \n";
            for(size_t vi = 0; vi < N_verts; vi++){
                Pos = geometry->inputVertexPositions[vi];
                Column[3*vi] = -Pos.z;
                Column[3*vi+2] = Pos.x;
            }
            for(size_t bi = 0; bi < N_beads; bi++){
                Pos = Beads[bi]->Pos;
                Column[3*(N_verts+bi)] = -Pos.z;
                Column[3*(N_verts+bi)+2] = Pos.x;
            }
            Jacobian_constraints.row(i) = Column;
            continue;
        }
        if(Constraints[i]=="Rz"){
            // std::cout<<"RZ for the win \n";
            for(size_t vi = 0; vi < mesh->nVertices(); vi++){
                Pos = geometry->inputVertexPositions[vi];
                Column[3*vi] = Pos.y;
                Column[3*vi+1] = -Pos.x;
            }
            for(size_t bi = 0; bi < N_beads; bi++){
                Pos = Beads[bi]->Pos;
                Column[3*(N_verts+bi)] = Pos.y;
                Column[3*(N_verts+bi)+1] = -Pos.x;
            }
            Jacobian_constraints.row(i) = Column;
            continue;
        }       

    }



}


SparseMatrix<double> E_Handler::Calculate_Hessian(){

    // Ok so the Hessian is the normal Hessian + lagrange multipliers * Hessian of the gradients
    // For beads you need to include the beads at the size of the matrix of the hessian
    // I recommend adding them at the bottom. 

    int N_verts = mesh->nVertices();
    int N_beads = Beads.size();
    // std::cout<<"THe number of beads is " << N_beads << "\n";
    // std::cout<<"Calculating energy\n";
    std::vector<double> Energy_constants_val;
    SparseMatrix<double> Hessian(3*(N_verts+N_beads),3*(N_verts+N_beads));
    std::string constraint;
    


    int bead_counter= 0;
    // std::cout<<"THe size of Energies is "<< Energies.size() << "\n";
    for(size_t i = 0; i < Energies.size(); i++){
        // std::cout<<"Energy is " << Energies[i] <<" \n";
        if(Energies[i] == "Bending"){
            Hessian += H_Bending(Energy_constants[i]);
        }
        if(Energies[i] == "Surface_tension"){
            Hessian += H_SurfaceTension(Energy_constants[i]);
        }
        if(Energies[i] =="Laplace"){
            Hessian += H_Laplace(Energy_constants[i]);
        }
        if(Energies[i] == "Bead"){
            // std::cout<<"Doing bead energy\n";
            // std::cout<<"The og hessian ahs size" << Hessian.rows() << " " << Hessian.cols() << "\n";
            // std::cout<<"The energy constant is "<< Beads[bead_counter]->Bead_I->Energy_constants[0] << " \n";
            // std::cout<<"FUNCTION NOT AVAILABLE BUT SHOULD LOOK LIKE this\n";
            Hessian += Beads[bead_counter]->Bead_I->Hessian();
            // std::cout<<"Done with bead hessian\n";
            bead_counter +=1;
        }
        if(Energies[i] =="Edge_reg"){
            // std::cout<<"There should be one energy constant\n";
            // std::cout<<"The energy constants are " << Energy_constants[i][0]<<" \n";
            Hessian += H_Edge_reg(Energy_constants[i]);
        }
        if(Energies[i] =="Face_reg"){
            Hessian += H_Face_reg(Energy_constants[i]);
        }

    }
    // std::cout<<"DOne calculating hessian\n";

    // Eigen::MatrixXd Hessian_E = Hessian.toDense();
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Hessian_E);

    // Eigen::VectorXd Eigenvalues = solver.eigenvalues();

    // double min_eig = Eigenvalues.minCoeff();
    // std::cout<<"The minimum eigenvalue is " << min_eig << "\n";
    
    // SparseMatrix<double> Id_reg = 
    // std::cout<<"Doing the volume const\n";
    // Eigen::SparseMatrix::Ide
    for(size_t i = 0; i < Constraints.size() ; i++){
        constraint = Constraints[i];
        // std::cout<<"The constraint is " << constraint<<"\n";
        if(constraint == "Volume"){
            // I need to find which are the energy constants
            // std::cout<<"Doing the volume constraint\n";       
            Hessian += -1.0*Lagrange_mult(i)*H_Volume(std::vector<double>{1.0});
            // std::cout<<"DONE WITH VOLUME CONSTRAINT\n";
        }
        if(constraint == "Area"){
            Hessian += -1.0*Lagrange_mult(i)*H_SurfaceTension(std::vector<double>{1.0});
        }

        
    }

    // std::cout<<"DOne with the Hessian\n";
    // I want this t
    

    return Hessian;
}


SparseMatrix<double> E_Handler::Calculate_Hessian_E(){

    // Ok so the Hessian is the normal Hessian + lagrange multipliers * Hessian of the gradients
    // For beads you need to include the beads at the size of the matrix of the hessian
    // I recommend adding them at the bottom. 
    // std::cout<<"Creating the Hessian E\n";
    int N_verts = mesh->nVertices();
    int N_beads = Beads.size();
    // std::cout<<"THe number of beads is " << N_beads << "\n";
    // std::cout<<"Calculating energy\n";
    std::vector<double> Energy_constants_val;
    SparseMatrix<double> Hessian(3*(N_verts+N_beads),3*(N_verts+N_beads));
    std::string constraint;
    


    int bead_counter= 0;
    // std::cout<<"THe size of Energies is "<< Energies.size() << "\n";
    for(size_t i = 0; i < Energies.size(); i++){
        // std::cout<<"Energy is " << Energies[i] <<" \n";
        if(Energies[i] == "Bending"){
            Hessian += H_Bending(Energy_constants[i]);
        }
        if(Energies[i] == "Surface_tension"){
            Hessian += H_SurfaceTension(Energy_constants[i]);
        }
        if(Energies[i] == "Bead"){
            // std::cout<<"Doing bead energy\n";
            // std::cout<<"The og hessian ahs size" << Hessian.rows() << " " << Hessian.cols() << "\n";
            // std::cout<<"The energy constant is "<< Beads[bead_counter]->Bead_I->Energy_constants[0] << " \n";
            // std::cout<<"FUNCTION NOT AVAILABLE BUT SHOULD LOOK LIKE this\n";
            Hessian += Beads[bead_counter]->Bead_I->Hessian_IP();
            // std::cout<<"Done with bead hessian\n";
            bead_counter +=1;
        }
        if(Energies[i] =="Laplace"){
            Hessian += H_Laplace(Energy_constants[i]);
        }

        if(Energies[i] =="Edge_reg"){
            // std::cout<<"There should be one energy constant\n";
            // std::cout<<"The energy constants are " << Energy_constants[i][0]<<" \n";
            Hessian += H_Edge_reg(Energy_constants[i]);
        }
        if(Energies[i] =="Face_reg"){
            Hessian += H_Face_reg(Energy_constants[i]);
        }

    }

    

    return Hessian;
}


SparseMatrix<double> E_Handler::Calculate_Hessian_Constraints(){

    // Ok so the Hessian is the normal Hessian + lagrange multipliers * Hessian of the gradients
    // For beads you need to include the beads at the size of the matrix of the hessian
    // I recommend adding them at the bottom. 

    int N_verts = mesh->nVertices();
    int N_beads = Beads.size();
    // std::cout<<"THe number of beads is " << N_beads << "\n";
    // std::cout<<"Calculating energy\n";
    std::vector<double> Energy_constants_val;
    SparseMatrix<double> Hessian(3*(N_verts+N_beads),3*(N_verts+N_beads));
    std::string constraint;
    

    for(size_t i = 0; i < Constraints.size() ; i++){
        constraint = Constraints[i];
        // std::cout<<"The constraint is " << constraint<<"\n";
        if(constraint == "Volume"){
            // I need to find which are the energy constants
            // std::cout<<"Doing the volume constraint\n";       
            Hessian += Lagrange_mult(i)*H_Volume(std::vector<double>{1.0});
            // std::cout<<"DONE WITH VOLUME CONSTRAINT\n";
        }
        if(constraint == "Area"){
            Hessian += Lagrange_mult(i)*H_SurfaceTension(std::vector<double>{1.0});
        }

        
    }

    // std::cout<<"DOne with the Hessian\n";
    // I want this t
    

    return Hessian;
}


void E_Handler::Do_nothing(){
    // This function does nothing, it is just a placeholder
    return;
}




//   for(size_t i = 0; i < Energies.size(); i++){
    
    
    // if(Energies[i]=="H1_Bead"){
    //   if(!H1_grad){
    //     construction_start = chrono::steady_clock::now();
    //     H1_grad = true;
        
    //     H1_mat = H1_operator(Energy_constants[i][1],Energy_constants[i][2],Energy_constants[i][3]);
    //     construction_end = chrono::steady_clock::now();
    //     time_construct += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();

    //     construction_start = chrono::steady_clock::now();
    //     solverH1.compute(H1_mat);
    //     construction_end = chrono::steady_clock::now();
    //     time_compute += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();
      
    //   }

    //   Energy_vals[i] = Beads[bead_count]->Energy(); 
    //   // std::cout<<"Bead count is" << bead_count <<"\n";
    //   Force_temp = Beads[bead_count]->Gradient();
    //   Beads[bead_count]->Bead_interactions();

    //   solve_start = chrono::steady_clock::now();
    //   Vector<double> RHS = Vector<double>(N_vert*3+N_constraints);
    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     RHS[3*index] = Force_temp[index].x;
    //     RHS[3*index+1] = Force_temp[index].y;
    //     RHS[3*index+2] = Force_temp[index].z;
    //   }
    //   for(int index = 0; index < N_constraints; index++){
    //     RHS[3*N_vert+index] = 0.0;
    //   }
      
    //   Vector<double> Solution = solverH1.solve(RHS);
    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     Force_temp[index] = Vector3({Solution.coeff(3*index),Solution.coeff(3*index+1),Solution.coeff(3*index+2)});
    //     // std::cout<< Force_temp[index].norm2() << " ";
    //   }
    //   solve_end = chrono::steady_clock::now();
    //   time_solve += std::chrono::duration_cast<std::chrono::milliseconds>(solve_end - solve_start).count();


    //   grad_value = 0;
    //   for(size_t j = 0; j < mesh->nVertices(); j++){
    //     grad_value+= Force_temp[j].norm2();
    //   }
    //   Gradient_norms.push_back(grad_value);
    //   Force+=Force_temp;
    //   bead_count +=1;
    //   continue;

    // }
    // if(Energies[i]=="H2_Bead"){
      
    //   if(!H2_grad){
    //     construction_start = chrono::steady_clock::now();
    //     H2_grad = true;
    //     N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];

    //     H2_mat = H2_operator(Energy_constants[i][1],Energy_constants[i][2],Energy_constants[i][3]);
    //     construction_end = chrono::steady_clock::now();
    //     time_construct += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();

    //     construction_start = chrono::steady_clock::now();
    //     solverH2.compute(H2_mat);
    //     construction_end = chrono::steady_clock::now();
    //     time_compute += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();
      
    //   }

    //   Energy_vals[i] = Beads[bead_count]->Energy(); 
    //   // std::cout<<"Bead count is" << bead_count <<"\n";
    //   Force_temp = Beads[bead_count]->Gradient();
    //   Beads[bead_count]->Bead_interactions();

    //   solve_start = chrono::steady_clock::now();
    //   Vector<double> RHS = Vector<double>(N_vert*3+N_constraints);
    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     RHS[3*index] = Force_temp[index].x;
    //     RHS[3*index+1] = Force_temp[index].y;
    //     RHS[3*index+2] = Force_temp[index].z;
    //   }
    //   for(int index = 0; index < N_constraints; index++){
    //     RHS[3*N_vert+index] = 0.0;
    //   }
      
    //   Vector<double> Solution = solverH2.solve(RHS);
    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     Force_temp[index] = Vector3({Solution.coeff(3*index),Solution.coeff(3*index+1),Solution.coeff(3*index+2)});
    //     // std::cout<< Force_temp[index].norm2() << " ";
    //   }
    //   solve_end = chrono::steady_clock::now();
    //   time_solve += std::chrono::duration_cast<std::chrono::milliseconds>(solve_end - solve_start).count();


    //   grad_value = 0;
    //   for(size_t j = 0; j < mesh->nVertices(); j++){
    //     grad_value+= Force_temp[j].norm2();
    //   }
    //   Gradient_norms.push_back(grad_value);
    //   Force+=Force_temp;
    //   bead_count +=1;
    //   // std::cout<<"Here\n";
    //   continue;


    // }

    // if(Energies[i]=="H2_Bending"){
    //   // std::cout<<"THe beding energy is hereee\n";
    //   if(!H2_grad){
    //     construction_start = chrono::steady_clock::now();
    //     H2_grad = true;
    //     N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];

    //     H2_mat = H2_operator(Energy_constants[i][1],Energy_constants[i][2],Energy_constants[i][3]);
    //     construction_end = chrono::steady_clock::now();
    //     time_construct += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();

    //     construction_start = chrono::steady_clock::now();
    //     solverH2.compute(H2_mat);
    //     construction_end = chrono::steady_clock::now();
    //     time_compute += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();
    //     // std::cout<<"Solving possible\n";
    //   }

    //   double KB = Energy_constants[i][0];
    //   Energy_vals[i] = E_Bending(0.0,KB);
    //   // Now i need to compute the original gradient
    //   Force_temp = KB*Bending(0.0);

    //   // I would say the solve considers the construction?

    //   N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];
    //   solve_start = chrono::steady_clock::now();
    //   Vector<double> RHS = Vector<double>(N_vert*3+N_constraints);
    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     RHS[3*index] = Force_temp[index].x;
    //     RHS[3*index+1] = Force_temp[index].y;
    //     RHS[3*index+2] = Force_temp[index].z;
    //     // std::cout<<"Force temp has magnitude" << Force_temp[index].norm2() <<"and barycentric " << Barycentric_area[index]<<"\n";
    //   }
    //   // std::cout<<"THe number of constraints is " << N_constraints<<" duh\n";
    //   for(int index = 0; index < N_constraints; index++){
    //     RHS[3*N_vert+index] = 0.0;
    //   }
    //   // std::cout<<"The RHS is readyy\n";
    //   // std::cout<<"Some components are " << RHS[0] <<" "<< RHS[1] <<" "<< RHS[2] <<" \n";
    //   // The RHS is readyy   
    //   // std::cout<<"Solving RHS\n";
      
    //   Vector<double> Solution = solverH2.solve(RHS);
    //   // std::cout<<"Solved\n";
    //   // std::cout<<"The solution is ready\n";
    //   // std::cout<<"Some components of the solution are " << Solution[0] <<" "<< Solution[1] <<" "<< Solution[2] <<" \n";
    //   // std::cout<<"The forces are ";
    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     Force_temp[index] = Vector3({Solution.coeff(3*index),Solution.coeff(3*index+1),Solution.coeff(3*index+2)});
    //     // std::cout<< Force_temp[index].norm2() << " ";
    //   }
    //   solve_end = chrono::steady_clock::now();
    //   time_solve += std::chrono::duration_cast<std::chrono::milliseconds>(solve_end - solve_start).count();
      
    //   // SO i have calculated the new force
    //   grad_value = 0;
    //   for(size_t j = 0; j < mesh->nVertices(); j++){
    //     grad_value+= Force_temp[j].norm2();
    //   }
    //   // std::cout<<"THe grad value is " << grad_value <<" \n";
    //   Gradient_norms.push_back(grad_value);
    //   Force+=Force_temp;

    //   continue;

    // }
    // if(Energies[i]=="H1_Bending"){
    //   // std::cout<<"THe beding energy is hereee\n";
    //   if(!H1_grad){
    //     construction_start = chrono::steady_clock::now();
    //     H1_grad = true;
    //     N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];

    //     // I need to build the sobolev operator
    //     // H2_mat = SparseMatrix<double>(N_vert*3+N_constraints,N_vert*3+N_constraints);
    //     H1_mat = H1_operator(Energy_constants[i][1],Energy_constants[i][2],Energy_constants[i][3]);
    //     construction_end = chrono::steady_clock::now();
    //     time_construct += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();

    //     // std::cout<<"THe solver is being solved\n"; 
    //     // std::cout<<"THe matrix has some nonzero components" << H2_mat.nonZeros() <<"\n";
    //     construction_start = chrono::steady_clock::now();
    //     solverH1.compute(H1_mat);
    //     construction_end = chrono::steady_clock::now();
    //     time_compute += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();
    //     // std::cout<<"Solving possible\n";
    //   }

    //   double KB = Energy_constants[i][0];
    //   Energy_vals[i] = E_Bending(0.0,KB);
    //   // Now i need to compute the original gradient
    //   Force_temp = KB*Bending(0.0);

    //   // double temp_grad = 0;
    //   // for(size_t j = 0; j < mesh->nVertices(); j++){
    //   //   temp_grad += Force_temp[j].norm2();
    //   // }
    //   // std::cout<<"Temp grad is " << temp_grad <<" \n";
    //   //Now i need to add the constraint

    //   // I would say the solve considers the construction?

    //   N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];
    //   solve_start = chrono::steady_clock::now();
    //   Vector<double> RHS = Vector<double>(N_vert*3+N_constraints);
    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     RHS[3*index] = Force_temp[index].x;
    //     RHS[3*index+1] = Force_temp[index].y;
    //     RHS[3*index+2] = Force_temp[index].z;
    //     // std::cout<<"Force temp has magnitude" << Force_temp[index].norm2() <<"and barycentric " << Barycentric_area[index]<<"\n";
    //   }
    //   // std::cout<<"THe number of constraints is " << N_constraints<<" duh\n";
    //   for(int index = 0; index < N_constraints; index++){
    //     RHS[3*N_vert+index] = 0.0;
    //   }
    //   // std::cout<<"The RHS is readyy\n";
    //   // std::cout<<"Some components are " << RHS[0] <<" "<< RHS[1] <<" "<< RHS[2] <<" \n";
    //   // The RHS is readyy   
    //   // std::cout<<"Solving RHS\n";
      
    //   Vector<double> Solution = solverH1.solve(RHS);
    //   // std::cout<<"Solved\n";
    //   // std::cout<<"The solution is ready\n";
    //   // std::cout<<"Some components of the solution are " << Solution[0] <<" "<< Solution[1] <<" "<< Solution[2] <<" \n";
    //   // std::cout<<"The forces are ";
    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     Force_temp[index] = Vector3({Solution.coeff(3*index),Solution.coeff(3*index+1),Solution.coeff(3*index+2)});
    //     // std::cout<< Force_temp[index].norm2() << " ";
    //   }
    //   solve_end = chrono::steady_clock::now();
    //   time_solve += std::chrono::duration_cast<std::chrono::milliseconds>(solve_end - solve_start).count();
      
    //   // SO i have calculated the new force
    //   grad_value = 0;
    //   for(size_t j = 0; j < mesh->nVertices(); j++){
    //     grad_value+= Force_temp[j].norm2();
    //   }
    //   // std::cout<<"THe grad value is " << grad_value <<" \n";
    //   Gradient_norms.push_back(grad_value);
    //   Force+=Force_temp;

    //   continue;

    // }
    // if(Energies[i]=="H2_Surface_tension"){
    //   if(!H2_grad){
    //     construction_start = chrono::steady_clock::now();
        
    //     H2_grad = true;
    //     N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];

    //     // I need to build the sobolev operator
    //     // H2_mat = SparseMatrix<double>(N_vert*3+N_constraints,N_vert*3+N_constraints);
    //     H2_mat = H2_operator(Energy_constants[i][1],Energy_constants[i][2],Energy_constants[i][3]);

    //     construction_end = chrono::steady_clock::now();
    //     time_construct += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();

    //     construction_start = chrono::steady_clock::now();
    //     solverH2.compute(H2_mat);
    //     construction_end = chrono::steady_clock::now();
    //     time_compute += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();
        
    //     // std::cout<<"Solving possible\n";
    //   }
    //   double sigma = Energy_constants[i][0];
    //   A = geometry->totalArea();
    //   Energy_vals[i] = A*sigma;
    //   Force_temp = sigma*SurfaceGrad();

    //   // I need to multiply by the mass dont I?


    //   N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];

    //   Vector<double> RHS = Vector<double>(N_vert*3+N_constraints);
    //   for (size_t i = 0; i < mesh->nVertices(); i++)
    //   {
    //     RHS[3*i] = Force_temp[i].x;
    //     RHS[3*i+1] = Force_temp[i].y;
    //     RHS[3*i+2] = Force_temp[i].z;
    //   }
    //   for(int index = 0; index < N_constraints; index++){
    //     RHS[3*N_vert+index] = 0.0;
    //   }

    //   solve_start = chrono::steady_clock::now();
    //   Vector<double> Solution = solverH2.solve(RHS);
    //   solve_end = chrono::steady_clock::now();
    //   time_solve += std::chrono::duration_cast<std::chrono::milliseconds>(solve_end - solve_start).count();

    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     Force_temp[index] = Vector3({Solution.coeff(3*index),Solution.coeff(3*index+1),Solution.coeff(3*index+2)});
    //   }
      
      
      
    //   grad_value  = 0.0;
    //   for(size_t j = 0; j < mesh->nVertices(); j++){
    //     grad_value+= Force_temp[j].norm2();
    //   }
    //   Gradient_norms.push_back(grad_value);
    //   Force+=Force_temp;
    //   continue;

    // }

    // if(Energies[i]=="H1_Surface_tension"){
    //   if(!H1_grad){
    //     construction_start = chrono::steady_clock::now();
        
    //     H1_grad = true;
    //     N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];

    //     // I need to build the sobolev operator
    //     // H2_mat = SparseMatrix<double>(N_vert*3+N_constraints,N_vert*3+N_constraints);
    //     H1_mat = H1_operator(Energy_constants[i][1],Energy_constants[i][2],Energy_constants[i][3]);

    //     construction_end = chrono::steady_clock::now();
    //     time_construct += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();

    //     construction_start = chrono::steady_clock::now();
    //     solverH1.compute(H1_mat);
    //     construction_end = chrono::steady_clock::now();
    //     time_compute += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();
        
    //     // std::cout<<"Solving possible\n";
    //   }
    //   double sigma = Energy_constants[i][0];
    //   A = geometry->totalArea();
    //   Energy_vals[i] = A*sigma;
    //   Force_temp = sigma*SurfaceGrad();

    //   // I need to multiply by the mass dont I?


    //   N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];

    //   Vector<double> RHS = Vector<double>(N_vert*3+N_constraints);
    //   for (size_t i = 0; i < mesh->nVertices(); i++)
    //   {
    //     RHS[3*i] = Barycentric_area[i]*Force_temp[i].x;
    //     RHS[3*i+1] = Barycentric_area[i]*Force_temp[i].y;
    //     RHS[3*i+2] = Barycentric_area[i]*Force_temp[i].z;
    //   }
    //   for(int index = 0; index < N_constraints; index++){
    //     RHS[3*N_vert+index] = 0.0;
    //   }

    //   solve_start = chrono::steady_clock::now();
    //   Vector<double> Solution = solverH1.solve(RHS);
    //   solve_end = chrono::steady_clock::now();
    //   time_solve += std::chrono::duration_cast<std::chrono::milliseconds>(solve_end - solve_start).count();

    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     Force_temp[index] = Vector3({Solution.coeff(3*index),Solution.coeff(3*index+1),Solution.coeff(3*index+2)});
    //   }
      
      
      
    //   grad_value  = 0.0;
    //   for(size_t j = 0; j < mesh->nVertices(); j++){
    //     grad_value+= Force_temp[j].norm2();
    //   }
    //   Gradient_norms.push_back(grad_value);
    //   Force+=Force_temp;
    //   continue;

    // }

//   }