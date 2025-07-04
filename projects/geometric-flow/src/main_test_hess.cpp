// GEOMETRIC FLOW

// #include <stdlib.h>
#include <unistd.h>


#include <omp.h>


#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <string>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/surface/remeshing.h"




#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"

#include <chrono>


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"


#include "Mem-3dg.h"
#include "Beads.h"
#include "Energy_Handler.h"
#include "math.h"


#include "libarcsim/include/cloth.hpp"
#include "libarcsim/include/collision.hpp"
#include "libarcsim/include/mesh.hpp"
#include "libarcsim/include/dynamicremesh.hpp"

#include "io.hpp"
#include "simulation.hpp"
#include "conf.hpp"
#include "log.hpp"


using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace std;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;
// so we can more easily pass these to different classes
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;



// Polyscope visualization handle, to quickly add data to the surface
// polyscope::SurfaceMesh* psMesh;

// Some global variables
float TIMESTEP = -4;
float kappa = 1.0;


float H0 = 1.0;
float V_bar= (4/3)*PI*10; 



float nu;

float c0;


float P0=10000.0;
float KA=1.0000;
float KB=0.000001;
float Kd=0.0;
double TS=0.001;

double Curv_adap=0.1;
double Min_rel_length=0.5;
double trgt_len;
double avg_remeshing;
bool edge_length_adj;
VertexData<Vector3> ORIG_VPOS; // original vertex positions
Vector3 CoM;                   // original center of mass

Mem3DG M3DG;
Bead Bead_1;

std::array<double, 3> BLUE = {0.11, 0.388, 0.89};
// glm::vec<3, float> ORANGE_VEC = {1, 0.65, 0};
std::array<double, 3> ORANGE = {1, 0.65, 0};



void showSelected() {
    // pass
}

// RemeshBoundaryCondition defaultRemeshOptions;
// void redraw() {
//     psMesh->updateVertexPositions(geometry->inputVertexPositions);
//     polyscope::requestRedraw();
// }

// void functionCallback() {

    
//     if (ImGui::Button("Tangential Vertex Smoothing " )) {

//         avg_remeshing = smoothByCircumcenter(*mesh,*geometry);
//         std::cout<< "The average amount the vertex were moved is "<< avg_remeshing << "\n";
//         geometry->normalize(CoM, false);
//         psMesh->remove();
//         psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath("dodecahedra"), geometry->inputVertexPositions,
//                                             mesh->getFaceVertexList(), polyscopePermutations(*mesh));

//         psMesh->setSurfaceColor({0.9607, 0.6627, 0.7215});
//         // psMesh->setSurfaceColor({1.0, 0.45, 0.0}); // orange        
//         polyscope::requestRedraw();
//         redraw();
//     }
//     if (ImGui::Button("Remesh " )) {

//         std::cout<< "The mean edge length is = "<< geometry->meanEdgeLength()<<"\n";
//         std::cout<< "The number of vertices is = "<< mesh->nVertices()<<"\n";
                
//         RemeshOptions Options;
//         Options.targetEdgeLength=0.21;
//         Options.curvatureAdaptation=0.1;
//         Options.maxIterations=10;
//         Options.minRelativeLength=0.2;
//         Options.smoothStyle=RemeshSmoothStyle::Circumcentric;
//         Options.boundaryCondition=RemeshBoundaryCondition::Tangential;
//         ManifoldSurfaceMesh& meshu= *mesh;
//         VertexPositionGeometry& geometryu = *geometry;
//         MutationManager Mutation_manager(meshu,*geometry);
        
//         remesh(*mesh,*geometry,Mutation_manager,Options);

//         geometry->normalize(CoM, false);
//         std::cout<< "The mean edge length after the remeshing is = "<< geometry->meanEdgeLength()<<"\n";
//         std::cout<< "The number of vertices afther the remeshing is  = "<< mesh->nVertices()<<"\n";
//         psMesh->remove();
//         psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath("dodecahedra"), geometry->inputVertexPositions,
//                                             mesh->getFaceVertexList(), polyscopePermutations(*mesh));

//         psMesh->setSurfaceColor({0.9607, 0.6627, 0.7215}); // not orange        
//         polyscope::requestRedraw();
//         redraw();
        
        
        

//     }



//     if (ImGui::Button("Reset")) {
//         geometry->inputVertexPositions = ORIG_VPOS;
//         psMesh->updateVertexPositions(ORIG_VPOS);
//         polyscope::requestRedraw();
//     }
    


//     ImGui::SliderFloat("Timestep", &TIMESTEP, -10.0, 10.0);
//     TS=pow(10,TIMESTEP);
    

// }

// void flipZ() {
//     // Rotate mesh 180 deg about up-axis on startup
//     glm::mat4x4 rot = glm::rotate(glm::mat4x4(1.0f), static_cast<float>(PI), glm::vec3(0, 1, 0));
//     for (Vertex v : mesh->vertices()) {
//         Vector3 vec = geometry->inputVertexPositions[v];
//         glm::vec4 rvec = {vec[0], vec[1], vec[2], 1.0};
//         rvec = rot * rvec;
//         geometry->inputVertexPositions[v] = {rvec[0], rvec[1], rvec[2]};
//     }
//     psMesh->updateVertexPositions(geometry->inputVertexPositions);
// }




void Save_mesh(std::string basic_name, size_t current_t) {
   // Build member variables: mesh, geometry
    Vector3 Pos;
    std::ofstream o(basic_name+"Mem_"+std::to_string(current_t)+".obj");
    o << "#This is a meshfile from a saved state\n" ;

    for(Vertex v : mesh->vertices()) {
    Pos=geometry->inputVertexPositions[v];
    o << "v " << Pos.x <<" "<< Pos.y << " "<< Pos.z <<"\n";

    }


    // I need to save the faces now

    for(Face f : mesh->faces()) {
    o<<"f";
    
    for(Vertex v: f.adjacentVertices()){
        o << " " <<v.getIndex()+1;
    }
    o<<"\n";
    
    }

    return ;
}

void Save_mesh(std::string basic_name,bool arcsim_remeshing, size_t current_t) {
   // Build member variables: mesh, geometry
    Vector3 Pos;
     std::ofstream o;
    if(arcsim_remeshing){
    o.open(basic_name+"Mem_arc_"+std::to_string(current_t)+".obj");
    }
    else{
    o.open(basic_name+"Mem_common_"+std::to_string(current_t)+".obj");
    }
    o << "#This is a meshfile from a saved state\n" ;

    for(Vertex v : mesh->vertices()) {
    Pos=geometry->inputVertexPositions[v];
    o << "v " << Pos.x <<" "<< Pos.y << " "<< Pos.z <<"\n";

    }


    // I need to save the faces now

    for(Face f : mesh->faces()) {
    o<<"f";
    
    for(Vertex v: f.adjacentVertices()){
        o << " " <<v.getIndex()+1;
    }
    o<<"\n";
    
    }

    return ;
}

arcsim::Mesh translate_to_arcsim(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry){

    arcsim::Mesh mesh1;
    // std::cout<< mesh1.verts.size()<<" number of vertices\n";

    // std::cout<<"Adding vertices?\n";
    for (size_t v = 0; v < mesh->nVertices(); v++) {
            // const Vert *vert0 = mesh0.verts[v];
            Vector3 pos_orig= geometry->inputVertexPositions[v];
            arcsim::Vec3 pos;
            pos[0]=pos_orig.x;
            pos[1]=pos_orig.y;
            pos[2]=pos_orig.z;
            arcsim::Vert *vert1 = new arcsim::Vert(pos, 1,0);
            mesh1.add(vert1);
            
    }

    // std::cout<<"Adding nodes?\n";
    for (size_t v = 0; v < mesh->nVertices(); v++) {
            // const Vert *vert0 = mesh0.verts[v];
            Vector3 pos_orig= geometry->inputVertexPositions[v];
            arcsim::Vec3 pos;
            pos[0]=pos_orig.x;
            pos[1]=pos_orig.y;
            pos[2]=pos_orig.z;
            // arcsim::Vert *vert1 = new arcsim::Vert(pos, 0, 1);
            // mesh1.add(vert1);
            arcsim::Node *node1 = new arcsim::Node(pos, pos, pos, 0, false);
            node1->preserve = false;
            node1->temp = false;
            node1->temp2 = false;
            node1->verts.resize(1);
            node1->verts[0] = mesh1.verts[v];
            
            mesh1.add(node1);
            // mesh1.add(new arcsim::Node(pos,
            //                      pos,
            //                      pos,
            //                      0,false));
            // mesh1.verts[v]->node=mesh1.nodes[v];
            // arcsim::include(mesh1.verts[v],mesh1.nodes[v]->verts);
            // mesh1.nodes[v]->verts.push_back(mesh1.verts[v]);
    }


    // std::cout<<"Adding edges?\n";
    for (size_t e = 0; e<mesh->nEdges();e++){
        Edge e_orig= mesh->edge(e);

        arcsim::Edge *edge1 = new arcsim::Edge(mesh1.nodes[e_orig.firstVertex().getIndex()],mesh1.nodes[e_orig.secondVertex().getIndex()],0);
        mesh1.add(edge1);

    }
    
    



    // std::cout<<"Adding faces?\n";
    for(size_t f = 0; f< mesh->nFaces(); f++){
    Face f_orig = mesh->face(f);
    Halfedge he = f_orig.halfedge();
    arcsim::Face *face1 = new arcsim::Face(mesh1.verts[he.vertex().getIndex()],mesh1.verts[he.next().vertex().getIndex()],mesh1.verts[he.next().next().vertex().getIndex()],0,0 );        

    mesh1.add(face1);

    }

    // std::cout<<"computing data?\n";
    arcsim::compute_ms_data(mesh1);



    mesh1.numComponents = 1;
    
    

    // std::cout<<"done translating\n";
    return mesh1;
}



std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
translate_to_geometry(arcsim::Mesh mesh){

SimplePolygonMesh simpleMesh;

    // std::cout<<"This is being called\n";
//   processLoadedMesh(simpleMesh, loadType);
    Vector3 v_pos;
    
    // for(int v =0 ; v<mesh.nodes.size();v++){
    //     arcsim::Vec3 pos_old = mesh.nodes[v]->x;
    //     v_pos.x=pos_old[0];
    //     v_pos.y=pos_old[1];
    //     v_pos.z=pos_old[2];
    //     simpleMesh.vertexCoordinates.push_back(v_pos);
    // }
    bool flag_warning=false;
    // int flag=0;
    vector<int> flags(0);

    // double avg_neigh = 0;
    // double verts = 0;
    // std::cout<<"Something\n";
    // while(flag_warning){

    // arcsim::update_indices(mesh);
    // std::cout<<"The number of vertices is "<< mesh.verts.size()<<" \n";
    // flag_warning=false;
    // flags = vector<int>();
    // simpleMesh = SimplePolygonMesh();
    for(size_t v = 0 ; v < mesh.verts.size(); v++){
        arcsim::Vec3 pos_old = mesh.nodes[v]->x;
        
        v_pos.x=pos_old[0];
        v_pos.y=pos_old[1];
        v_pos.z=pos_old[2];
        // avg_neigh+=mesh.verts[v]->adjf.size();
        // verts+=1;
        // if(mesh.verts[v]->adjf.size()<=2){
        //     std::cout<<"The number of neighbors is "<< mesh.verts[v]->adjf.size()<<"\n";
        
    
            
        //     // for(int f_index = 0 ; f_index < mesh.verts[v]->adjf.size();f_index++)
        //     // { 
        //     // std::cout<<"Deleting face \n";
        //     // mesh.remove(mesh.verts[v]->adjf[f_index]);    
        //     // }
        //     // std::cout<<"Now we delete the vert\n";
        //     // mesh.remove(mesh.verts[v]);
            
        //     // // mesh.remove(mesh.verts[v]->adjf[0]);

        //     std::cout<<"The vertex index to not consider are"<< v<<" out of a total of"<< mesh.verts.size()<<"\n";
        //     flag_warning=true;
        //     // flag=v;
        //     flags.push_back(v);
        //     continue;
        // }
        
        
        simpleMesh.vertexCoordinates.push_back(v_pos);
    }

    // }
    int id1;
    int id2;
    int id3;
    // int flag_idx=0;
    // int number_of_flags=0;

    if(flag_warning){
    std::cout<<"THe number of flags is "<<flags.size()<<"\n";
    std::cout<<"THe flags are \n";
    for (size_t flag = 0 ; flag < flags.size(); flag++){
        std::cout<< flags[flag]<<"\t ";
    }
    std::cout<<" \n";
    }
    // std::cout<<"hihi\n";
    bool non_manifold = false;
    for(size_t f = 0 ; f < mesh.faces.size();f++){
        
        std::vector<size_t> polygon(3);
        
        id1 = mesh.faces[f]->v[0]->index;
        id2 = mesh.faces[f]->v[1]->index;
        id3 = mesh.faces[f]->v[2]->index;

        // int less_id1 = 0;
        // int less_id2 = 0;
        // int less_id3 = 0;

        // for(size_t flag = 0 ; flag < flags.size(); flag++){
        // if( id1 == flags[flag]|| id2 == flags[flag] || id3 == flags[flag]){
        //     non_manifold=true;
        // }
        // if(id1>flags[flag]&& flag_warning) less_id1+=1;
        // if(id2>flags[flag]&& flag_warning) less_id2+=1;
        // if(id3>flags[flag]&& flag_warning) less_id3+=1;
        // }
        // if(non_manifold){
        //     non_manifold=false;
        //     continue;
        // }
        // id1 = id1 - less_id1;
        // id2 = id2 - less_id2;
        // id3 = id3 - less_id3; 

        // if(id1==6075 || id2==6075 || id3 ==6075) std::cout<<" 2. This is is being called\n";

        // std::cout<<id1 <<" "<< id2 << " "<< id3 << "\n";
        polygon[0] = id1;
        polygon[1] = id2;
        polygon[2] = id3;
        
        simpleMesh.polygons.push_back(polygon);
        
    }
    
    // std::cout<<"Does this happen after loading the data to create the mesh?\n";
    // std::cout<<" THe information in the mesh is, "<< simpleMesh.vertexCoordinates.size()<<"number of vertices\n";
  auto lvals = makeManifoldSurfaceMeshAndGeometry(simpleMesh.polygons, simpleMesh.vertexCoordinates);
    if(flag_warning){
    std::cout<<"The problem is not the translation\n";
    }


 return std::tuple<std::unique_ptr<ManifoldSurfaceMesh>,
                    std::unique_ptr<VertexPositionGeometry>>(std::move(std::get<0>(lvals)),  // mesh
                                                             std::move(std::get<1>(lvals))); // geometry
}



int main(int argc, char** argv) {

    
    
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    auto start2 = chrono::steady_clock::now();
    auto end2 = chrono::steady_clock::now();
    srand((unsigned int) time(0));



    std::cout<< "Current path is " << argv[0]<<"\n";
    std::string basic_name = "../Results/Testing_Hessian/";

    // OK SO LETS get some testing going shal we 

    Eigen::Vector<double,9> Positions = Eigen::Vector<double,9>::Random(9);
      
    std::cout<<"This is a vector with 9 random components " << Positions.transpose() <<" \n" << "THe area of the triangle is" << geometry->Triangle_area(Positions) <<" \n"    ;


    // The next step is to calculate the derivatives wrt all the directions 
    

    double prev_val;
    double fwd_val;
    Eigen::Vector<double,9> Finite_Gradient = Eigen::Vector<double,9>::Zero(9);
    Eigen::Vector<double,9> Gradient = geometry->gradient_triangle_area(Positions);

    for( int dim = 0; dim < Positions.size(); dim++){

        Eigen::Vector<double,9> Positions2 = Positions;
        Positions2[dim]+=1e-6;
        fwd_val = geometry->Triangle_area(Positions2);
        Positions2[dim]-=2e-6;
        prev_val = geometry->Triangle_area(Positions2);
        Finite_Gradient[dim] = (fwd_val-prev_val)/(2e-6);

    }

    // std::cout<<"Now if we compare the area gradients we get"<<" \n";
    // for(int dim = 0; dim < Positions.size(); dim++){
    //     std::cout<< "The gradient in the direction " << dim << " is " << Gradient[dim] << " and the finite difference is " << Finite_Gradient[dim] << "\n";
    // }
    std::cout<<"The summed absolute difference between the area gradients is " << (Gradient-Finite_Gradient).cwiseAbs().sum() << "\n";

    // Lets see then 

    // std::cout<<"THe current area is " << geometry->Triangle_area(Positions) << "\n";
    // std::cout<<"IF we move with finite difference " << geometry->Triangle_area(Positions+Finite_Gradient*1e-6) << "\n";
    // std::cout<<"If we move with the gradient " << geometry->Triangle_area(Positions+Gradient*1e-6) << "\n";



    // Now we need to see if the cross product matrix works too.

    Eigen::Vector3d p1 = { Positions[0], Positions[1], Positions[2] };
    Eigen::Vector3d p2 = { Positions[3], Positions[4], Positions[5] };

    Eigen::Vector3d p1p2 = p1.cross(p2);
    Eigen::Vector3d p1p2Mat = geometry->Cross_product_matrix(p1)*p2;

    // std::cout<<"We are testing the cross product matrix, ";
    // std::cout<<"THe vectors are\n "<< p1p2 <<" \n" << "and \n " << p1p2Mat << " \n";
    // std::cout<<"the difference between the two vectors is " << (p1p2-p1p2Mat).norm() << "\n";



    // std::cout<<"Lets test the gradient of the edge-length now \n";

    Eigen::Vector<double,6> Edge_Positions = Eigen::Vector<double,6>::Random(6);

    // std::cout<< "THe edge length is" << geometry->Edge_length(Edge_Positions) <<" \n"    ;

    Eigen::Vector<double,6> Edge_Gradient = geometry->gradient_edge_length(Edge_Positions);

    Eigen::Vector<double,6> Finite_Gradient_Edge = Eigen::Vector<double,6>::Zero(6);

    for( int dim = 0; dim < Edge_Positions.size(); dim++){

        Eigen::Vector<double,6> Edge_Positions2 = Edge_Positions;
        Edge_Positions2[dim]+=1e-6;
        fwd_val = geometry->Edge_length(Edge_Positions2);
        Edge_Positions2[dim]-=2e-6;
        prev_val = geometry->Edge_length(Edge_Positions2);
        Finite_Gradient_Edge[dim] = (fwd_val-prev_val)/(2e-6);

    }

    // std::cout<< "Now if we compare the edge length gradients we get"<<" \n";
    // for(int dim = 0; dim < Edge_Positions.size(); dim++){
    //     std::cout<< "The gradient in the direction " << dim << " is " << Edge_Gradient[dim] << " and the finite difference is " << Finite_Gradient_Edge[dim] << "\n";
    // }
    std::cout<<"The summed absolute difference between the gradients of the edge length is " << (Edge_Gradient-Finite_Gradient_Edge).cwiseAbs().sum() << "\n";



    // We keep on testing the gradients

    Eigen::Vector<double,12> Dihedral_Positions = Eigen::Vector<double,12>::Random(12);

    // std::cout<< "THe dihedral angle is" << geometry->Dihedral_angle(Dihedral_Positions) <<" \n" ;

    Eigen::Vector<double,12> Dihedral_Gradient = geometry->gradient_dihedral_angle(Dihedral_Positions);
    
    Eigen::Vector<double,12> Finite_Gradient_Dihedral = Eigen::Vector<double,12>::Zero(12);

    for(size_t dim = 0; dim < Dihedral_Positions.size(); dim++){

        Eigen::Vector<double,12> Dihedral_Positions2 = Dihedral_Positions;
        Dihedral_Positions2[dim]+=1e-6;
        fwd_val = geometry->Dihedral_angle(Dihedral_Positions2);
        Dihedral_Positions2[dim]-=2e-6;
        prev_val = geometry->Dihedral_angle(Dihedral_Positions2);
        Finite_Gradient_Dihedral[dim] = (fwd_val-prev_val)/(2e-6);

    }
    
    std::cout<<"The summed absolute difference between the dihedral angle gradients is "<< (Dihedral_Gradient-Finite_Gradient_Dihedral).cwiseAbs().sum() << "\n";

    Eigen::Matrix<double, 6,6> Hessian_Edge = geometry->hessian_edge_length(Edge_Positions);

    // So i have the hessian 
    Eigen::Matrix<double, 6,6> Hessian_Edge_finite = Eigen::Matrix<double,6,6>::Zero(6,6);

    // Ok so i have the base, we need to calculate finite differences, we got the same number of variations as we have dims right?

    for(size_t dim = 0; dim < Edge_Positions.size(); dim++){

        Eigen::Vector<double,6> Edge_Positions2 = Edge_Positions;
        Edge_Positions2[dim]+=1e-6;
        Eigen::Vector<double,6> Gradient_fwd = geometry->gradient_edge_length(Edge_Positions2);
        Edge_Positions2[dim]-=2e-6;
        Eigen::Vector<double,6> Gradient_prev = geometry->gradient_edge_length(Edge_Positions2);
        Hessian_Edge_finite.col(dim) = (Gradient_fwd - Gradient_prev)/(2e-6);


    }


    // We will test our new volume function


    Eigen::Vector<double,9> Volume_Positions = Eigen::Vector<double,9>::Random(9);

    Eigen::Vector<double,9> Volume_gradient = geometry->gradient_volume(Volume_Positions);

    Eigen::Vector<double,9> Volume_gradient_finite = Eigen::Vector<double,9>::Zero(9);
    for(size_t dim = 0; dim < Volume_Positions.size(); dim++){

        Eigen::Vector<double,9> Volume_Positions2 = Volume_Positions;
        Volume_Positions2[dim]+=1e-6;
        fwd_val = geometry->Volume(Volume_Positions2);
        Volume_Positions2[dim]-=2e-6;
        prev_val = geometry->Volume(Volume_Positions2);
        Volume_gradient_finite[dim] = (fwd_val-prev_val)/(2e-6);

    }  

    std::cout<<"The summed absolute difference between the volume gradients is " << (Volume_gradient-Volume_gradient_finite).cwiseAbs().sum() << "\n";




    // std::cout<<"THe difference between the Hessians is \n";
    // for(size_t dim = 0; dim < EEdge_Positions.size(); dim++){
        std::cout<<"THe absolute difference between the Edge Hessians is " << (Hessian_Edge_finite-Hessian_Edge).sum() <<" \n";
        // std::cout<< Hessian_Edge_finite.col(dim) - Hessian_Edge.col(dim) << " \n";
    // }
    // std::cout<<"The hessian for the edge length is \n" << Hessian_Edge << " \n";

    // std::cout<<"THe hessian for the edge length with finite differences is \n" << Hessian_Edge_finite << " \n";




    Eigen::Matrix<double, 9, 9> Hessian_Area = geometry->hessian_triangle_area(Positions);
    Eigen::Matrix<double, 9,9> Hessian_Area_finite = Eigen::Matrix<double,9,9>::Zero(9,9);

    for(size_t dim = 0; dim < Positions.size(); dim++){

        Eigen::Vector<double,9> Positions2 = Positions;
        Positions2[dim]+=1e-6;
        Eigen::Vector<double,9> Gradient_fwd = geometry->gradient_triangle_area(Positions2);
        Positions2[dim]-=2e-6;
        Eigen::Vector<double,9> Gradient_prev = geometry->gradient_triangle_area(Positions2);
        Hessian_Area_finite.col(dim) = (Gradient_fwd - Gradient_prev)/(2e-6);
    }

    // std::cout<<"THe difference between the Hessians is \n";e
    std::cout<<"THe absolute difference between the Area Hessians is " << (Hessian_Area_finite-Hessian_Area).sum() <<" \n";
    // std::cout<<"The hessian of the Area is\n " << Hessian_Area << " \n";
    // std::cout<<"The hessian of the Area with finite differences is \n" << Hessian_Area_finite << " \n";


    Eigen::Matrix<double, 12, 12> Hessian_Dihedral = geometry->hessian_dihedral_angle(Dihedral_Positions);

    Eigen::Matrix<double, 12,12> Hessian_Dihedral_finite = Eigen::Matrix<double,12,12>::Zero(12,12);

    // Its time to finite difference this 

    for(size_t dim = 0; dim < Dihedral_Positions.size(); dim++){

        Eigen::Vector<double,12> Dihedral_Positions2 = Dihedral_Positions;
        Dihedral_Positions2[dim]+=1e-6;
        Eigen::Vector<double,12> Gradient_fwd = geometry->gradient_dihedral_angle(Dihedral_Positions2);
        Dihedral_Positions2[dim]-=2e-6;
        Eigen::Vector<double,12> Gradient_prev = geometry->gradient_dihedral_angle(Dihedral_Positions2);
        Hessian_Dihedral_finite.col(dim) = (Gradient_fwd - Gradient_prev)/(2e-6);

    }
    // std::cout<<"THe difference between the Hessians is \n";
    std::cout<<"THe absolute difference between the Dihedral Hessians is " << (Hessian_Dihedral_finite-Hessian_Dihedral).sum() <<" \n";
    // std::cout<<"The hessian of the Dihedral angle is\n " << Hessian_Dihedral << " \n";
    // std::cout<<"The hessian of the Dihedral angle with finite differences is \n" << Hessian_Dihedral_finite << " \n";


    Eigen::Matrix<double, 9,9> Hessian_Volume = geometry->hessian_volume(Volume_Positions);
    Eigen::Matrix<double, 9,9> Hessian_Volume_finite = Eigen::Matrix<double,9,9>::Zero(9,9);

    // Its time to finite difference this

    for(size_t dim = 0; dim < Volume_Positions.size(); dim++){

        Eigen::Vector<double,9> Volume_Positions2 = Volume_Positions;
        Volume_Positions2[dim]+=1e-6;
        Eigen::Vector<double,9> Gradient_fwd = geometry->gradient_volume(Volume_Positions2);
        Volume_Positions2[dim]-=2e-6;
        Eigen::Vector<double,9> Gradient_prev = geometry->gradient_volume(Volume_Positions2);
        Hessian_Volume_finite.col(dim) = (Gradient_fwd - Gradient_prev)/(2e-6);

    }

    std::cout<<"The absolute difference between the Volume Hessians is " << (Hessian_Volume_finite-Hessian_Volume).sum() <<" \n";



    std::cout<<"\n\n We start testing on a mesh now \n\n";

    std::string filepath = "../../../input/sphere.obj";
    // std::string filepath = "../Results/Mem3DG_Cell_Shape_KB_evol_flip/nu_0.625_c0_0.000_KA_10.000_KB_0.010000_init_cond_2_Nsim_11/Membrane_2067500.obj";
    // std::string filepath = "../../../input/20_icosahedron.obj";
    
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);



    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    std::vector<std::string> Energies(0);
    std::vector<std::vector<double>> Energy_constants(0);
    Energies.push_back("Surface_Tension");
    Energy_constants.push_back({10.0,0.0}); // just a dummy value for the surface tension
    
    E_Handler Sim_handler;
    Sim_handler = E_Handler(mesh, geometry, Energies, Energy_constants);

    Sim_handler.boundary = false;

    VertexData<Vector3> Force1(*mesh);
    VertexData<Vector3> Force2(*mesh);

    start = chrono::steady_clock::now();
    Force1 = Sim_handler.F_SurfaceTension(Energy_constants[0]);
    end = chrono::steady_clock::now();
    start2 = chrono::steady_clock::now();
    Force2 = Sim_handler.F_SurfaceTension_2(Energy_constants[0]);
    end2 = chrono::steady_clock::now();

    std::cout<<"The time for the original method is " << chrono::duration_cast<chrono::microseconds>(end - start).count() << " microseconds\n";
    std::cout<<"The time for the new method is " << chrono::duration_cast<chrono::microseconds>(end2 - start2).count() << " microseconds\n";

    // std::cout<<"The difference between the two methods is " << (Tension1-Tension2) << "\n";
    double abs_diff = 0.0;
    for(Vertex v: mesh->vertices()){
        abs_diff += (Force1[v]-Force2[v]).norm2();
    }
    std::cout<<"The absolute difference between the two tensions is " << abs_diff << "\n";

    start = chrono::steady_clock::now();
    Force1 = Sim_handler.F_Bending(Energy_constants[0]);
    end = chrono::steady_clock::now();

    start2 = chrono::steady_clock::now();
    Force2 = Sim_handler.F_Bending_2(Energy_constants[0]);
    end2 = chrono::steady_clock::now();

    std::cout<<"The time for the original bending method is " << chrono::duration_cast<chrono::microseconds>(end - start).count() << " microseconds\n";
    std::cout<<"The time for the new  bending method is " << chrono::duration_cast<chrono::microseconds>(end2 - start2).count() << " microseconds\n";

    abs_diff = 0.0;
    double abs_tot=0.0;
    for(Vertex v: mesh->vertices()){

        // std::cout<<"The relative difference is " << (Force1[v]-Force2[v]).norm()/Force1[v].norm() << "\n";
        // std::cout<<"Force 1 is " << Force1[v] << " and Force 2 is " << Force2[v] << "\n";
        abs_diff += (Force1[v]-Force2[v]).norm2();
        abs_tot +=Force1[v].norm2();
    }
    std::cout<<"The absolute difference between the two bending forces is " << abs_diff << "\n";
    // std::cout<<"The relative difference between the two bending forces is " << abs_diff/abs_tot << "\n";


    Energy_constants.push_back({1.0,4.5});
    
    start = chrono::steady_clock::now();
    Force1 = Sim_handler.F_Volume_constraint(Energy_constants[1]);
    end = chrono::steady_clock::now();
    start2 = chrono::steady_clock::now();
    Force2 = Sim_handler.F_Volume_constraint_2(Energy_constants[1]);
    end2 = chrono::steady_clock::now();
    std::cout<<"The time for the original volume constraint method is " << chrono::duration_cast<chrono::microseconds>(end - start).count() << " microseconds\n";
    std::cout<<"The time for the new volume constraint method is " << chrono::duration_cast<chrono::microseconds>(end2 - start2).count() << " microseconds\n";

    abs_diff = 0.0;
    for(Vertex v: mesh->vertices()){
        abs_diff += (Force1[v]-Force2[v]).norm2();
    }
    std::cout<<"The absolute difference between the two volume constraint forces is " << abs_diff << "\n";
    





    // OK so we want to do the same calculation now for the Hessian;

    const int Nverts = mesh->nVertices();

    // SparseMatrix<double> Hessian_bending = Sim_handler.H_Bending(Energy_constants[0]);

    // Eigen::MatrixXd Hessian_finite_diff(3*Nverts,3*Nverts);

    // VertexData<Vector3> Gradient_fwd(*mesh);
    // VertexData<Vector3> Gradient_prev(*mesh);

    // Gradient_fwd = Sim_handler.F_Bending(Energy_constants[0]);

    // // Gradient_fwd.array()

    // Eigen::VectorXd Difference(3*Nverts);


    // size_t dim = 0;
    // for(Vertex v: mesh->vertices()){
    //     // SO the idea here is  
    //     for(size_t coord = 0; coord < 3; coord++){
            
    //         geometry->inputVertexPositions[v][coord]+=1e-6; //move forward
    //         geometry->refreshQuantities();
    //         Gradient_fwd = Sim_handler.F_Bending_2(Energy_constants[0]);
    //         geometry->inputVertexPositions[v][coord]-=2e-6; //move backward
    //         geometry->refreshQuantities();
    //         Gradient_prev = Sim_handler.F_Bending_2(Energy_constants[0]);
    //         geometry->inputVertexPositions[v][coord]+=1e-6; //restore
    //         geometry->refreshQuantities();
    //         // Create the vector that goes in the column
    //         for(Vertex v2: mesh->vertices()){
    //             Difference[3*v2.getIndex()] = (Gradient_fwd[v2]-Gradient_prev[v2]).x/(2e-6);
    //             Difference[3*v2.getIndex()+1] = (Gradient_fwd[v2]-Gradient_prev[v2]).y/(2e-6);
    //             Difference[3*v2.getIndex()+2] = (Gradient_fwd[v2]-Gradient_prev[v2]).z/(2e-6); 
    //         }


    //         Hessian_finite_diff.col(dim) = Difference;
    //         dim+=1;


    //     }
    // }

    // // At this point i have both Hessians
    // // std::cout<<"The finite difference hessian is \n" << Hessian_finite_diff <<"\n";

    // std::cout<<"\n\n Testing Hessians now \n";
    // // std::cout<<"The sparse matrix i guess not so sparse\n" << Hessian_bending <<" \n";

    // // OK so this is the part where i compare my hessian with the finite difference one :P 

    // Eigen::MatrixXd DifferenceHessians = Hessian_finite_diff + Hessian_bending;

    // // std::cout<<"The matrix difference between the hessians is \n" << DifferenceHessians <<"\n";

    // // Ok lets do this 
    // std::cout<<"THe difference between the bending energy Hessians is " << DifferenceHessians.cwiseAbs().sum() <<" \n";

    // // Ok lets do our thing now we want the eigenvalues of this hessian

    // // Eigen::EigenSolver<Eigen::SPARSE> EIGENSOLVER;

    
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Hessian_finite_diff);

    // // std::cout << "The eigenvalues of the  hessian are :" 
    // //  << std::endl << es.eigenvalues() << std::endl;

    // //  Eigen::MatrixXcd EIGENVECS = es.eigenvectors();

    // // for(int i = 0 ; i < EIGENVECS.cols(); i++){
    // //     // NOW WE DO THE THINGY
    // //     double eig = es.eigenvalues()(i).real();
    // //     std::cout<<"The eigenvalue is "<< eig <<" \n";
    // //     std::cout<<"With eigenvector " << EIGENVECS.col(i).transpose() <<"\n";
    // // }
    
    // Eigen::MatrixXd Bending_H_analitical = Hessian_bending.toDense();
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ES(Bending_H_analitical);


    // // std::cout << "The eigenvalues of the hessian are :" 
    // //  << std::endl << ES.eigenvalues() << std::endl;
    

    // // Eigen::VectorXd Eigenvals_mat = EIGENSOLVER.eigenvalues().real();

    // // std::cout<<"The eigenvalues are " << Eigenvals_mat.transpose() <<"\n";





    // SparseMatrix<double> Hessian_tension = Sim_handler.H_SurfaceTension(Energy_constants[0]);

    // Eigen::MatrixXd Hessian_tension_finite_diff(3*Nverts,3*Nverts);





    // dim = 0;
    // for(Vertex v: mesh->vertices()){
    //     // SO the idea here is  
    //     for(size_t coord = 0; coord < 3; coord++){
            
    //         geometry->inputVertexPositions[v][coord]+=1e-6; //move forward
    //         geometry->refreshQuantities();
    //         Gradient_fwd = Sim_handler.F_SurfaceTension_2(Energy_constants[0]);
    //         geometry->inputVertexPositions[v][coord]-=2e-6; //move backward
    //         geometry->refreshQuantities();
    //         Gradient_prev = Sim_handler.F_SurfaceTension_2(Energy_constants[0]);
    //         geometry->inputVertexPositions[v][coord]+=1e-6; //restore
    //         geometry->refreshQuantities();
    //         // Create the vector that goes in the column
    //         for(Vertex v2: mesh->vertices()){
    //             Difference[3*v2.getIndex()] = (Gradient_fwd[v2]-Gradient_prev[v2]).x/(2e-6);
    //             Difference[3*v2.getIndex()+1] = (Gradient_fwd[v2]-Gradient_prev[v2]).y/(2e-6);
    //             Difference[3*v2.getIndex()+2] = (Gradient_fwd[v2]-Gradient_prev[v2]).z/(2e-6); 
    //         }


    //         Hessian_tension_finite_diff.col(dim) = Difference;
    //         dim+=1;


    //     }
    // }

    // // std::cout<<"The finite difference hessian is \n" << Hessian_tension_finite_diff <<"\n";
    // // std::cout<<"The analitical hessian is \n"<< Hessian_tension <<"\n";

    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ES_TFINITE(Hessian_tension_finite_diff);


    // // std::cout<< " THe eigenvalues of the finite difference hessian of the tension are " << std::endl << ES_TFINITE.eigenvalues() << std::endl;
    


    // Eigen::MatrixXd Tension_H_analitical = Hessian_tension.toDense();
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ES2(Tension_H_analitical);


    // // std::cout << "The eigenvalues of the hessian are :" 
    // //  << std::endl << ES2.eigenvalues() << std::endl;

    // // std::cout<<"The sparse matrix i guess not so sparse\n" << Hessian_tension <<" \n";

    // DifferenceHessians = Hessian_tension_finite_diff + Hessian_tension;

    // std::cout<<"THe difference between the surface tension energy Hessians is " << DifferenceHessians.cwiseAbs().sum() <<" \n";


    // // Lets do a mixture

    // // Eigen::MatrixXd B_T_H_analitical = Bending_H_analitical + Tension_H_analitical;
    // // Eigen::EigenSolver<Eigen::MatrixXd> ES3(B_T_H_analitical, false);
    // // std::cout << "The eigenvalues of the bending + tension hessian are :" 
    // //  << std::endl << ES3.eigenvalues() << std::endl;

    // //  Eigen::MatrixXcd EIGENVECSs = ES3.eigenvectors();

    // // for(int i = 0 ; i < EIGENVECSs.cols(); i++){
    //     // NOW WE DO THE THINGY
    //     // double eig = ES3.eigenvalues()(i).real();
    //     // if(eig < 1e-8 && eig > -1e-8){
    //         // std::cout<<"The eigenvalue is "<< eig <<" \n";
    //         // std::cout<<"With eigenvector " << EIGENVECSs.col(i) <<"\n";
    //     // }
    // // }


    // SparseMatrix<double> Hessian_Volume_mesh = Sim_handler.H_Volume(Energy_constants[0]);
    // Eigen::MatrixXd Hessian_finite_diff_vol(3*Nverts,3*Nverts);

    // dim = 0 ;
    // for(Vertex v: mesh->vertices()){
    //     for(size_t coord = 0; coord <3; coord++){
    //         geometry->inputVertexPositions[v][coord]+=1e-6; //move forward
    //         geometry->refreshQuantities();
    //         Gradient_fwd = Sim_handler.F_Volume(Energy_constants[0]);
    //         geometry->inputVertexPositions[v][coord]-=2e-6; //move backward
    //         geometry->refreshQuantities();
    //         Gradient_prev = Sim_handler.F_Volume(Energy_constants[0]);
    //         geometry->inputVertexPositions[v][coord]+=1e-6; //restore
    //         geometry->refreshQuantities();

    //         for(Vertex v2: mesh->vertices()){
    //             Difference[3*v2.getIndex()] = (Gradient_fwd[v2]-Gradient_prev[v2]).x/(2e-6);
    //             Difference[3*v2.getIndex()+1] = (Gradient_fwd[v2]-Gradient_prev[v2]).y/(2e-6);
    //             Difference[3*v2.getIndex()+2] = (Gradient_fwd[v2]-Gradient_prev[v2]).z/(2e-6); 
    //         }

    //         Hessian_finite_diff_vol.col(dim) = Difference;
    //         dim+=1;

    //     }
    // }
    // DifferenceHessians = Hessian_finite_diff_vol + Hessian_Volume_mesh;




    // std::cout<<"The difference between the volume energy Hessians is " << DifferenceHessians.cwiseAbs().sum() <<" \n";


    // std::ofstream Data_hess("../Results/Test_Hessian_2/Hessian_vol.txt");
    // Data_hess << Hessian_Volume_mesh;
    // Data_hess.close();
    // Data_hess = std::ofstream("../Results/Test_Hessian_2/Hessian_vol_finite.txt");
    // Data_hess << Hessian_finite_diff_vol;
    // Data_hess.close();

    // std::cout<<"THe volume hessian is \n" << Hessian_Volume_mesh<<" \n"; 
    
    // std::cout<<"The hessian of the finite difference is \n" << Hessian_finite_diff_vol <<" \n";

    arcsim::Cloth Cloth_1;
    arcsim::Cloth::Remeshing remeshing_params;

    remeshing_params.aspect_min = 0.2;
    remeshing_params.refine_angle = 0.7;
    remeshing_params.refine_compression =1.0;
    remeshing_params.refine_velocity = 1.0;
    remeshing_params.size_max = 0.5;
    remeshing_params.size_min = 0.001;
    remeshing_params.total_op = -1;

    // Now we do the fire test :P

    
    Energies.resize(0);
    Energy_constants.resize(0);
    std::vector<double> Constants(0);
    std::vector<std::string> Constraints(0);
    Eigen::VectorXd Lagrange_mults(4);
    Lagrange_mults(0) = 0.0;
    Lagrange_mults(1) = 0.0;
    Lagrange_mults(2) = 0.0;
    Lagrange_mults(3) = 0.0;
    

    Energies.push_back("Bending");
    Constants.push_back(10.0);
    Constants.push_back(0.0);
    Energy_constants.push_back(Constants);
    Energies.push_back("Area_constraint");
    Constants.resize(0);
    Constants.push_back(1000.0);
    double A_bar = geometry->totalArea();
    std::cout<<"The target area is "  << A_bar<<" and the area is" <<  geometry->totalArea() <<"  \n";
    Constants.push_back(A_bar);
    Energy_constants.push_back(Constants);

    Constraints.push_back("Volume");

    Constraints.push_back("CMx");
    Constraints.push_back("CMy");
    Constraints.push_back("CMz");
    // Energies.push_back()
    std::string Integration = "Newton";

    M3DG = Mem3DG(mesh,geometry);
    Sim_handler = E_Handler(mesh,geometry,Energies,Energy_constants);
    Sim_handler.Lagrange_mult = Lagrange_mults;
    Sim_handler.Constraints = Constraints;
    
    M3DG.Sim_handler = &Sim_handler;
    M3DG.recentering = false;
    M3DG.boundary = false;
    std::ofstream Some_ofstream;

    
    // arcsim::Mesh remesher_mesh = translate_to_arcsim(mesh,geometry);
    // Cloth_1.mesh = remesher_mesh;
    // Cloth_1.remeshing = remeshing_params;
    // arcsim::compute_masses(Cloth_1);
    // arcsim::compute_ws_data(Cloth_1.mesh);
    
    // arcsim::dynamic_remesh(Cloth_1);

    // std::tie(mesh_uptr, geometry_uptr) = translate_to_geometry(Cloth_1.mesh);
    // arcsim::delete_mesh(Cloth_1.mesh);
    // mesh = mesh_uptr.release();
    // geometry = geometry_uptr.release();
    // M3DG.mesh = mesh;
    // M3DG.geometry = geometry;
    // Sim_handler.mesh = mesh;
    // Sim_handler.geometry = geometry;
    
    std::cout<<"\t \t The volume is " << geometry->totalVolume() <<" \n";

    Sim_handler.Trgt_vol = geometry->totalVolume();
    // Sim_handler.Trgt_vol = 4.0;
    std::cout<<"The trgt vol is " << Sim_handler.Trgt_vol << "\n";

    std::cout<<"INtegrating newton\n";
    Save_mesh("../Results/Test_Hessian_3/",0);




    std::cout<<"Lets explore something\n";
    VertexData<Vector3> grad_vol = Sim_handler.F_Volume(std::vector<double>{-1.0});

    std::cout<<"The current volume is " << geometry->totalVolume() << " ";
    geometry->inputVertexPositions = geometry->inputVertexPositions + grad_vol*1e-3;
    geometry->refreshQuantities();
    std::cout<<" after moving forward is " << geometry->totalVolume() << " ";
    geometry->inputVertexPositions = geometry->inputVertexPositions - grad_vol*2e-3;
    geometry->refreshQuantities();
    std::cout<<" after moving backwards is " << geometry->totalVolume() << " ";

    geometry->inputVertexPositions = geometry->inputVertexPositions + grad_vol*1e-3;



    for(size_t step = 1; step < 20; step++){
    // M3DG.integrate(Some_ofstream,0.0,Energies,false);
    M3DG.integrate_Newton(Some_ofstream,0.0,Energies,false,Constraints);
    std::cout<<"THe area is " << geometry->totalArea()<< " and the volume is  " << geometry->totalVolume() << "\n";
    if(false){
        arcsim::Mesh remesher_mesh2 = translate_to_arcsim(mesh,geometry);
        Cloth_1.mesh=remesher_mesh2;
        Cloth_1.remeshing=remeshing_params;
        arcsim::compute_masses(Cloth_1);
        arcsim::compute_ws_data(Cloth_1.mesh);
        arcsim::dynamic_remesh(Cloth_1);

        delete mesh;
        delete geometry;

        std::tie(mesh_uptr, geometry_uptr) = translate_to_geometry(Cloth_1.mesh);
        arcsim::delete_mesh(Cloth_1.mesh);

        mesh = mesh_uptr.release();
        geometry = geometry_uptr.release();

        M3DG.mesh = mesh;
        M3DG.geometry = geometry;
        Sim_handler.mesh = mesh;
        Sim_handler.geometry = geometry;

    }


    // mkdir("../Results/Test_Hessian/", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    Save_mesh("../Results/Test_Hessian_3/",step);
    }



    std::string filename_basic = basic_name+"Output_data.txt";

    std::ofstream Sim_data;
    

    

    delete mesh;
    delete geometry;

    return EXIT_SUCCESS;
}



