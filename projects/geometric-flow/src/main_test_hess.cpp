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
#include "Interaction.h"


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


    
    std::cout<<"Lets do the gradient of the angle now\n";

    Eigen::Vector<double,9> Angle_positions;

    Angle_positions = Eigen::Vector<double,9>::Random(9);
    std::cout<<"Angle positions defined \n";
    std::cout<<"The angle is " << geometry->Angle(Angle_positions) << "\n";
    Eigen::Vector<double,9> Angle_gradient = geometry->gradient_angle(Angle_positions);
    Eigen::Vector<double,9> Finite_Gradient_Angle = Eigen::Vector<double,9>::Zero(9);

    for( int dim = 0; dim < Angle_positions.size(); dim++){

        Eigen::Vector<double,9> Angle_positions2 = Angle_positions;
        Angle_positions2[dim]+=1e-6;
        fwd_val = geometry->Angle(Angle_positions2);
        Angle_positions2[dim]-=2e-6;
        prev_val = geometry->Angle(Angle_positions2);
        Finite_Gradient_Angle[dim] = (fwd_val-prev_val)/(2e-6);

    }

    std::cout<<"Now if we compare the angle gradients we get"<<" \n";
    
    std::cout<<"The analitical gradient is " << Angle_gradient.transpose() << "\n";
    std::cout<<"The finite difference gradient is " << Finite_Gradient_Angle.transpose() << "\n";
    std::cout<<"The summed absolute difference between the angle gradients is " << (Angle_gradient-Finite_Gradient_Angle).cwiseAbs().sum() << "\n";

    // OK so we want to do the Hessian fo the angle 

    Eigen::Matrix<double, 9,9> Angle_hessian = geometry->hessian_angle(Angle_positions);
    Eigen::Matrix<double,9,9> Finite_Hessian_Angle = Eigen::Matrix<double, 9,9>::Zero(9,9);

    for( int dim = 0; dim < Angle_positions.size(); dim++){
        Eigen::Vector<double,9> Angle_positions2 = Angle_positions;
        Angle_positions2[dim]+=1e-6;
        Eigen::Vector<double,9> Gradient_fwd = geometry->gradient_angle(Angle_positions2);
        Angle_positions2[dim]-=2e-6;
        Eigen::Vector<double,9> Gradient_prev = geometry->gradient_angle(Angle_positions2);
        Finite_Hessian_Angle.col(dim) = (Gradient_fwd-Gradient_prev)/(2e-6);

    }

    std::cout<<"THe absolute difference between the angle Hessians is "<< (Angle_hessian - Finite_Hessian_Angle).cwiseAbs().sum() <<" \n";

    std::cout<<"THe angle hessian is \n" << Angle_hessian <<"\n";
    std::cout<<"THe finite dif hessian is \n" << Finite_Hessian_Angle<<"\n"; 




    // return 9;





    Eigen::Vector<double, 12> Cotan_w_positions;

    Cotan_w_positions = Eigen::Vector<double, 12>::Random(12);

    std::cout<<"The cotangent weight is " << geometry->Cotan_weight(Cotan_w_positions) << "\n";
    Eigen::Vector<double, 12> Cotan_w_gradient = geometry->gradient_cotan_weight(Cotan_w_positions);
    Eigen::Vector<double, 12> Finite_Gradient_Cotan_w = Eigen::Vector<double, 12>::Zero(12);

    for( int dim = 0; dim < Cotan_w_positions.size(); dim++){

        Eigen::Vector<double,12> Cotan_w_positions2 = Cotan_w_positions;
        Cotan_w_positions2[dim]+=1e-6;
        fwd_val = geometry->Cotan_weight(Cotan_w_positions2);
        Cotan_w_positions2[dim]-=2e-6;
        prev_val = geometry->Cotan_weight(Cotan_w_positions2);
        Finite_Gradient_Cotan_w[dim] = (fwd_val-prev_val)/(2e-6);

    }

    std::cout<<"Now if we compare the cotangent weight gradients we get"<<" \n";
    std::cout<<"The analitical gradient is " << Cotan_w_gradient.transpose() << "\n";
    std::cout<<"The finite difference gradient is " << Finite_Gradient_Cotan_w.transpose() << "\n";
    std::cout<<"The summed absolute difference between the cotangent weight gradients is " << (Cotan_w_gradient-Finite_Gradient_Cotan_w).cwiseAbs().sum() << "\n";


    Eigen::Matrix<double,12,12> Cotan_hessian = geometry->hessian_cotan_weight(Cotan_w_positions);
    Eigen::Matrix<double, 12,12> Finite_Hessian_Cotan = Eigen::Matrix<double,12,12>::Zero(9,9);

    for( int dim = 0; dim < Cotan_w_positions.size(); dim++){
        Eigen::Vector<double,12> Cotan_w_positions2 = Cotan_w_positions;
        Cotan_w_positions2[dim]+=1e-6;
        Eigen::Vector<double,12> Gradient_fwd = geometry->gradient_cotan_weight(Cotan_w_positions2);
        Cotan_w_positions2[dim]-=2e-6;
        Eigen::Vector<double,12> Gradient_prev = geometry->gradient_cotan_weight(Cotan_w_positions2);
        Finite_Hessian_Cotan.col(dim) = (Gradient_fwd-Gradient_prev)/(2e-6);
    }

    std::cout<<"The absolute difference between the cotan weights is "<< (Cotan_hessian - Finite_Hessian_Cotan).cwiseAbs().sum()<<"\n";

    std::cout<<"The cotan hessian is \n" << Cotan_hessian <<"\nThe finite dif hessian is \n"<< Finite_Hessian_Cotan<<"\n";
    // geometry->hessian_cotan_weight(Cotan_w_positions);

    // return 2;


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



    // We test this new energy

    Eigen::Vector<double,12> Edge_positions_reg = Eigen::Vector<double,12>::Random(12);
    Eigen::Vector<double,5> Edge_lenghts_red;

    Eigen::Vector3d p1_reg = { Edge_positions_reg[0], Edge_positions_reg[1], Edge_positions_reg[2] };
    Eigen::Vector3d p2_reg = { Edge_positions_reg[3], Edge_positions_reg[4], Edge_positions_reg[5] };
    Eigen::Vector3d p3_reg = { Edge_positions_reg[6], Edge_positions_reg[7], Edge_positions_reg[8] };
    Eigen::Vector3d p4_reg = { Edge_positions_reg[9], Edge_positions_reg[10], Edge_positions_reg[11] };
    Edge_lenghts_red[0] = (p1_reg-p2_reg).norm();
    Edge_lenghts_red[1] = (p1_reg-p3_reg).norm();
    Edge_lenghts_red[2] = (p1_reg-p4_reg).norm();
    Edge_lenghts_red[3] = (p2_reg-p3_reg).norm();
    Edge_lenghts_red[4] = (p2_reg-p4_reg).norm();

    std::cout<<"The edge lengths are " << Edge_lenghts_red.transpose() << "\n";
    std::cout<<"The edge length regularization energy is " << geometry->Ej_edge_regular(Edge_positions_reg) << "\n";

    Eigen::Vector<double,12> Edge_Gradient_reg = geometry->gradient_edge_regular(Edge_positions_reg);
    Eigen::Vector<double,12> Finite_Gradient_Edge_reg = Eigen::Vector<double,12>::Zero(12);

    for( int dim = 0; dim < Edge_positions_reg.size(); dim++){

        Eigen::Vector<double,12> Edge_positions_reg2 = Edge_positions_reg;
        Edge_positions_reg2[dim]+=1e-6;

        fwd_val = geometry->Ej_edge_regular(Edge_positions_reg2);
        Edge_positions_reg2[dim]-=2e-6;

        prev_val = geometry->Ej_edge_regular(Edge_positions_reg2);
        Finite_Gradient_Edge_reg[dim] = (fwd_val-prev_val)/(2e-6);

    }

    std::cout<<"Now if we compare the edge regularization gradients we get"<<" \n";
    std::cout<<"The analitical gradient is " << Edge_Gradient_reg.transpose() << "\n";
    std::cout<<"The finite difference gradient is " << Finite_Gradient_Edge_reg.transpose() << "\n";
    std::cout<<"The summed absolute difference between the edge regularization gradients is " << (Edge_Gradient_reg-Finite_Gradient_Edge_reg).cwiseAbs().sum() << "\n";


    // We will test the hessian of the Edge regularizer

    Eigen::Matrix<double, 12,12> Hessian_Edge_reg = geometry->hessian_edge_regular(Edge_positions_reg);
    Eigen::Matrix<double, 12,12> Hessian_Edge_reg_finite = Eigen::Matrix<double,12,12>::Zero(12,12);
    // So i have the base, we need to calculate finite differences, we got the same number of variations as we have dims right?
    for(long int dim = 0; dim < Edge_positions_reg.size(); dim++){

        Eigen::Vector<double,12> Edge_positions_reg2 = Edge_positions_reg;
        Edge_positions_reg2[dim]+=1e-6;
        Eigen::Vector<double,12> Gradient_fwd = geometry->gradient_edge_regular(Edge_positions_reg2);
        Edge_positions_reg2[dim]-=2e-6;
        Eigen::Vector<double,12> Gradient_prev = geometry->gradient_edge_regular(Edge_positions_reg2);
        Hessian_Edge_reg_finite.col(dim) = (Gradient_fwd - Gradient_prev)/(2e-6);
    }

    std::cout<<"THe absolute difference between the Edge regularization Hessians is " << (Hessian_Edge_reg_finite-Hessian_Edge_reg).cwiseAbs().sum() <<" \n";
    // std::cout<<"The hessian of the Edge regularization is \n" << Hessian_Edge_reg << " \n";
    // std::cout<<"The hessian of the Edge regularization with finite differences is \n" << Hessian_Edge_reg_finite << " \n";





    // We keep on testing the gradients

    Eigen::Vector<double,12> Dihedral_Positions = Eigen::Vector<double,12>::Random(12);

    // std::cout<< "THe dihedral angle is" << geometry->Dihedral_angle(Dihedral_Positions) <<" \n" ;

    Eigen::Vector<double,12> Dihedral_Gradient = geometry->gradient_dihedral_angle(Dihedral_Positions);
    
    Eigen::Vector<double,12> Finite_Gradient_Dihedral = Eigen::Vector<double,12>::Zero(12);

    for(long int dim = 0; dim < Dihedral_Positions.size(); dim++){

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

    for(long int dim = 0; dim < Edge_Positions.size(); dim++){

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
    for(long int dim = 0; dim < Volume_Positions.size(); dim++){

        Eigen::Vector<double,9> Volume_Positions2 = Volume_Positions;
        Volume_Positions2[dim]+=1e-6;
        fwd_val = geometry->Volume(Volume_Positions2);
        Volume_Positions2[dim]-=2e-6;
        prev_val = geometry->Volume(Volume_Positions2);
        Volume_gradient_finite[dim] = (fwd_val-prev_val)/(2e-6);

    }  

    std::cout<<"The summed absolute difference between the volume gradients is " << (Volume_gradient-Volume_gradient_finite).cwiseAbs().sum() << "\n";




    // std::cout<<"THe difference between the Hessians is \n";
    // for(long int dim = 0; dim < EEdge_Positions.size(); dim++){
        std::cout<<"THe absolute difference between the Edge Hessians is " << (Hessian_Edge_finite-Hessian_Edge).sum() <<" \n";
        // std::cout<< Hessian_Edge_finite.col(dim) - Hessian_Edge.col(dim) << " \n";
    // }
    // std::cout<<"The hessian for the edge length is \n" << Hessian_Edge << " \n";

    // std::cout<<"THe hessian for the edge length with finite differences is \n" << Hessian_Edge_finite << " \n";




    Eigen::Matrix<double, 9, 9> Hessian_Area = geometry->hessian_triangle_area(Positions);
    Eigen::Matrix<double, 9,9> Hessian_Area_finite = Eigen::Matrix<double,9,9>::Zero(9,9);

    for(long int dim = 0; dim < Positions.size(); dim++){

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

    for(long int dim = 0; dim < Dihedral_Positions.size(); dim++){

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

    for(long int dim = 0; dim < Volume_Positions.size(); dim++){

        Eigen::Vector<double,9> Volume_Positions2 = Volume_Positions;
        Volume_Positions2[dim]+=1e-6;
        Eigen::Vector<double,9> Gradient_fwd = geometry->gradient_volume(Volume_Positions2);
        Volume_Positions2[dim]-=2e-6;
        Eigen::Vector<double,9> Gradient_prev = geometry->gradient_volume(Volume_Positions2);
        Hessian_Volume_finite.col(dim) = (Gradient_fwd - Gradient_prev)/(2e-6);

    }

    std::cout<<"The absolute difference between the Volume Hessians is " << (Hessian_Volume_finite-Hessian_Volume).sum() <<" \n";


    // Lets test our new functinos

    Eigen::Vector<double,6 > R_positions = Eigen::Vector<double,6>::Random(6);
    Eigen::Vector<double,6> R_gradient = geometry->gradient_r(R_positions);
    Eigen::Vector<double,6> R_gradient_finite = Eigen::Vector<double,6>::Zero(6);
    for(long int dim = 0; dim < R_positions.size(); dim++){

        Eigen::Vector<double,6> R_positions2 = R_positions;
        R_positions2[dim]+=1e-6;
        fwd_val = geometry->r(R_positions2);
        R_positions2[dim]-=2e-6;
        prev_val = geometry->r(R_positions2);
        R_gradient_finite[dim] = (fwd_val-prev_val)/(2e-6);

    }
    std::cout<<"The summed absolute difference between the r gradients is " << (R_gradient-R_gradient_finite).cwiseAbs().sum() << "\n";
    std::cout<<"THe r gradient is " << R_gradient.transpose() << "\n";
    std::cout<<"THe r finite gradient is " << R_gradient_finite.transpose() << "\n";
    Eigen::Matrix<double, 6,6> Hessian_R = geometry->hessian_r(R_positions);
    Eigen::Matrix<double, 6,6> Hessian_R_finite = Eigen::Matrix<double,6,6>::Zero(6,6);
    // Its time to finite difference this
    for(long int dim = 0; dim < R_positions.size(); dim++){

        Eigen::Vector<double,6> R_positions2 = R_positions;
        R_positions2[dim]+=1e-6;
        Eigen::Vector<double,6> Gradient_fwd = geometry->gradient_r(R_positions2);
        R_positions2[dim]-=2e-6;
        Eigen::Vector<double,6> Gradient_prev = geometry->gradient_r(R_positions2);
        Hessian_R_finite.col(dim) = (Gradient_fwd - Gradient_prev)/(2e-6);

    }
    std::cout<<"The absolute difference between the Hessian of r is " << (Hessian_R_finite-Hessian_R).sum() <<" \n";
    // std::cout<<"The Hessian of r is \n" << Hessian_R << " \n";
    // std::cout<<"The Hessian of r with finite differences is \n" << Hessian_R_finite << " \n";

    // Now we do the next test for the triple product
    Eigen::Vector<double,13> Triple_positions = Eigen::Vector<double,13>::Random(13);
        
    for(int index = 0  ; index <3; index++){
        Triple_positions[12] = index;
        // std::cout<<"The index is" << Triple_positions[12] <<" \n";
        Eigen::Vector<double,12> Triple_gradient = geometry->gradient_triple_product(Triple_positions);
        Eigen::Vector<double,12> Triple_gradient_finite = Eigen::Vector<double,12>::Zero(12);
        
        std::cout<<"The triple product is " << geometry->Triple_product(Triple_positions)<<" \n";


        for(long int dim = 0; dim < 12; dim++){

            Eigen::Vector<double,13> Triple_positions2 = Triple_positions;
            // std::cout<<"THe index is " << Triple_positions2[12 ] << " \n";
            Triple_positions2[dim]+=1e-6;
            fwd_val = geometry->Triple_product(Triple_positions2);
            Triple_positions2[dim]-=2e-6;
            prev_val = geometry->Triple_product(Triple_positions2);
            Triple_gradient_finite[dim] = (fwd_val-prev_val)/(2e-6);

        }
        std::cout<<"The summed absolute difference between the triple product gradients is " << (Triple_gradient-Triple_gradient_finite).cwiseAbs().sum() << "\n";
        // std::cout<<"THe triple product gradient is " << Triple_gradient.transpose() << "\n";
        // std::cout<<"THe triple product finite gradient is " << Triple_gradient_finite.transpose() << "\n";

        // Lets do the hessian too
        Eigen::Matrix<double, 12,12> Hessian_Triple = geometry->hessian_triple_product(Triple_positions);
        Eigen::Matrix<double, 12,12> Hessian_Triple_finite = Eigen::Matrix<double,12,12>::Zero(12,12);

        for(long int dim = 0; dim < 12; dim++){
            Eigen::Vector<double,13> Triple_positions2 = Triple_positions;
            Triple_positions2[dim]+=1e-6;
            Eigen::Vector<double,12> fwd_val = geometry->gradient_triple_product(Triple_positions2);
            Triple_positions2[dim]-=2e-6;
            Eigen::Vector<double,12> prev_val = geometry->gradient_triple_product(Triple_positions2);
            Hessian_Triple_finite.col(dim) = (fwd_val-prev_val)/(2e-6);

        }

        std::cout<<"The absolute difference between the Hessian of triple prod is " << (Hessian_Triple-Hessian_Triple_finite).sum() <<" \n";
    

    }

    // we DO THE NEW EDGE REG

    Eigen::Vector<double,9> Positions_edge_reg = Eigen::Vector<double,9>::Random(9);
    Eigen::Vector<double,9> Edge_reg_gradient = geometry->gradient_edge_regular(Positions_edge_reg);
    Eigen::Vector<double,9> Edge_reg_gradient_finite = Eigen::Vector<double,9>::Zero(9);
    for(long int dim = 0; dim < Positions_edge_reg.size(); dim++){

        Eigen::Vector<double,9> Positions_edge_reg2 = Positions_edge_reg;
        Positions_edge_reg2[dim]+=1e-6;
        fwd_val = geometry->Ej_edge_regular(Positions_edge_reg2);
        Positions_edge_reg2[dim]-=2e-6;
        prev_val = geometry->Ej_edge_regular(Positions_edge_reg2);
        Edge_reg_gradient_finite[dim] = (fwd_val-prev_val)/(2e-6);

    }
    std::cout<<"The summed absolute difference between the edge regularization gradients is " << (Edge_reg_gradient-Edge_reg_gradient_finite).cwiseAbs().sum() << "\n";
    std::cout<<"THe edge regularization gradient is " << Edge_reg_gradient.transpose() << "\n";
    std::cout<<"THe edge regularization finite gradient is " << Edge_reg_gradient_finite.transpose() << "\n";

    Eigen::Matrix<double, 9,9> Edge_reg_Hessian = geometry->hessian_edge_regular(Positions_edge_reg);
    Eigen::Matrix<double, 9,9> Edge_reg_Hessian_finite = Eigen::Matrix<double,9,9>::Zero(9,9);
    // Its time to finite difference this
    for(long int dim = 0; dim < Positions_edge_reg.size(); dim++){

        Eigen::Vector<double,9> Positions_edge_reg2 = Positions_edge_reg;
        Positions_edge_reg2[dim]+=1e-6;
        Eigen::Vector<double,9> Gradient_fwd = geometry->gradient_edge_regular(Positions_edge_reg2);
        Positions_edge_reg2[dim]-=2e-6;
        Eigen::Vector<double,9> Gradient_prev = geometry->gradient_edge_regular(Positions_edge_reg2);
        Edge_reg_Hessian_finite.col(dim) = (Gradient_fwd - Gradient_prev)/(2e-6);

    }
    std::cout<<"The absolute difference between the Hessian of edge regularization is " << (Edge_reg_Hessian_finite-Edge_reg_Hessian).sum() <<" \n";
    std::cout<<"The Hessian of edge regularization is \n" << Edge_reg_Hessian  << " \n";
    std::cout<<"The Hessian of edge regularization with finite differences is \n" << Edge_reg_Hessian_finite << " \n";


    std::cout<<"\n\n We start testing on a mesh now \n\n";

    std::string filepath = "../../../input/4_tetrahedron.obj";
    // std::string filepath = "../Results/Mem3DG_Cell_Shape_KB_evol_flip/nu_0.625_c0_0.000_KA_10.000_KB_0.010000_init_cond_2_Nsim_11/Membrane_2067500.obj";
    // std::string filepath = "../../../input/20_icosahedron.obj";
    
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);



    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    std::vector<std::string> Energies(0);
    std::vector<std::vector<double>> Energy_constants(0);
    Energies.push_back("Surface_Tension");
    Energy_constants.push_back({1.0,1.0}); // just a dummy value for the surface tension
    
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
    


    // I want to test my Laplace Energy (:


    double E_Lap = Sim_handler.E_Laplace(Energy_constants[0]);
    double E_Lap_fwd;
    double E_lap_bkwd;
    std::cout<<"The Laplace Energy is " << E_Lap << "\n";

    VertexData<Vector3> Gradient_Laplace = Sim_handler.F_Laplace(Energy_constants[0]);
    VertexData<Vector3> Gradient_Laplace_finite(*mesh);
    VertexData<Vector3> Difference(*mesh);

    long int dim = 0;
    for(Vertex v: mesh->vertices()){
        // SO the idea here is  
        for(size_t coord = 0; coord < 3; coord++){
            
            geometry->inputVertexPositions[v][coord]+=1e-6; //move forward
            geometry->refreshQuantities();
            E_Lap_fwd = Sim_handler.E_Laplace(Energy_constants[0]);
            geometry->inputVertexPositions[v][coord]-=2e-6; //move backward
            geometry->refreshQuantities();
            E_lap_bkwd = Sim_handler.E_Laplace(Energy_constants[0]);
            geometry->inputVertexPositions[v][coord]+=1e-6; //restore
            geometry->refreshQuantities();

            Gradient_Laplace_finite[v][coord] = (E_Lap_fwd - E_lap_bkwd)/(2e-6);
            
            }
        
        Difference[v] = Gradient_Laplace[v] + Gradient_Laplace_finite[v];
        // std::cout<<"The gradient f            or vertex " << v.getIndex() << " is " << Gradient_Laplace[v] << "\n";
        // std::cout<<"The finite diff gradient for vertex " << v.getIndex() << " is " << Gradient_Laplace_finite[v] << "\n";
        std::cout<<"The difference for vertex           " << v.getIndex() << " is " << Difference[v] << "\n";

    }
    // std::cout<<"The gradient of the Laplace energy is \n" << Gradient_Laplace << "\n";
    // std::cout<<"The finite difference gradient of the Laplace energy is \n" << Gradient_Laplace_finite << "\n";
    // std::cout<<"The difference between gradients is \n" << Difference << "\n";

    // return 1;

    double E_reg = Sim_handler.E_Edge_reg(Energy_constants[0]);
    double E_reg_fwd;
    double E_reg_bkwd;

    std::cout<<"The Edge Regularization Energy is " << E_reg << "\n";
    VertexData<Vector3> Gradient_Edge_reg = Sim_handler.F_Edge_reg(Energy_constants[0]);
    VertexData<Vector3> Gradient_Edge_reg_finite(*mesh);

    for(Vertex v: mesh->vertices()){
        // SO the idea here is  
        for(size_t coord = 0; coord < 3; coord++){
            
            geometry->inputVertexPositions[v][coord]+=1e-6; //move forward
            geometry->refreshQuantities();
            E_reg_fwd = Sim_handler.E_Edge_reg(Energy_constants[0]);
            geometry->inputVertexPositions[v][coord]-=2e-6; //move backward
            geometry->refreshQuantities();
            E_reg_bkwd = Sim_handler.E_Edge_reg(Energy_constants[0]);
            geometry->inputVertexPositions[v][coord]+=1e-6; //restore
            geometry->refreshQuantities();

            Gradient_Edge_reg_finite[v][coord] = (E_reg_fwd - E_reg_bkwd)/(2e-6);
            
        }
        // std::cout<<"The gradient for vertex " << v.getIndex() << " is " << Gradient_Edge_reg[v] << "\n";
        // std::cout<<"The finite diff gradient for vertex " << v.getIndex() << " is " << Gradient_Edge_reg_finite[v] << "\n";
    }




    // std::cout<<"The gradient of the Edge Regularization energy is \n" << Gradient_Edge_reg << "\n";
    // std::cout<<"The finite difference gradient of the Edge Regularization energy is \n" << Gradient_Edge_reg_finite << "\n";
    // std::cout<<"The difference between gradients is \n" << (Gradient_Edge_reg - Gradient_Edge_reg_finite) << "\n";
    std::cout<<"Lets check the hessian of the edge regularization in the mesh\n";

    SparseMatrix<double> Hessian_Edge_reg_mesh = Sim_handler.H_Edge_reg(Energy_constants[0]);
    Eigen::MatrixXd Hessian_Edge_reg_finite_mesh(3*mesh->nVertices(),3*mesh->nVertices());
    VertexData<Vector3> Gradient_fwd(*mesh);
    VertexData<Vector3> Gradient_prev(*mesh);
    
    Eigen::VectorXd Difference_grad(3*mesh->nVertices());
    dim = 0;

  
    for(Vertex v: mesh->vertices()){
        // SO the idea here is  
        for(size_t coord = 0; coord < 3; coord++){
            
            geometry->inputVertexPositions[v][coord]+=1e-6; //move forward
            geometry->refreshQuantities();
            Gradient_fwd = Sim_handler.F_Edge_reg(Energy_constants[0]);
            geometry->inputVertexPositions[v][coord]-=2e-6; //move backward
            geometry->refreshQuantities();
            Gradient_prev = Sim_handler.F_Edge_reg(Energy_constants[0]);
            geometry->inputVertexPositions[v][coord]+=1e-6; //restore
            geometry->refreshQuantities();
            // Create the vector that goes in the column
            for(Vertex v2: mesh->vertices()){
                Difference_grad[3*v2.getIndex()] = (Gradient_fwd[v2]-Gradient_prev[v2]).x/(2e-6);
                Difference_grad[3*v2.getIndex()+1] = (Gradient_fwd[v2]-Gradient_prev[v2]).y/(2e-6);
                Difference_grad[3*v2.getIndex()+2] = (Gradient_fwd[v2]-Gradient_prev[v2]).z/(2e-6); 
            }


            Hessian_Edge_reg_finite_mesh.col(dim) = Difference_grad;
            dim+=1;


        }
    }
    
     std::cout<<"The finite difference hessian is \n" << Hessian_Edge_reg_finite_mesh <<"\n";

    // std::cout<<"\n\n Testing Hessians now \n";
    std::cout<<"The sparse matrix i guess not so sparse\n" << Hessian_Edge_reg_mesh <<" \n";

    // OK so this is the part where i compare my hessian with the finite difference one :P 

    // Eigen::MatrixXd DifferenceHessians = Hessian_finite_diff + Hessian_bending;
    Eigen::MatrixXd DifferenceHessians = Hessian_Edge_reg_finite_mesh + Hessian_Edge_reg_mesh;
    std::cout<<"The matrix difference between the hessians is \n" << DifferenceHessians <<"\n";                                                                                     


    std::cout<<"Lets check the Hessian of laplace\n";

    SparseMatrix<double> Hessian_Laplace_mesh = Sim_handler.H_Laplace(Energy_constants[0]);
    Eigen::MatrixXd Hessian_Laplace_finite_mesh(3*mesh->nVertices(),3*mesh->nVertices());
    // VertexData<Vector3> Gradient_fwd(*mesh);
    // VertexData<Vector3> Gradient_prev(*mesh);
    // Eigen::VectorXd Difference_grad_laplace(3*mesh->nVertices());
    dim = 0;
    for(Vertex v: mesh->vertices()){
        for(size_t coord = 0; coord <3; coord++){
            geometry->inputVertexPositions[v][coord]+=1e-6;
            geometry->refreshQuantities();
            Gradient_fwd = Sim_handler.F_Laplace(Energy_constants[0]);

            geometry->inputVertexPositions[v][coord]-=2e-6;
            geometry->refreshQuantities();
            Gradient_prev = Sim_handler.F_Laplace(Energy_constants[0]);

            geometry->inputVertexPositions[v][coord]+=1e-6;
            geometry->refreshQuantities();

            for(Vertex v2: mesh->vertices()){
                Difference_grad[3*v2.getIndex()] = (Gradient_fwd[v2]-Gradient_prev[v2]).x/(2e-6);
                Difference_grad[3*v2.getIndex()+1] = (Gradient_fwd[v2]-Gradient_prev[v2]).y/(2e-6);
                Difference_grad[3*v2.getIndex()+2] = (Gradient_fwd[v2]-Gradient_prev[v2]).z/(2e-6); 
            }

            Hessian_Laplace_finite_mesh.col(dim) = Difference_grad;
            dim+=1;


        }

    }
    
    std::cout<<"The finite difference hessian is \n"<< Hessian_Laplace_finite_mesh <<"\nThe sparse matrix i guess not so sparse\n" << Hessian_Laplace_mesh<<"\n";

    DifferenceHessians = Hessian_Laplace_mesh + Hessian_Laplace_finite_mesh;

    // std::cout<<"The matrix difference between the hessians is \n"<< DifferenceHessians<<"\n";

    SparseMatrix<double> Bending_hessian = Sim_handler.H_Bending(Energy_constants[0]);

    std::cout<<"THe bending hessian is \n" << Bending_hessian << "\n";




    // Lets check the face reg gradient

    Sim_handler.update_face_reference();

    // 


    double E_face = Sim_handler.E_Face_reg(Energy_constants[0]);
    double E_face_fwd;
    double E_face_bkwd;

    std::cout<<"The Face Regularization Energy is " << E_face << "\n";
    VertexData<Vector3> Gradient_Face_reg = Sim_handler.F_Face_reg(Energy_constants[0]);
    VertexData<Vector3> Gradient_Face_reg_finite(*mesh,{0.0,0.0,0.0});
    VertexData<Vector3> Difference_Face_reg(*mesh);

    for(Vertex v: mesh->vertices()){
        // SO the idea here is  
        for(size_t coord = 0; coord < 3; coord++){
            
            geometry->inputVertexPositions[v][coord]+=1e-6; //move forward
            geometry->refreshQuantities();
            E_face_fwd = Sim_handler.E_Face_reg(Energy_constants[0]);
            geometry->inputVertexPositions[v][coord]-=2e-6; //move backward
            geometry->refreshQuantities();
            E_face_bkwd = Sim_handler.E_Face_reg(Energy_constants[0]);
            geometry->inputVertexPositions[v][coord]+=1e-6; //restore
            geometry->refreshQuantities();

            Gradient_Face_reg_finite[v][coord] = (E_face_fwd - E_face_bkwd)/(2e-6);
            
        }
        Difference_Face_reg[v] = Gradient_Face_reg[v] + Gradient_Face_reg_finite[v];
        std::cout<<"The difference for vertex           " << v.getIndex() << " is " << Difference_Face_reg[v] << "\n";
        // std::cout<<"The gradient for vertex             " << v.getIndex() << " is " << Gradient_Face_reg[v] << "\n";
        // std::cout<<"The finite diff gradient for vertex " << v.getIndex() << " is " << Gradient_Face_reg_finite[v] << "\n";

    }

    // Lets do the hessian now

    SparseMatrix<double> Hessian_Face_reg = Sim_handler.H_Face_reg(Energy_constants[0]);
    std::cout << "The Hessian of face regularization is \n" << Hessian_Face_reg << "\n";

    // Compare the Hessian of face regularization with the finite difference Hessian

    Eigen::MatrixXd Hessian_Face_reg_finite_mesh(3*mesh->nVertices(),3*mesh->nVertices());

    dim = 0;
    for(Vertex v: mesh->vertices()){
        for(size_t coord = 0; coord <3; coord++){
            geometry->inputVertexPositions[v][coord]+=1e-6;
            geometry->refreshQuantities();
            Gradient_fwd = Sim_handler.F_Face_reg(Energy_constants[0]);

            geometry->inputVertexPositions[v][coord]-=2e-6;
            geometry->refreshQuantities();
            Gradient_prev = Sim_handler.F_Face_reg(Energy_constants[0]);

            geometry->inputVertexPositions[v][coord]+=1e-6;
            geometry->refreshQuantities();

            for(Vertex v2: mesh->vertices()){
                Difference_grad[3*v2.getIndex()] = (Gradient_fwd[v2]-Gradient_prev[v2]).x/(2e-6);
                Difference_grad[3*v2.getIndex()+1] = (Gradient_fwd[v2]-Gradient_prev[v2]).y/(2e-6);
                Difference_grad[3*v2.getIndex()+2] = (Gradient_fwd[v2]-Gradient_prev[v2]).z/(2e-6); 
            }

            Hessian_Face_reg_finite_mesh.col(dim) = Difference_grad;
            dim+=1;


        }

    }


    Eigen::MatrixXd Difference_Hessian_Face_reg = Hessian_Face_reg + Hessian_Face_reg_finite_mesh;
    std::cout << "The matrix difference between the face regularization Hessians is \n" << Difference_Hessian_Face_reg << "\n";




    return 1;




    // We want to test the Energies now

    Interaction* Bead_I;
    double rc = 2.0;
    std::vector<double> params(0);
    params.push_back(1.0); // epsilon
    params.push_back(1.0); // sigma
    params.push_back(rc); // cutoff

    // Bead_I = new Constant_Normal(params);
    // Bead_I = new Frenkel_Normal(params);
    Bead_I = new Frenkel_Normal(params);

    double r = 1.7325;
    double r_2;
    // I want do debug the Frenkel Normal interaction

    double E_fwd;
    double E_prev;

    double Grad_E = Bead_I->dE_r(r,params);

    r_2 = r;
    r_2 += 1e-6; // move forward
    E_fwd = Bead_I->E_r(r_2,params);
    r_2 -= 2e-6; // move backward
    E_prev = Bead_I->E_r(r_2,params);

    double Finite_diff_grad = (E_fwd - E_prev)/(2e-6);
    std::cout<<"The gradient of the Frenkel Normal interaction is " << Grad_E << "\n";
    std::cout<<"The finite difference gradient of the Frenkel Normal interaction is " << Finite_diff_grad << "\n";

    // First thing is to 

    double DD_E = Bead_I->ddE_r(r,params);

    double d_E_fwd;
    double d_E_prev;
    r_2 = r;
    r_2 += 1e-6; // move forward
    d_E_fwd = Bead_I->dE_r(r_2,params);
    r_2 -= 2e-6; // move backward
    d_E_prev = Bead_I->dE_r(r_2,params);
    double Finite_diff_dE = (d_E_fwd - d_E_prev)/(2e-6);
    std::cout<<"The second derivative of the Frenkel Normal interaction is " << DD_E << "\n";
    std::cout<<"The finite difference second derivative of the Frenkel Normal interaction is " << Finite_diff_dE << "\n";



    // Lets do an experiment here 

    double E_over_r = Bead_I->E_r(r,params)/r;

    double dE_over_r = Bead_I->dE_r(r,params)/r - Bead_I->E_r(r,params)/(r*r);

    double dE_over_r_finite;

    r_2 = r;
    r_2 += 1e-6; // move forward
    E_fwd = Bead_I->E_r(r_2,params)/r_2;
    r_2 -= 2e-6; // move backward
    E_prev = Bead_I->E_r(r_2,params)/r_2;
    dE_over_r_finite = (E_fwd - E_prev)/(2e-6);

    std::cout<<"The energy over r is " << E_over_r << "\n";
    std::cout<<"The derivative of the energy over r is " << dE_over_r << "\n";
    std::cout<<"The finitative of the energy over r is " << dE_over_r_finite << "\n";




    // I need to create a Bead now

    std::cout<<"Creating bead\n";
    int N_beads =1;
    Bead Bead_1(mesh,geometry,Vector3{2.0,0.0,0.0},params,Bead_I,0,N_beads);
    std::cout<<"Calculating bead energy\n";
    Bead_I->mesh = mesh;
    Bead_I->geometry = geometry;
    Bead_I->Bead_1 = &Bead_1;

    Bead_1.interaction = "Shifted_LJ_Normal_nopush";
    Bead_1.rc = 2.0;
    Bead_1.sigma = 1.0;
    Bead_1.strength = 1.0;

    VertexData<Vector3> Force_old_method = Bead_1.Gradient();
    Vector3 old_Bead_F = Bead_1.Total_force;
    Bead_1.Total_force = Vector3{0.0,0.0,0.0};

    std::cout<<"The bead energy is" << Bead_1.Bead_I->Tot_Energy() << " \n";
    std::cout<<"The bead energyold" << Bead_1.Energy()<<"\n";
    VertexData<Vector3> Inter_Force = Bead_1.Bead_I->Gradient();
    VertexData<Vector3> Inter_Force_finite(*mesh,{0.0,0.0,0.0});

    double fwd_val_x;
    double prev_val_x;
    double fwd_val_y;
    double prev_val_y;
    double fwd_val_z;
    double prev_val_z;
    // Now we want to calculate the finite difference gradient of the bead interaction
    for(Vertex v: mesh->vertices()){
        // SO the idea here is  
            // std::cout<<"THe vertex is " << v.getIndex() << "\n";
            // std::cout<<"Position is " << geometry->inputVertexPositions[v] << "\n";

            
            geometry->inputVertexPositions[v].x+=1e-7; //move forward
            geometry->refreshQuantities();
            fwd_val_x = Bead_1.Bead_I->Tot_Energy();
            geometry->inputVertexPositions[v].x-=2e-7; //move backward
            geometry->refreshQuantities();
            prev_val_x = Bead_1.Bead_I->Tot_Energy();
            geometry->inputVertexPositions[v].x+=1e-7; //restore
            geometry->refreshQuantities();

            // std::cout<< setprecision(9)<<"The contribution of vertex to the bead energy is " << fwd_val_x << "and " << prev_val_x << "\n";

            geometry->inputVertexPositions[v].y+=1e-7; //move forward
            geometry->refreshQuantities();
            fwd_val_y = Bead_1.Bead_I->Tot_Energy();
            geometry->inputVertexPositions[v].y-=2e-7; //move backward
            geometry->refreshQuantities();
            prev_val_y = Bead_1.Bead_I->Tot_Energy();
            geometry->inputVertexPositions[v].y+=1e-7; //restore
            geometry->refreshQuantities();

            // std::cout<< setprecision(9) <<"The contribution of vertex to the bead energy is " << fwd_val_y << "and " << prev_val_y << "\n";


            geometry->inputVertexPositions[v].z+=1e-7; //move forward
            geometry->refreshQuantities();
            fwd_val_z = Bead_1.Bead_I->Tot_Energy();
            geometry->inputVertexPositions[v].z-=2e-7; //move backward
            geometry->refreshQuantities();
            prev_val_z = Bead_1.Bead_I->Tot_Energy();
            geometry->inputVertexPositions[v].z+=1e-7; //restore
            geometry->refreshQuantities();

            // std::cout<< setprecision(9)<<"The contribution of vertex to the bead energy is " << fwd_val_y << "and " << prev_val_y << "\n";


            Inter_Force_finite[v].x = (fwd_val_x - prev_val_x)/(2e-7);
            Inter_Force_finite[v].y = (fwd_val_y - prev_val_y)/(2e-7);
            Inter_Force_finite[v].z = (fwd_val_z - prev_val_z)/(2e-7);

            

        // std::cout<<"The gradient for vertex " << v.getIndex() << " is " << Inter_Force[v] << "\n";
        // std::cout<<"The finite diff gradient for vertex " << v.getIndex() << " is " << Inter_Force_finite[v] << "\n";
    }
    // std::cout<<"We calculated the gradient for the bead interaction\n";

    for(Vertex v: mesh->vertices()){
        // std::cout<<"The position is" << geometry->inputVertexPositions[v] << "\n";
        std::cout<<"The gradient for vertex " << v.getIndex() << " is " << Inter_Force[v] << "\n";
        // std::cout<<"The olgradit for vertex " << v.getIndex() << " is " << Force_old_method[v] << "\n";
        std::cout<<"The finigrad for vertex " << v.getIndex() << " is " << Inter_Force_finite[v] << "\n";
    }
    std::cout<<"The gradient for the Bead is "  << Bead_1.Total_force << "\n";

    Bead_1.Pos.x += 1e-6;
    fwd_val_x = Bead_1.Bead_I->Tot_Energy();
    Bead_1.Pos.x -= 2e-6; //move backward
    prev_val_x = Bead_1.Bead_I->Tot_Energy();
    Bead_1.Pos.x += 1e-6; //restore

    Bead_1.Pos.y += 1e-6;
    fwd_val_y = Bead_1.Bead_I->Tot_Energy();
    Bead_1.Pos.y -= 2e-6; //move backward
    prev_val_y = Bead_1.Bead_I->Tot_Energy();
    Bead_1.Pos.y += 1e-6; //restore

    Bead_1.Pos.z += 1e-6;
    fwd_val_z = Bead_1.Bead_I->Tot_Energy();
    Bead_1.Pos.z -= 2e-6; //move backward
    prev_val_z = Bead_1.Bead_I->Tot_Energy();
    Bead_1.Pos.z += 1e-6; //restore

    Vector3 Finite_diff_Bead_grad;
    Finite_diff_Bead_grad.x = (fwd_val_x - prev_val_x)/(2e-6);
    Finite_diff_Bead_grad.y = (fwd_val_y - prev_val_y)/(2e-6);
    Finite_diff_Bead_grad.z = (fwd_val_z - prev_val_z)/(2e-6);
    std::cout<<"THe oldgradi for the Bead is " << old_Bead_F << "\n";
    std::cout<<"The finigrad for the Bead is " << Finite_diff_Bead_grad << "\n";

    for(Vertex v: mesh->vertices()){
        std::cout<<"The position is" << geometry->inputVertexPositions[v] << "\n";
    }


    


    // Now we have a full Potential tested. WHat was the next step? 
    // Right, the Hessian
    
    // We are testing the Hessian now, this is insane

    std::cout<<"Testing the Hessian now\n";

    std::cout<<"The bead id is" << Bead_1.Bead_id << "\n";

    Eigen::MatrixXd Hessian_Bead = Bead_1.Bead_I->Hessian();

    // std::cout<<"The Hessian of the Bead interaction is \n" << Hessian_Bead << "\n";

    Eigen::MatrixXd Hessian_Bead_finite = Eigen::MatrixXd::Zero(3*mesh->nVertices()+3*N_beads,3*mesh->nVertices()+3*N_beads);

    // Its time to finite difference this
    

    dim = 0;
    Eigen::VectorXd Difference_grad_B(3*mesh->nVertices()+3*N_beads);
    Vector3 Force_fwd;
    Vector3 Force_prev;
    for(Vertex v: mesh->vertices()){
        // SO the idea here is  
        for(size_t coord = 0; coord < 3; coord++){
            
            geometry->inputVertexPositions[v][coord]+=1e-6; //move forward
            geometry->refreshQuantities();
            Bead_1.Total_force= Vector3{0.0,0.0,0.0}; // reset the force
            Gradient_fwd = Bead_1.Bead_I->Gradient();
            Force_fwd = Bead_1.Total_force;

            Bead_1.Total_force= Vector3{0.0,0.0,0.0}; // reset the force
            geometry->inputVertexPositions[v][coord]-=2e-6; //move backward
            geometry->refreshQuantities();
            Gradient_prev = Bead_1.Bead_I->Gradient();
            Force_prev = Bead_1.Total_force;
            geometry->inputVertexPositions[v][coord]+=1e-6; //restore
            geometry->refreshQuantities();
            // Create the vector that goes in the column
            for(Vertex v2: mesh->vertices()){
                Difference_grad_B[3*v2.getIndex()] = (Gradient_fwd[v2]-Gradient_prev[v2]).x/(2e-6);
                Difference_grad_B[3*v2.getIndex()+1] = (Gradient_fwd[v2]-Gradient_prev[v2]).y/(2e-6);
                Difference_grad_B[3*v2.getIndex()+2] = (Gradient_fwd[v2]-Gradient_prev[v2]).z/(2e-6); 
            }
            Difference_grad_B[3*mesh->nVertices()+3*Bead_1.Bead_id+0] = (Force_fwd.x - Force_prev.x)/(2e-6);
            Difference_grad_B[3*mesh->nVertices()+3*Bead_1.Bead_id+1] = (Force_fwd.y - Force_prev.y)/(2e-6);
            Difference_grad_B[3*mesh->nVertices()+3*Bead_1.Bead_id+2] = (Force_fwd.z - Force_prev.z)/(2e-6);


            Hessian_Bead_finite.col(dim) = Difference_grad_B;
            dim+=1;


        }
        }
        // And u also need to do the one for the bead
        for(size_t coord = 0; coord < 3; coord++){
            
            Bead_1.Pos[coord]+=1e-6; //move forward
            geometry->refreshQuantities();
            Bead_1.Total_force= Vector3{0.0,0.0,0.0}; // reset the force

            Gradient_fwd = Bead_1.Bead_I->Gradient();
            Force_fwd = Bead_1.Total_force;
            Bead_1.Total_force= Vector3{0.0,0.0,0.0}; // reset the force
            Bead_1.Pos[coord]-=2e-6; //move backward
            geometry->refreshQuantities();
            Gradient_prev = Bead_1.Bead_I->Gradient();
            Force_prev = Bead_1.Total_force;

            Bead_1.Total_force =  Vector3{0.0,0.0,0.0};; // reset the force
            Bead_1.Pos[coord]+=1e-6; //restore
            geometry->refreshQuantities();
            // Create the vector that goes in the column
            for(Vertex v2: mesh->vertices()){
                Difference_grad_B[3*v2.getIndex()] = (Gradient_fwd[v2]-Gradient_prev[v2]).x/(2e-6);
                Difference_grad_B[3*v2.getIndex()+1] = (Gradient_fwd[v2]-Gradient_prev[v2]).y/(2e-6);
                Difference_grad_B[3*v2.getIndex()+2] = (Gradient_fwd[v2]-Gradient_prev[v2]).z/(2e-6); 
            }
            Difference_grad_B[3*mesh->nVertices()+3*Bead_1.Bead_id+0] = (Force_fwd.x - Force_prev.x)/(2e-6);
            Difference_grad_B[3*mesh->nVertices()+3*Bead_1.Bead_id+1] = (Force_fwd.y - Force_prev.y)/(2e-6);
            Difference_grad_B[3*mesh->nVertices()+3*Bead_1.Bead_id+2] = (Force_fwd.z - Force_prev.z)/(2e-6);

            Hessian_Bead_finite.col(dim) = Difference_grad_B;
            dim+=1;
        }


    
    

    // std::cout<<"THe Hessian of the Bead interaction with finite differences is \n" << Hessian_Bead_finite << "\n";
    std::cout<<"The summ of the Hessians is" << (Hessian_Bead_finite+Hessian_Bead).cwiseAbs().sum() << "\n";

    Energies.resize(0);
    Energy_constants.resize(0);


    Energies.push_back("Bending");
    Energy_constants.push_back(std::vector<double>{1.0,0.0});
    
    Energies.push_back("Bead");
    Energy_constants.push_back(params);

    Energies.push_back("Edge_reg");
    Energy_constants.push_back(std::vector<double>{1.0});

    

    Sim_handler = E_Handler(mesh, geometry, Energies, Energy_constants);
    Sim_handler.Add_Bead(&Bead_1);
    Sim_handler.boundary = false;

    M3DG.mesh = mesh;
    M3DG.geometry = geometry;
    M3DG.Add_bead(&Bead_1);
    Sim_handler.mesh = mesh;
    Sim_handler.geometry = geometry;
    M3DG.Sim_handler = &Sim_handler;

    std::ofstream Some_ofstream;
    std::vector<std::string> Constraints(0);
    Constraints.push_back("Volume");
    Eigen::VectorXd Lagrange_mults(1);
    Lagrange_mults(0) = 0.0;
    Sim_handler.Lagrange_mult = Lagrange_mults;
    Sim_handler.Trgt_vol = geometry->totalVolume();


    Sim_handler.Constraints = Constraints;
    double step = M3DG.integrate_Newton(Some_ofstream,0.0,Energies,false,Constraints,std::vector<std::string>{"Files"});
    std::cout<<"Took a step of " << step <<" wow";
    return 0; //Finishhh



    // OK so we want to do the same calculation now for the Hessian;

    const int Nverts = mesh->nVertices();

    // SparseMatrix<double> Hessian_bending = Sim_handler.H_Bending(Energy_constants[0]);

    // Eigen::MatrixXd Hessian_finite_diff(3*Nverts,3*Nverts);

    // VertexData<Vector3> Gradient_fwd(*mesh);
    // VertexData<Vector3> Gradient_prev(*mesh);

    // Gradient_fwd = Sim_handler.F_Bending(Energy_constants[0]);

    // // Gradient_fwd.array()

    // Eigen::VectorXd Difference(3*Nverts);


    // long int dim = 0;
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
    Constraints.resize(0);
    // Eigen::VectorXd Lagrange_mults(4);
    Lagrange_mults.resize(4);
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
    // std::ofstream Some_ofstream;

    
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
    std::cout<<" after moving backwards is " << geometry->totalVolume() << " \n";

    geometry->inputVertexPositions = geometry->inputVertexPositions + grad_vol*1e-3;

    VertexData<Vector3> grad_area = Sim_handler.F_SurfaceTension_2(std::vector<double>{-1.0});
    std::cout<<"The current area is " << geometry->totalArea() << " ";
    geometry->inputVertexPositions = geometry->inputVertexPositions + grad_area*1e-3;
    geometry->refreshQuantities();
    std::cout<<" after moving forward is " << geometry->totalArea() << " ";
    geometry->inputVertexPositions = geometry->inputVertexPositions - grad_area*2e-3;
    geometry->refreshQuantities();
    std::cout<<" after moving backwards is " << geometry->totalArea() << " \n";
    geometry->inputVertexPositions = geometry->inputVertexPositions + grad_area*1e-3;



    for(size_t step = 1; step < 20; step++){
    // M3DG.integrate(Some_ofstream,0.0,Energies,false);
    M3DG.integrate_Newton(Some_ofstream,0.0,Energies,false,Constraints,Constraints);
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



