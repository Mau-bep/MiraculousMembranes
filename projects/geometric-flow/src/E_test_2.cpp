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
  
    // std::cout<<"Adding vertices?\n";
    for (int v = 0; v < mesh->nVertices(); v++) {
            // const Vert *vert0 = mesh0.verts[v];
            Vector3 pos_orig= geometry->inputVertexPositions[v];
            arcsim::Vec3 pos;
            pos[0]=pos_orig.x;
            pos[1]=pos_orig.y;
            pos[2]=pos_orig.z;
            arcsim::Vert *vert1 = new arcsim::Vert(pos, 1);
            mesh1.add(vert1);
            
    }

    // std::cout<<"Adding nodes?\n";
    for (int v = 0; v < mesh->nVertices(); v++) {
            // const Vert *vert0 = mesh0.verts[v];
            Vector3 pos_orig= geometry->inputVertexPositions[v];
            arcsim::Vec3 pos;
            pos[0]=pos_orig.x;
            pos[1]=pos_orig.y;
            pos[2]=pos_orig.z;
            // arcsim::Vert *vert1 = new arcsim::Vert(pos, 0, 1);
            // mesh1.add(vert1);
            mesh1.add(new arcsim::Node(pos,arcsim::Vec3(0)));
            mesh1.verts[v]->node=mesh1.nodes[v];
            mesh1.nodes[v]->verts.push_back(mesh1.verts[v]);
    }


    // std::cout<<"Adding edges?\n";
    for (int e = 0; e<mesh->nEdges();e++){
        Edge e_orig= mesh->edge(e);

        arcsim::Edge *edge1 = new arcsim::Edge(mesh1.nodes[e_orig.firstVertex().getIndex()],mesh1.nodes[e_orig.secondVertex().getIndex()],0);
        mesh1.add(edge1);

    }
    
    



    // std::cout<<"Adding faces?\n";
    for(int f=0; f< mesh->nFaces(); f++){
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
    
    for(int v = 0 ; v<mesh.verts.size(); v++){
        arcsim::Vec3 pos_old = mesh.nodes[v]->x;
        v_pos.x=pos_old[0];
        v_pos.y=pos_old[1];
        v_pos.z=pos_old[2];
        simpleMesh.vertexCoordinates.push_back(v_pos);
    }
    for(int f = 0 ; f<mesh.faces.size();f++){
        
        std::vector<size_t> polygon(3);
        polygon[0] = mesh.faces[f]->v[0]->index;
        polygon[1] = mesh.faces[f]->v[1]->index;
        polygon[2] = mesh.faces[f]->v[2]->index;
        simpleMesh.polygons.push_back(polygon);
    }
    // std::cout<<"Does this happen after loading the data to create the mesh?\n";
    // std::cout<<" THe information in the mesh is, "<< simpleMesh.vertexCoordinates.size()<<"number of vertices\n";
  auto lvals = makeManifoldSurfaceMeshAndGeometry(simpleMesh.polygons, simpleMesh.vertexCoordinates);



 return std::tuple<std::unique_ptr<ManifoldSurfaceMesh>,
                    std::unique_ptr<VertexPositionGeometry>>(std::move(std::get<0>(lvals)),  // mesh
                                                             std::move(std::get<1>(lvals))); // geometry

}



int main(int argc, char** argv) {

    
    nu=std::stod(argv[1]);
    c0=std::stod(argv[2]);
    KA=std::stod(argv[3]);
    KB=std::stod(argv[4]);
    bool evolve=false;
    size_t target_index=120;
    // I will do it so i can give this values
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    
    bool with_bead=false;

    bool Save_data=false;
    TS=5*pow(10,-3);
    double time;
    double dt_sim;
    bool test_remesher = true;


    bool preserving_vol=true;


    std::cout<< "Current path is " << argv[0]<<"\n";

    std::string filepath = "../../../input/sphere.obj";


    if(preserving_vol){
        filepath = "../../../input/bunny.obj";
    }
    // std::string filepath = "../../../input/8_octahedron.obj";
    // std::string filepath = "../../../input/Simple_cil_regular.obj";
    
    // std::string filepath = "../../../input/bloodcell.obj";
    // std::string filepath = "../input/sphere.obj"; //this is for debug
    
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    std::cout<<"The mesh is correctly loaded\n";
    



    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    
    trgt_len=geometry->meanEdgeLength();
    V_bar=geometry->totalVolume();
    
    EdgeData<int> No_remesh(*mesh,0);
    No_remesh[0]=1;




    ORIG_VPOS = geometry->inputVertexPositions;
    CoM = geometry->centerOfMass();
    double radius=1.0;
    double Interaction_str=1.0;
    Bead_1 = Bead(mesh,geometry,Vector3({5.8,0.0,0.0}),radius,Interaction_str);
    Bead_1.interaction="pulling";
    M3DG = Mem3DG(mesh,geometry,Bead_1);



    // Add visualization options.
    // psMesh->setSmoothShade(false);
    
    // psMesh->setSurfaceColor({0.9607, 0.6627, 0.7215});// not orange
    
    std::stringstream nustream;
    std::stringstream c0stream;
    std::stringstream KAstream;
    std::stringstream KBstream;
    
    
    std::stringstream Curv_adapstream;
    std::stringstream Min_rel_lengthstream;


    nustream << std::fixed << std::setprecision(3) << nu;
    c0stream << std::fixed << std::setprecision(3) << c0;
    KAstream << std::fixed << std::setprecision(3) << KA;
    KBstream << std::fixed << std::setprecision(6) << KB;



    Curv_adapstream << std::fixed << std::setprecision(2) << Curv_adap;
    Min_rel_lengthstream << std::fixed << std::setprecision(2) <<Min_rel_length;
    std::string basic_name;
    std::string first_dir;



    double avg_edge= geometry->meanEdgeLength();


    if(test_remesher){
        first_dir = "../Results/Tests_remesher/";
        int status1 = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        basic_name ="../Results/Tests_remesher/nu_"+nustream.str()+"_c0_"+c0stream.str()+"_KA_"+KAstream.str()+"_KB_"+KBstream.str()+"/";
        int status = mkdir(basic_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        std::cout<<"The status is "<< status <<" \n";
        std::cout<<"Testing remesher \n";
        std::cout << "The avg edge length is = " << std::fixed << std::setprecision(10) << geometry->meanEdgeLength() << std::endl;

        arcsim::Mesh remesher_mesh;
        remesher_mesh = translate_to_arcsim(mesh,geometry);


        std::tie(mesh_uptr, geometry_uptr) = translate_to_geometry(remesher_mesh);
        mesh = mesh_uptr.release();
        geometry = geometry_uptr.release();

    

        arcsim::Cloth Cloth_1;
        Cloth_1.mesh=remesher_mesh;
        arcsim::Cloth::Remeshing remeshing_params;
        remeshing_params.aspect_min=0.2;
        remeshing_params.refine_angle=0.1;
        remeshing_params.refine_compression=0.05;
        remeshing_params.refine_velocity=1.0;
        remeshing_params.size_max=avg_edge*2.0;
        remeshing_params.size_min=avg_edge*0.1;

        Cloth_1.remeshing=remeshing_params;
   
        for(int i =0 ; i<Cloth_1.mesh.nodes.size(); i++){
            if(arcsim::is_seam_or_boundary( Cloth_1.mesh.nodes[i])) std::cout<<"THis node is a seam or boundary? wtf\n\n";
        }


        // arcsim::dynamic_remesh(Cloth_1);


        // arcsim::save_obj(Cloth_1.mesh, basic_name+"dynamicallyremeshed_girlie.obj");
        // std::cout<<"Saved mesh at "<< basic_name<<"dynamicallyremeshed_girlie.obj\n";
        // std::cout<<"After remeshing the mesh looks like this: \n , there are " << Cloth_1.mesh.verts.size()<<" vertices and "<< Cloth_1.mesh.edges.size()<< "edges\n";
        
        // delete mesh;
        // delete geometry;
        
        // std::tie(mesh_uptr, geometry_uptr) = translate_to_geometry(Cloth_1.mesh);
        // mesh = mesh_uptr.release();
        // geometry = geometry_uptr.release();


        auto start = chrono::steady_clock::now();
        auto end = chrono::steady_clock::now();
        double avg_time;
        double counter_iter;



    // Ok now is then the fun begins

        // I want to do normal flow.
        bool arcsim_remesh;
        if(c0>0){
        arcsim_remesh=true;
        }
        else{
        arcsim_remesh=false;
        }
        bool preserving_vol=true;
        double targ_vol=geometry->totalVolume();
        std::ofstream Sim_data;
        if(preserving_vol){
            std::string filename = basic_name+"Output_data.txt";
            
            Sim_data=std::ofstream(filename);
            bool Save=false;
             Sim_data<<"T_Volume time Volume Area E_vol E_sur grad_norm backtrackstep\n";
        }

        int n_vert_old;
        int n_vert_new;
        int counter;
        for(size_t current_t = 0 ; current_t<50000; current_t ++){
            // std::cout<<"Current t is "<< current_t <<"\n";
            // fIRST THING IS REMESHIng 
            if(arcsim_remesh){
                start = chrono::steady_clock::now();
                // std::cout<<"translating?\n";
                arcsim::Mesh remesher_mesh2 = translate_to_arcsim(mesh,geometry);
                Cloth_1.mesh=remesher_mesh2;
                Cloth_1.remeshing=remeshing_params;
                // std::cout<<"remeshing\n";
                arcsim::dynamic_remesh(Cloth_1);
                delete mesh;
                delete geometry;
                // std::cout<<"translating back?\n";
                std::tie(mesh_uptr, geometry_uptr) = translate_to_geometry(Cloth_1.mesh);
                arcsim::delete_mesh(Cloth_1.mesh);

                // std::cout<<"defining pointers\n";
                mesh = mesh_uptr.release();
                geometry = geometry_uptr.release();
                end=chrono::steady_clock::now();

                if(current_t==0){
                avg_time=chrono::duration_cast<chrono::milliseconds>(end-start).count();
                }
                else{
                avg_time = avg_time + (chrono::duration_cast<chrono::milliseconds>(end-start).count()- avg_time )/(current_t+1);
                }
                // std::cout<<"reinitializing M3DG?\n";
                // M3DG = Mem3DG(mesh,geometry);
            }
            else{
                start = chrono::steady_clock::now();
                n_vert_old=mesh->nVertices();
                n_vert_new=1;

            counter=0;
            while(n_vert_new!=n_vert_old && counter<10){
                // std::cout<<"Remeshing\n";
                n_vert_old=n_vert_new;
                RemeshOptions Options;
                Options.targetEdgeLength=trgt_len;
                Options.curvatureAdaptation=Curv_adap;
                Options.maxIterations=10;
                Options.minRelativeLength=0.1;
                Options.smoothStyle=RemeshSmoothStyle::Circumcentric;
                Options.boundaryCondition=RemeshBoundaryCondition::Tangential;
                // Options.remesh_list=true;
                // Options.No_remesh_list=No_remesh;
                MutationManager Mutation_manager(*mesh,*geometry);
                remesh(*mesh,*geometry,Mutation_manager,Options);
                n_vert_new=mesh->nVertices();
                counter=counter+1; 
            }
            end=chrono::steady_clock::now();
            if(current_t==0){
                avg_time=chrono::duration_cast<chrono::milliseconds>(end-start).count();
                }
                else{
                avg_time = avg_time + (chrono::duration_cast<chrono::milliseconds>(end-start).count()- avg_time )/(current_t+1);
                }

            }


            // I am doing it wrong, the idea is to move them depending on the position, not the normals.
            if(preserving_vol){ 
                bool Save=false;
                if(current_t%250==0){
                    Save=true;
                }
                double Volume;
                M3DG=Mem3DG(mesh,geometry);
                Volume = geometry->totalVolume();
                dt_sim = M3DG.integrate(TS,targ_vol,P0,KA,Sim_data,time,Save);
                time+=dt_sim;

            }
            else{
            
            Vector3 rhat;
            
            
            for(int i =0 ; i<mesh->nVertices(); i++){
                rhat = geometry->inputVertexPositions[i].unit();
                geometry->inputVertexPositions[i]-=1e-2*rhat;

            }
            geometry->refreshQuantities();
            }

            // M3DG = Mem3DG(mesh,geometry);
            // std::cout<<"Is the problem after this?";
            
            // VertexData<Vector3> Volume_grad = M3DG.OsmoticPressure(1.0);
            // std::cout<<"What is killing?\n";
            // geometry->inputVertexPositions+=1e-2*M3DG.OsmoticPressure(1.0);
            // std::cout<<"refreshing maybe?\n";
            geometry->refreshQuantities();
            
            if(current_t%100==0){
                // std::cout<<"saving?\n";
                std::cout<<"The average time for the "; 
                if(arcsim_remesh) std::cout<<"arcsim "; 
                else std::cout<<" normal ";  
                std::cout<<"remesher is"<< avg_time<< " milliseconds\n";
            Save_mesh(basic_name,arcsim_remesh,current_t);    
            //  Save_mesh(basic_name,current_t);    
           
            }


        }


    //     for(size_t current_t=0; current_t <1000;current_t++){
   
    // if(current_t%200==0){
    //     // filename = basic_name+"Vol_Gradient_evaluation_"+std::to_string(current_t) + ".txt";
    //     // std::ofstream Gradient_data(filename);
    //     // Gradient_data.open(filename);
    //     // std::ofstream o(basic_name+std::to_string(current_t)+".obj");
    //     Gradient_data<< "Volume grad\n";
        
        

    //     Gradient_data.close();
    //     // double A_bar=4*PI*pow(3*V_bar/(4*PI*nu_evol),2.0/3.0);



    //     Volume= geometry->totalVolume();
    //     Area=geometry->totalArea();
    //     nu_obs=3*Volume/(4*PI*pow( Area/(4*PI) ,1.5 ));
    //     H0=sqrt(4*PI/Area)*c0/2.0;
        
    //     if(current_t==0){
    //         nu_0=nu_obs;
    //     }
   

    // }
     
  
  
    // Sim_data.close();

 
    
    
    // }













        std::cout<<"Done testing remeshing\n";



    }





    if(with_bead){

    first_dir="../Results/Tests_bead/";
    int status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"If this name is 0 the directory was created succesfully "<< status ;

    // first_dir="./Tests/";
    // status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"\nIf this name is 0 the directory was created succesfully "<< status ;
    
    basic_name ="../Results/Tests_bead/nu_"+nustream.str()+"_c0_"+c0stream.str()+"_KA_"+KAstream.str()+"_KB_"+KBstream.str()+"/";
    status = mkdir(basic_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"\nIf this number is 0 the directory was created succesfully "<< status<<"\n" ;
    }
    else{

    first_dir="../Results/Tests_general/";
    int status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"If this name is 0 the directory was created succesfully "<< status ;

    // first_dir="./Tests/";
    // status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"\nIf this name is 0 the directory was created succesfully "<< status ;
    
    basic_name ="../Results/Tests_general/nu_"+nustream.str()+"_c0_"+c0stream.str()+"_KA_"+KAstream.str()+"_KB_"+KBstream.str()+"/";
    status = mkdir(basic_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        
    }

    std::string filename_basic = basic_name+"Output_data.txt";

    std::ofstream Sim_data;
    

    std::string filename2 = basic_name + "Bead_data.txt";
    std::ofstream Bead_data(filename2);
    bool Save_bead_data=false;
        
    if(with_bead){
        
        Bead_data<<"####### This data is taken every 250 steps just like the mesh radius is " << radius<<" \n";
    
    }

    // Here i want to run my video
    size_t n_vert;
    size_t n_vert_old;
    size_t n_vert_new;
    double Volume;
    double Area;
    double nu_obs;
    double nu_evol;
    double nu_0;

    size_t counter=0;
    time=0.01;
    dt_sim=0.0;

    bool pulling;
    // Here i have loaded the mesh and i am ready to start testing


    // I may need to run this and save things so i can se the evolution...
    std::ofstream Gradient_data;
    std::ofstream Gradient_data_area;
    std::ofstream Gradient_data_bending;
    std::ofstream Gradient_data_bending_norms;
    std::ofstream Gradient_data_tot_area;

    std::string filename;

    if(with_bead){
    std::ofstream Gradient_data_bead;

    std::ofstream Gradient_data_bead_dx;
    


    filename = basic_name+"Bead_Gradient_evaluation_dx_"+std::to_string(0) + ".txt";
    Gradient_data_bead_dx.open(filename);
    Gradient_data_bead_dx << "Bead_gradient dx evaluation\n";
    M3DG.Grad_Bead_dx(Gradient_data_bead_dx,true);

    Gradient_data_bead_dx.close();

    }


    // Since this is the test secion, i want to check something
    double Total_A=0;
    double Total_A_dual_bar=0;
    double Total_A_dual_circ=0;
    Total_A=geometry->totalArea();
    double max_x=0;
    double max_y=0;
    // double 
    for( Vertex v : mesh->vertices()){
        if(max_x<geometry->inputVertexPositions[v].x){
            max_x=geometry->inputVertexPositions[v].x;
        }    
        Total_A_dual_bar+=geometry->barycentricDualArea(v);
        Total_A_dual_circ+=geometry->circumcentricDualArea(v);
    }

    std::cout<<"THe areas a are: Total "<<Total_A <<" Barycentric "<< Total_A_dual_bar <<" and Circumcentric "<< Total_A_dual_circ <<" \n";
    
    Volume=geometry->totalVolume();
    std::cout<<"The radius (according to area) is " << sqrt(Total_A/(4*PI))<<"\n";
    std::cout<<"THe radius (according to volume) is" << pow(  3*Volume/(4*PI)   ,1.0/3.0)<<"\n"; 
    std::cout<<"THe maximum x coordinate is "<< max_x <<"\n";


    // Lets check what is the dihedral angle in a sphere, to get an idea.

    int count=0;
    double dihedral=0.0;
    double angle;
    double maxdihedral=0.0;
    double mindihedral=1.0;
    for (Halfedge he : mesh->halfedges()){
        count+=1;
        angle=abs(geometry->dihedralAngle(he));
        dihedral+=angle;
        if(angle>maxdihedral){
            maxdihedral=angle;
        }
        if(angle<mindihedral){
            mindihedral=angle;
        }
        // std::cout<<"Angle is"<<angle <<"\n";
        
    }

    std::cout<<"\n\n The mean dihedral angle ( in radians ) is" << dihedral/count<<"\n";
    std::cout<<"The mean dihedral angle ( in degrees ) is" << 360*dihedral/(count*2*3.141592)<<"\n";
    std::cout<<"The minimum dihedral angle (in radians) is " << mindihedral <<" and the maximum is "<< maxdihedral <<"\n\n";

    // Lets check the normal difference

    count=0;
    double avg_diff=0;
    double maxdiff=0;
    double mindiff=100;
    for(Edge e : mesh->edges()){
        Halfedge he = e.halfedge();
        Vertex v1 = he.vertex();
        Vertex v2 = he.next().vertex();

        double norm2 = (geometry->vertexNormalAngleWeighted(v1)-geometry->vertexNormalAngleWeighted(v2)).norm2();
        avg_diff+=norm2;
        count+=1;
        if(norm2>=maxdiff){
            maxdiff=norm2;
        }
        if(norm2<mindiff){
            mindiff=norm2;
        }

    }
    
    std::cout<<"\n\n The mean difference in normals is "<< avg_diff/count<<" \n";
    std::cout<<" The minimum normal difference is " << mindiff <<" \n";
    std::cout<< " The  maxmimum difference in normals is"<< maxdiff << " \n\n\n";
    
    std::cout<<"This normalized reads like this \n";
    std::cout<<"\n\n The mean difference in normals is "<< avg_diff/(count*maxdiff)<<" \n";
    std::cout<<" The minimum normal difference is " << mindiff/maxdiff <<" \n";
    std::cout<< " The  maxmimum difference in normals is"<< maxdiff/maxdiff << " \n\n\n";
    




    if(evolve){
    for(size_t current_t=0; current_t <1000;current_t++){
    if(current_t==0){
    Sim_data.open(filename_basic);
    }
    else{
    Sim_data.open(filename_basic,std::ios_base::app);
    }
   
    if(false){
        // if(current_t%10==0 ){
        n_vert_old=mesh->nVertices();
        n_vert_new=1;

        counter=0;
        while(n_vert_new!=n_vert_old && counter<10){
        // std::cout<<"Remeshing\n";
        n_vert_old=n_vert_new;
        RemeshOptions Options;
        Options.targetEdgeLength=trgt_len;
        Options.curvatureAdaptation=Curv_adap;
        Options.maxIterations=10;
        Options.minRelativeLength=Min_rel_length;
        Options.smoothStyle=RemeshSmoothStyle::Circumcentric;
        Options.boundaryCondition=RemeshBoundaryCondition::Tangential;
        // Options.remesh_list=true;
        // Options.No_remesh_list=No_remesh;
        MutationManager Mutation_manager(*mesh,*geometry);
        remesh(*mesh,*geometry,Mutation_manager,Options);
        n_vert_new=mesh->nVertices();
        counter=counter+1; 
        }

        dihedral=0.0;
        mindihedral=1.0;
        maxdihedral=0.0;
        for (Halfedge he : mesh->halfedges()){
        count+=1;
        angle=abs(geometry->dihedralAngle(he));
        dihedral+=angle;
        if(angle>maxdihedral){
            maxdihedral=angle;
        }
        if(angle<mindihedral){
            mindihedral=angle;
        }
        // std::cout<<"Angle is"<<angle <<"\n";
        
    }

        std::cout<<"The mean dihedral angle ( in radians ) is" << dihedral/count<<"\n";
        std::cout<<"The mean dihedral angle ( in degrees ) is" << 360*dihedral/(count*2*3.141592)<<"\n";
        std::cout<<"The minimum dihedral angle (in radians) is " << mindihedral <<" and the maximum is "<< maxdihedral <<"\n";



        }

        n_vert=mesh->nVertices();
        std::cout<< "THe number of vertices is "<< n_vert <<"\n";    
        std::cout << "The avg edge length is = " << std::fixed << std::setprecision(10) << geometry->meanEdgeLength() << std::endl;

        Volume= geometry->totalVolume();
        Area=geometry->totalArea();
        nu_obs=3*Volume/(4*PI*pow( Area/(4*PI) ,1.5 ));
        // H0=sqrt(4*PI/Area)*c0/2.0;
        std::cout<< "The volume is "<< Volume << "\n";

        std::cout<< "The reduced volume is "<< nu_obs << "\n";
        if(current_t==0){
            nu_0=nu_obs;
        }

    if(current_t%1==0){

            Save_mesh(basic_name,current_t);
            if(with_bead)
            {
                Save_bead_data=false;
            }
        }
    
    if(current_t%200==0){
        filename = basic_name+"Vol_Gradient_evaluation_"+std::to_string(current_t) + ".txt";
        // std::ofstream Gradient_data(filename);
        Gradient_data.open(filename);
        // std::ofstream o(basic_name+std::to_string(current_t)+".obj");
        Gradient_data<< "Volume grad\n";
        M3DG.Grad_Vol(Gradient_data,P0,V_bar,true);

        Gradient_data.close();
        double A_bar=4*PI*pow(3*V_bar/(4*PI*nu_evol),2.0/3.0);

        filename = basic_name+"Area_Gradient_evaluation_"+std::to_string(current_t) + ".txt";
        Gradient_data_area.open(filename);
        Gradient_data_area<< "Area grad\n";
        M3DG.Grad_Area(Gradient_data_area,A_bar,KA,true);

        Gradient_data_area.close();

        filename = basic_name+"Bending_Gradient_evaluation_"+std::to_string(current_t) + ".txt";
        Gradient_data_bending.open(filename);
        Gradient_data_bending<< "Bending grad\n";
        double H_bar=sqrt(4*PI/A_bar)*c0/2.0;
        M3DG.Grad_Bending(Gradient_data_bending,H_bar,KB,true);
        Gradient_data_bending.close();


        filename = basic_name+"Bending_Gradient_evaluation_2_"+std::to_string(current_t) + ".txt";
        Gradient_data_bending_norms.open(filename);
        Gradient_data_bending_norms<< "Bending grad norm diff\n";
        // double H_bar=sqrt(4*PI/A_bar)*c0/2.0;
        M3DG.Grad_Bending_2(Gradient_data_bending_norms,H_bar,KB);
        Gradient_data_bending_norms.close();
        


        filename = basic_name+"Area_tot_Gradient_evaluation_"+std::to_string(current_t) + ".txt";
        Gradient_data_tot_area.open(filename);
        Gradient_data_tot_area<<"Area tot grad\n";
        M3DG.Grad_tot_Area(Gradient_data_tot_area,true);
        Gradient_data_tot_area.close();


        // filename = basic_name+"Bead_Gradient_evaluation_"+std::to_string(current_t) + ".txt";
        // Gradient_data_bead.open(filename);
        // Gradient_data_bead<< "Bead grad\n";
        // // double H_bar=sqrt(4*PI/A_bar)*c0/2.0;
        // M3DG.Grad_Bead(Gradient_data_bead,true,false);
        // Gradient_data_bead.close();


        Volume= geometry->totalVolume();
        Area=geometry->totalArea();
        nu_obs=3*Volume/(4*PI*pow( Area/(4*PI) ,1.5 ));
        H0=sqrt(4*PI/Area)*c0/2.0;
        
        if(current_t==0){
            nu_0=nu_obs;
        }
   

    }
     
    nu_evol= time<50 ? nu_0 + (nu-nu_0)*time/50 : nu; 
    if(with_bead){
        dt_sim=M3DG.integrate(TS,V_bar,nu_evol,c0,P0,KA,KB,Kd,Sim_data,time,Save_bead_data,Bead_data,Save_data,pulling);
        
    }
    else{
    dt_sim=M3DG.integrate(TS,V_bar,nu_evol,c0,P0,KA,KB,Kd,Sim_data,time,Save_data);
    }

    if(dt_sim==-1)
        {
        std::cout<<"Sim broke or timestep very small\n";
        break;
        }
    else{
        time+=dt_sim;
            
        }

    

    Sim_data.close();

    // filename = basic_name+"Bending_Gradient_analitical_terms.txt";
    // std::ofstream Gradient_data_bending_2(filename);
    // Gradient_data_bending_2<< "Bending gradients per component\n";
    // // double H_bar=sqrt(4*PI/A_bar)*c0/2.0;
    // M3DG.Bending_test(Gradient_data_bending_2,H_bar,KB);

    // Gradient_data_bending_2.close();

    
    
    
    }
}




    delete mesh;
    delete geometry;

    return EXIT_SUCCESS;
}



