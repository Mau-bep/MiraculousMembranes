// GEOMETRIC FLOW

// #include <stdlib.h>
#include <unistd.h>

#include <sys/resource.h>

// #include <omp.h>


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


std::vector<arcsim::Mesh> Saved_meshes(3);
std::vector<arcsim::Mesh> Saved_after_remesh(3);





// Polyscope visualization handle, to quickly add data to the surface

// polyscope::SurfaceMesh* psMesh;

// Some global variables
float TIMESTEP = -4;
float kappa = 1.0;


float H0 = 1.0;
float V_bar= (4/3)*PI*10; 



float nu;
float Interaction_str;
float c0;


float P0=100000.0;
float KA=1.0000;
float KB=0.0001;
float sigma=0.1;
double TS=0.0001;

double Curv_adap=1.0;
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

// RemeshBoundaryCondition defaultRemeshOptions;
    



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

void showSelected() {
    // pass
}

// void redraw() {
//     psMesh->updateVertexPositions(geometry->inputVertexPositions);
//     polyscope::requestRedraw();
// }

long get_mem_usage()
{
    struct rusage usage;
    int ret;
    ret = getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss; // in KB
}


bool measure_angles(){

    // This function will check if the angles get close to Pi

    double Pi = 3.1415926536;

    double angle;

    for(Vertex v : mesh->vertices()){
        // For every vertex, check all the angles

        for(Halfedge he : v.outgoingHalfedges()){
            // 
            angle = geometry->angle(he.corner());
            
            if(angle> Pi-Pi/20){
                
                std::cout<<"The angle on vertex "<< v.getIndex() << "is way too close to PI with a value of" << angle <<"\n";
                return true;
            }

        }

    }





    return false;
}

void area_ratios(){
    // 
    double max_ratio=-1;
    double min_ratio=100;
    double mean_min_ratio=0;
    double ratio;
    for(Edge e : mesh->edges()){
        Face f1 = e.halfedge().face();
        Face f2 = e.halfedge().twin().face();
        // Ok
        ratio = geometry->faceArea(f1)/geometry->faceArea(f2);
        ratio = min(ratio,1.0/ratio);
        mean_min_ratio += ratio;
        if(ratio<min_ratio) min_ratio = ratio;

    }
    std::cout<<"The average min ratio is"<< mean_min_ratio/mesh->nEdges()<< "and the min ratio is"<< min_ratio <<"\n";

}




void Save_mesh(std::string basic_name, size_t current_t) {
   // Build member variables: mesh, geometry
    Vector3 Pos;

  
    std::ofstream o(basic_name+"membrane_"+std::to_string(current_t)+".obj");
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
        if(mesh.verts[v]->adjf.size()<=2){
            std::cout<<"The number of neighbors is "<< mesh.verts[v]->adjf.size()<<"\n";
        
    
            
            // for(int f_index = 0 ; f_index < mesh.verts[v]->adjf.size();f_index++)
            // { 
            // std::cout<<"Deleting face \n";
            // mesh.remove(mesh.verts[v]->adjf[f_index]);    
            // }
            // std::cout<<"Now we delete the vert\n";
            // mesh.remove(mesh.verts[v]);
            
            // // mesh.remove(mesh.verts[v]->adjf[0]);

            std::cout<<"The vertex index to not consider are"<< v<<" out of a total of"<< mesh.verts.size()<<"\n";
            flag_warning=true;
            // flag=v;
            flags.push_back(v);
            continue;
        }
        
        
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

        int less_id1 = 0;
        int less_id2 = 0;
        int less_id3 = 0;

        for(size_t flag = 0 ; flag < flags.size(); flag++){
        if( id1 == flags[flag]|| id2 == flags[flag] || id3 == flags[flag]){
            non_manifold=true;
        }
        if(id1>flags[flag]&& flag_warning) less_id1+=1;
        if(id2>flags[flag]&& flag_warning) less_id2+=1;
        if(id3>flags[flag]&& flag_warning) less_id3+=1;
        }
        if(non_manifold){
            non_manifold=false;
            continue;
        }
        id1 = id1 - less_id1;
        id2 = id2 - less_id2;
        id3 = id3 - less_id3; 

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


    
    nu=1.0;
    Curv_adap=std::stod(argv[1]);
    // c0=std::stod(argv[2]);
    // KB=std::stod(argv[4]);
    Interaction_str=std::stod(argv[2]);
    int Init_cond = std::stoi(argv[3]);
    int Nsim = std::stoi(argv[4]);
    KA=std::stod(argv[5]);
    double radius=std::stod(argv[6]);
    KB = std::stod(argv[7]);
    Min_rel_length = 0.1;
    c0=0.0;

    bool continue_sim = false;

    bool pulling = false;
    bool arcsim = true;
    // I will do it so i can give this values
 
    arcsim::Cloth Cloth_1;
    Cloth_1.dump_info = false;
    arcsim::Cloth::Remeshing remeshing_params;

    int saved_mesh_idx = 0;
    bool Saving_last_states=true;

    std::vector<Vector3> Bead_pos_saved(3);



  
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();




    auto start_time_control= chrono::steady_clock::now();
    auto end_time_control= chrono::steady_clock::now();
    
    double remeshing_elapsed_time=0;
    double integrate_elapsed_time=0;
    double saving_mesh_time=0;
    



    TS=pow(10,-3);


    std::cout<< "Current path is " << argv[0];
    std::string filepath;
    if(Init_cond==1){
        filepath = "../../../input/sphere.obj";
    }
    if(Init_cond==2){
        std::cout<<"Icosphere condition\n";
        filepath = "../../../input/Icosphere.obj"; 
    }
    if(Init_cond==3){
        filepath = "../../../input/big_sphere.obj";
    }
    // std::string filepath = "../../../input/sphere_dense_40k.obj";
    // Load mesh
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    
    trgt_len = 0.0725;
    trgt_len = geometry->meanEdgeLength();

    

   


    // 
    // trgt_len=geometry->meanEdgeLength();
    V_bar=geometry->totalVolume();
    
    

    
      if(arcsim){
        std::cout<<"Settin remesher params";
        remeshing_params.aspect_min=0.2;
        remeshing_params.refine_angle=0.65;
        remeshing_params.refine_compression=1e-4;
        remeshing_params.refine_velocity=1.0;
        remeshing_params.size_max=trgt_len*3.0;
        remeshing_params.size_min=trgt_len*0.2;

        std::cout<<"Minimum edge length allowed is "<< trgt_len*0.2<<" muak\n";

    
    }
    arcsim::Mesh remesher_mesh = translate_to_arcsim(mesh,geometry);
    Cloth_1.mesh = remesher_mesh;
    Cloth_1.remeshing = remeshing_params; 
    arcsim::dynamic_remesh(Cloth_1);
    delete mesh;
    delete geometry;
    std::tie(mesh_uptr, geometry_uptr) = translate_to_geometry(Cloth_1.mesh);
    arcsim::delete_mesh(Cloth_1.mesh);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    

    ORIG_VPOS = geometry->inputVertexPositions;
    CoM = geometry->centerOfMass();

    // M3DG.Add_bead(&Bead_1);

    
    std::stringstream nustream;
    std::stringstream c0stream;
    std::stringstream KAstream;
    std::stringstream KBstream;

    std::stringstream radiusstream;
    std::stringstream Interactionstrstream;
    
    
    // std::stringstream H0stream;
    // std::stringstream kappastream;
    // std::stringstream sigmastream;
  

    


    nustream << std::fixed << std::setprecision(3) << nu;
    c0stream << std::fixed << std::setprecision(3) << c0;
    KAstream << std::fixed << std::setprecision(3) << KA;
    KBstream << std::fixed << std::setprecision(6) << KB;
    radiusstream << std::fixed << std::setprecision(3) << radius;
    Interactionstrstream << std::fixed << std::setprecision(6) << Interaction_str;


    
    

    std::string first_dir="../Results/Mem3DG_Bead_Reciprocal_arcsim_Phase/";
    int status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"If this name is 0 the directory was created succesfully "<< status ;

    std::string basic_name=first_dir+"nu_"+nustream.str()+"_radius_"+radiusstream.str()+"_KA_"+KAstream.str()+"_KB_"+KBstream.str()+"_strength_"+Interactionstrstream.str()+"_Init_cond_"+std::to_string(Init_cond)+"_Nsim_"+std::to_string(Nsim)+"/";
    status = mkdir(basic_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    std::cout<<"\nIf this number is 0 the directory was created succesfully "<< status<<"\n" ;

    std::string filename = basic_name+"Output_data.txt";

    std::string filename3 = basic_name+"Output_data.txt";

    std::ofstream Sim_data(filename3);
    Sim_data<<"T_Volume T_Area time Volume Area E_vol E_sur E_bend grad_norm backtrackstep\n";
    Sim_data.close();


    if(continue_sim){
        // If i want to restart is a bit tricky.
        std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(basic_name+"Final_state.obj");

        mesh = mesh_uptr.release();
        geometry = geometry_uptr.release();

        // I need to load the positions from a text file 

        std::ifstream inputFile(basic_name+"Saved_bead_info.txt");

        string line;

        std::getline(inputFile,line);
        stringstream ss(line);
        string pos;
        std::vector<double> coords(0);
        
        std::cout<<"coords are \t ";
        while(getline(ss,pos,' ')){
            coords.push_back(std::stod(pos));
            std::cout<<pos <<"\t";
        }
        std::cout<<"\n";
        // Now i have a line that is x y z 
        
        Bead_1 = Bead(mesh,geometry,Vector3({coords[0],coords[1],coords[2]}),radius,Interaction_str);
    

    }
    else{

        double x_furthest=0.0;
        for(Vertex v : mesh->vertices()){
            if(geometry->inputVertexPositions[v].x>x_furthest){
                x_furthest=geometry->inputVertexPositions[v].x;
            }

        }
    
        x_furthest+=radius*1.9;
        Bead_1 = Bead(mesh,geometry,Vector3({x_furthest,0.0,0.0}),radius,Interaction_str);
    
    
    }
    std::vector<Bead> Beads;


    

    Bead_1.interaction="Shifted_LJ_Normal_nopush";
    Bead_1.rc=radius*2.0;
    Beads.push_back(Bead_1);
    






    M3DG = Mem3DG(mesh,geometry);
    for( size_t i = 0 ; i < Beads.size() ; i++) M3DG.Add_bead(&Beads[i]);






    std::vector<std::string> Bead_filenames;
    std::ofstream Bead_datas;

    







    

    for(size_t i = 0 ; i < Beads.size(); i++){
        Bead_filenames.push_back(basic_name+ "Bead_"+std::to_string(i)+"_data.txt");
        Bead_datas = std::ofstream(Bead_filenames[i]);
        
        Bead_datas<<"####### This data is taken every 250 steps just like the mesh radius is " << radius<<" \n";
        Bead_datas.close();
    }


    std::string filename2 = basic_name + "Bead_data.txt";
    
    std::ofstream Bead_data(filename2);

    bool Save_bead_data=false;
    bool Save_output_data=false;
    bool small_Ts;
    Bead_data<<"####### This data is taken every 500 steps just like the mesh dump, radius is " << radius<<" \n";
    Bead_data.close();
    // Here i want to run my video
    size_t n_vert;
    size_t n_vert_old;
    size_t n_vert_new;
    double Volume;
    double Area;
    double nu_obs;
    double nu_evol;
    double nu_0;
    double c_null;
    size_t counter=0;
    double time=0.01;
    double dt_sim=0.0;
    int sys_time = 0;
    bool some_flag=false;

    bool seam = false;
    start = chrono::steady_clock::now();

    long currentMem = get_mem_usage();
    std::cout<<"THe code is using " << currentMem <<" in KB";
 
    bool saving_state = true;
    for(size_t current_t = 0; current_t <= 300000; current_t++ ){
        // for(size_t non_used_var=0;non_used_var<100;)
        // MemF.integrate(TS,sigma,kappa,H0,P,V0);
        
        
        if(arcsim){
            start_time_control=chrono::steady_clock::now();
            // Creo que aqui puedo guardar la mesh, 
            

            // if(measure_angles()) std::cout<<"The close to pi angle was due to integration\n";
            arcsim::Mesh remesher_mesh2 = translate_to_arcsim(mesh,geometry);
            Cloth_1.mesh=remesher_mesh2;
            if(Saving_last_states){
                saved_mesh_idx = (saved_mesh_idx+1)%3;
                arcsim::delete_mesh(Saved_meshes[saved_mesh_idx]);

                Bead_pos_saved[saved_mesh_idx] = Beads[0].Pos;
                Saved_meshes[saved_mesh_idx] = arcsim::deep_copy(remesher_mesh2);




            }
            Cloth_1.remeshing=remeshing_params;

            // std::cout<<"\t \t \t Current ts is "<<current_t<<"\n";
            // if(current_t>2300) std::cout<<"\t "+std::to_string(current_t) + "\t";
            arcsim::dynamic_remesh(Cloth_1);
            


            if(Saving_last_states){
                arcsim::delete_mesh(Saved_after_remesh[saved_mesh_idx]);
                Saved_after_remesh[saved_mesh_idx] = arcsim::deep_copy(Cloth_1.mesh);
            }

            
            // if( current_t>366){
            //     std::cout<<"Saving mesh\n";
            //     Save_mesh(basic_name,current_t); 
            //     arcsim::save_obj(Cloth_1.mesh, basic_name + "Debugging_after_"+ std::to_string(current_t)+".obj" );
            //     Cloth_1.dump_info = false;
            // }
            // if(seam ){
            //    std::cout<<"Seam\n";
            //     break;
            // }
            // Bead_1 = M3DG.Bead_1;
            small_Ts = M3DG.small_TS;
            sys_time = M3DG.system_time;

            delete mesh;
            delete geometry;
            // std::cout<<"translating back?\n";
            std::tie(mesh_uptr, geometry_uptr) = translate_to_geometry(Cloth_1.mesh);
            arcsim::delete_mesh(Cloth_1.mesh);

                // std::cout<<"defining pointers\n";
            mesh = mesh_uptr.release();
            geometry = geometry_uptr.release();
            // Bead_1 = Bead(mesh,geometry,Bead_1.Pos,radius,Interaction_str);
            // Bead_1.rc = 2.0;
            // if(measure_angles()) std::cout<<"The close to pi angle was due to remeshing\n";
            // area_ratios();

            
            
            
            
            end_time_control=chrono::steady_clock::now();
            remeshing_elapsed_time+=chrono::duration_cast<chrono::milliseconds>(end_time_control-start_time_control).count();
            M3DG= Mem3DG(mesh,geometry);
            for(size_t i = 0 ; i<Beads.size(); i++){
                Beads[i].Reasign_mesh(mesh,geometry);
                M3DG.Add_bead(&Beads[i]);
            }
            M3DG.small_TS = small_Ts;
            M3DG.system_time = sys_time;
        }
        else{

        start_time_control=chrono::steady_clock::now();
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
        MutationManager Mutation_manager(*mesh,*geometry);
        remesh(*mesh,*geometry,Mutation_manager,Options);
        n_vert_new=mesh->nVertices();
        counter=counter+1; 
        }
        end_time_control=chrono::steady_clock::now();
        remeshing_elapsed_time+=chrono::duration_cast<chrono::milliseconds>(end_time_control-start_time_control).count();
        
        }

        // psMesh->remove();
        
        // psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath("dodecahedra"), geometry->inputVertexPositions,
        //                                     mesh->getFaceVertexList(), polyscopePermutations(*mesh));

        // psMesh->setSurfaceColor({0.9607, 0.6627, 0.7215});
        // psMesh->setEdgeWidth(1.0);

        // bool should_save  = false;
        

        
        // if(counter == 100) saving_state = false;
        // if(current_t < 1000)  {
        //     counter=99;
        //     saving_state = true;
        
        // }
        // if(current_t%1000==0)
        // {
        //     saving_state = true;
        //     counter = 0;
        //  }
        
        if(current_t%100==0  ){
        // {
            // currentMem = get_mem_usage();
            // std::cout<<"THe code is using " << currentMem <<" in KB\n";
            // std::cout<<"WIth "<< mesh->nFaces()<<" faces\n";
            // Bead_data.close();
            // Sim_data.close();
            
            start_time_control=chrono::steady_clock::now();
            // if(current_t%100==0){

            Save_mesh(basic_name,current_t);
            // counter+=1;
            // should_save = false;
            // }
            end_time_control = chrono::steady_clock::now();
            saving_mesh_time+=chrono::duration_cast<chrono::milliseconds>(end_time_control-start_time_control).count();
            Save_bead_data=true;
            Bead_data = std::ofstream(filename2,std::ios_base::app);
            Save_output_data=true;
            Sim_data = std::ofstream(filename3,std::ios_base::app);
            

        
        }
        
        // if(current_t%100==0){

        //     Save_output_data=true;
        // }
        if(current_t%1000==0) {

            end=chrono::steady_clock::now();
            n_vert=mesh->nVertices();
            std::cout<< "THe number of vertices is "<< n_vert <<"\n";    
            std::cout<< "the avg edge lenghth is " <<geometry->meanEdgeLength()<<"\n";
            std::cout << "The avg edge length is = " << std::fixed << std::setprecision(10) << geometry->meanEdgeLength() << std::endl;

            Volume= geometry->totalVolume();
            Area=geometry->totalArea();
            nu_obs=3*Volume/(4*PI*pow( Area/(4*PI) ,1.5 ));
            H0=sqrt(4*PI/Area)*c0/2.0;


            // c_null=2*H0 *pow(Area/(4*PI),0.5);
            std::cout<< "The reduced volume is "<< nu_obs << "\n";
            if(current_t==0){
                nu_0=nu_obs;
            }

            std::cout<< "The spontaneous curvature is " << H0<< "\n";
            std::cout << "Current t is " << current_t <<"\n";
            std::cout << "The system time is " << time <<"\n\n";
            
            std::cout<<"A thousand iterations took "<<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" miliseconds\n";

            
            // polyscope::screenshot(basic_name+std::to_string(current_t)+".jpg",true);
   
            start = chrono::steady_clock::now();
            if(M3DG.small_TS){
                break;
            }
    

        }
        // std::cout<<"Redeclaring M3DG\n";
        // std::cout<<"Bead_1 position changed? "<< Bead_1.Pos << " \n";
        
        // std::cout<<"1\n";
        nu_evol= time<50 ? nu_0 + (nu-nu_0)*time/50 : nu; 
        // std::cout<<"2\n";
        
        start_time_control = chrono::steady_clock::now();
        // std::cout<<"3\n";
        // std::cout<<"Integrating\n";
        
        dt_sim=M3DG.integrate(TS,V_bar,nu_evol,c0,P0,KA,KB,sigma,Sim_data, time,Save_bead_data,Bead_filenames,Save_output_data,pulling);
        Bead_data.close();
        Sim_data.close();
        // Bead_data = std::ofstream(filename2,std::ios_base::app);
        // Sim_data = std::ofstream(filename3,std::ios_base::app);
        // std::cout<<"4\n";
        
        end_time_control = chrono::steady_clock::now();
        // std::cout<<"5\n";
        
        integrate_elapsed_time += chrono::duration_cast<chrono::milliseconds>(end_time_control-start_time_control).count(); 
        Save_output_data=false;
        Save_bead_data=false;


        if(current_t%1000==0){
            std::cout<< "Remeshing has taken a total of "<< remeshing_elapsed_time <<" milliseconds\n" << "Saving the mesh has taken a total of "<< saving_mesh_time<< "milliseconds \n Integrating the forces has taken a total of "<< integrate_elapsed_time <<" milliseconds \n\n"; 
        }

        if(dt_sim==-1){
            std::cout<<"Sim broke or timestep very small\n";
            break;
        }
        else{
            // std::cout<<"Adding time\n";
            time+=dt_sim;
            // std::cout<<"SUccesfully\n";
        }





    }
    Sim_data.close();
    Bead_data.close();

    Vector3 Pos;
    std::cout<<"Wrapping up\n";


    if(continue_sim){
    std::ofstream beads_saved(basic_name+"Saved_bead_info_cont.txt");
    beads_saved << std::setprecision(std::numeric_limits<double>::max_digits10);


    for( int k = 0 ; k < 3 ; k++){
        std::cout<<Bead_pos_saved[k].x << " " << Bead_pos_saved[k].y << " " << Bead_pos_saved[k].z <<" \n";
        std::cout<<"Saved mesh idx is " << saved_mesh_idx<<" \n";
        arcsim::save_obj(Saved_meshes[saved_mesh_idx],basic_name+"Saved_continued_final_frame_"+std::to_string(3-k)+".obj");
        

        beads_saved<< Bead_pos_saved[saved_mesh_idx].x<<" "<<Bead_pos_saved[saved_mesh_idx].y<<" "<<Bead_pos_saved[saved_mesh_idx].z<<" \n";

        
        saved_mesh_idx = (saved_mesh_idx -1 +3)%3;

    }
    }
    else{
     std::ofstream beads_saved(basic_name+"Saved_bead_info.txt");
    beads_saved << std::setprecision(std::numeric_limits<double>::max_digits10);


    for( int k = 0 ; k < 3 ; k++){
        std::cout<<Bead_pos_saved[k].x << " " << Bead_pos_saved[k].y << " " << Bead_pos_saved[k].z <<" \n";
        std::cout<<"Saved mesh idx is " << saved_mesh_idx<<" \n";
        arcsim::save_obj(Saved_meshes[saved_mesh_idx],basic_name+"Saved_final_frame_"+std::to_string(5-2*k)+".obj");
        arcsim::save_obj(Saved_after_remesh[saved_mesh_idx],basic_name+"Saved_final_frame_"+std::to_string(6-2*k)+".obj");

        beads_saved<< Bead_pos_saved[saved_mesh_idx].x<<" "<<Bead_pos_saved[saved_mesh_idx].y<<" "<<Bead_pos_saved[saved_mesh_idx].z<<" \n";

        
        saved_mesh_idx = (saved_mesh_idx -1 +3)%3;

    }   
    }


    std::ofstream o(basic_name+"Final_state.obj");
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
    



    delete mesh;
    delete geometry;

    return EXIT_SUCCESS;
}




