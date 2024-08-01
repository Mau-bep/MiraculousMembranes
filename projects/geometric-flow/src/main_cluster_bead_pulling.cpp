// GEOMETRIC FLOW

// #include <stdlib.h>
#include <unistd.h>


// #include <omp.h>


#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <string>

#include <limits>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/surface/remeshing.h"



#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"

#include "happly.h"

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
float Interaction_str;
float c0;


float P0=100000.0;
float KA=1.0000;
float KB=0.0001;
float Kd=0.0;
double TS=0.0001;

double Curv_adap=0.1;
double Min_rel_length=0.05;
double trgt_len;
double prev_force;
bool Save_ss;
bool Stop_increasing;
double prev_strength;
double prev_E;
Vector3 Mem_Force;

Vector3 Bead_force;
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

void Save_mesh(std::string basic_name,std::string nametag , size_t current_t) {
   // Build member variables: mesh, geometry
    Vector3 Pos;
    std::ofstream o(basic_name+nametag+"_"+std::to_string(current_t)+".obj");
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
            mesh1.add(new arcsim::Node(pos,
                                 pos,
                                 pos,
                                 0,false));
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
    for(int v = 0 ; v<mesh.verts.size(); v++){
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

            std::cout<<"The vertex index to not consider are"<< v<<"\n";
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
    for (int flag = 0 ; flag< flags.size(); flag++){
        std::cout<< flags[flag]<<"\t ";
    }
    }
    // std::cout<<"hihi\n";
    bool non_manifold = false;
    for(int f = 0 ; f<mesh.faces.size();f++){
        
        std::vector<size_t> polygon(3);
        
        id1 = mesh.faces[f]->v[0]->index;
        id2 = mesh.faces[f]->v[1]->index;
        id3 = mesh.faces[f]->v[2]->index;

        int less_id1 = 0;
        int less_id2 = 0;
        int less_id3 = 0;

        for(int flag =0 ; flag < flags.size(); flag++){
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
    // std::cout<<"Done with translating";
    // std::cout<<"Does this happen after loading the data to create the mesh?\n";
    // std::cout<<" THe information in the mesh is, "<< simpleMesh.vertexCoordinates.size()<<"number of vertices\n";
  auto lvals = makeManifoldSurfaceMeshAndGeometry(simpleMesh.polygons, simpleMesh.vertexCoordinates);



 return std::tuple<std::unique_ptr<ManifoldSurfaceMesh>,
                    std::unique_ptr<VertexPositionGeometry>>(std::move(std::get<0>(lvals)),  // mesh
                                                             std::move(std::get<1>(lvals))); // geometry
}






// void init_physics(const string &json_file, bool is_reloading, arcsim::Simulation &sim){
        
    // vector<arcsim::Mesh *> base_meshes(sim.obstacles.size());
    // arcsim::load_json(json_file,sim);
// }

int main(int argc, char** argv) {

    
    double rc=std::stod(argv[1]);
    nu=1.0;
    double pulling_force=0.0;
    // c0=std::stod(argv[2]);
    // KA=std::stod(argv[3]);
    // KB=std::stod(argv[4]);
    Interaction_str=std::stod(argv[2]);
    int Init_cond = std::stoi(argv[3]);
    int Nsim = std::stoi(argv[4]);
    bool pulling =true;
    size_t SS_index=0;
    c0=0.0;
    KB=std::stod(argv[5]);;
    KA = std::stod(argv[6]);
    // KA= 100000;
    // KB=0.01;


    

    // I will do it so i can give this values
 
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    

    auto start_saving = chrono::steady_clock::now();
    auto end_saving = chrono::steady_clock::now();


    auto start_report = chrono::steady_clock::now();
    auto end_report = chrono::steady_clock::now();


    TS=pow(10,-4);


    std::cout<< "Current path is " << argv[0];
    std::string filepath;
    if(Init_cond==1){
        filepath = "../../../input/sphere.obj";
    }
    if(Init_cond==2){
        filepath = "../../../input/sphere_dense_40k.obj"; 
    }
    if(Init_cond==3){
        filepath = "../../../input/big_sphere.obj";
    }
    if(Init_cond==4){
        filepath = "../Results/Mem3DG_Bead_Pulling_rc_calib/nu_1.000_rc_2.000_KA_500.000_KB_1.000000_strength_30.000000_Init_cond_1_Nsim_667/membrane_28000.obj";
    }
    // std::string filepath = "../../../input/sphere_dense_40k.obj";
    // Load mesh
    std::cout<<"The mesh is"<< filepath<<"\n";
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    
    trgt_len=geometry->meanEdgeLength();
    V_bar=geometry->totalVolume();
    
  
    // arcsim::Mesh remesher_mesh;
    // remesher_mesh = translate_to_arcsim(mesh,geometry);


    // std::tie(mesh_uptr, geometry_uptr) = translate_to_geometry(remesher_mesh);
    // mesh = mesh_uptr.release();
    // geometry = geometry_uptr.release();

    

    // arcsim::Simulation sim;
    // sim.cloths.resize(1);
    // arcsim::Cloth Cloth_1;
    // Cloth_1.mesh=remesher_mesh;
    // sim.cloths[0]=Cloth_1;
    // vector<arcsim::AccelStruct *> obs_accs;
    // I need to add more data to the Cloth for it to work
    // arcsim::dynamic_remesh(Cloth_1,false,obs_accs,true);

    

    
    



    ORIG_VPOS = geometry->inputVertexPositions;
    CoM = geometry->centerOfMass();
    double x_furthest=0.0;
    for(Vertex v : mesh->vertices()){
        if(geometry->inputVertexPositions[v].x>x_furthest){
            x_furthest=geometry->inputVertexPositions[v].x;
        }

    }
    x_furthest+=1.4;
    double radius=1.0;
    
    Bead_1 = Bead(mesh,geometry,Vector3({x_furthest,0.0,0.0}),radius,Interaction_str,rc);
    
    // Bead_1.pulling_speed=pullin_force;
    Bead_1.interaction="Shifted_LJ_Normal_nopush";
    M3DG = Mem3DG(mesh,geometry,Bead_1);
    
    // Add visualization options.
    
    std::stringstream nustream;
    std::stringstream rcstream;
    std::stringstream KAstream;
    std::stringstream KBstream;
    std::stringstream Interactionstrstream;
    
    // std::stringstream H0stream;
    // std::stringstream kappastream;
    // std::stringstream sigmastream;
    std::stringstream Curv_adapstream;
    std::stringstream Min_rel_lengthstream;
    


    nustream << std::fixed << std::setprecision(3) << nu;
    rcstream << std::fixed << std::setprecision(3) << rc;
    KAstream << std::fixed << std::setprecision(3) << KA;
    KBstream << std::fixed << std::setprecision(6) << KB;
    Interactionstrstream << std::fixed << std::setprecision(6) << Interaction_str;



    Curv_adapstream << std::fixed << std::setprecision(2) << Curv_adap;
    Min_rel_lengthstream << std::fixed << std::setprecision(2) <<Min_rel_length;
    
    

    std::string first_dir="../Results/Mem3DG_Bead_Pulling_rc_august_arcsim/";
    int status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"If this name is 0 the directory was created succesfully "<< status ;

    std::string basic_name=first_dir+"nu_"+nustream.str()+"_rc_"+rcstream.str()+"_KA_"+KAstream.str()+"_KB_"+KBstream.str()+"_strength_"+Interactionstrstream.str()+"_Init_cond_"+std::to_string(Init_cond)+"_Nsim_"+std::to_string(Nsim)+"/";
    status = mkdir(basic_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    std::cout<<"\nIf this number is 0 the directory was created succesfully "<< status<<"\n" ;

    std::string filename = basic_name+"Output_data.txt";

    std::string filename3 = basic_name+"Output_data.txt";


    std::ofstream Sim_data(filename3);
    Sim_data<<"T_Volume T_Area time Volume Area E_vol E_sur E_bend grad_norm backtrackstep\n";
    Sim_data.close();

    std::string filename2 = basic_name + "Bead_data.txt";
    
    std::ofstream Bead_data(filename2);

    bool Save_bead_data=false;
    bool crash = false;
    bool Save_output_data=false;
    Bead_data<<"####### This data is taken every 250 steps just like the mesh radius is " << radius<<" \n";
    
    std::string filename4 = basic_name + "Bead_data_SS.txt";
    std::ofstream Bead_data_SS(filename4);
    Bead_data_SS<<"#### This data is taken everytime i do an increase in the interacion strength\n";
    Bead_data_SS.close();


    bool flag_SS=true;
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
    double time=0.0;
    double dt_sim=0.0;
    bool arcsim = true;
    arcsim::Cloth Cloth_1;
    arcsim::Cloth::Remeshing remeshing_params;
    double remeshing_elapsed_time=0;
    double integrate_elapsed_time=0;
    double saving_mesh_time=0;

    trgt_len=geometry->meanEdgeLength();
    if(arcsim){
        std::cout<<"Setting remesher params\n";
        remeshing_params.aspect_min=0.2;
        remeshing_params.refine_angle=0.75;
        remeshing_params.refine_compression=0.01;
        remeshing_params.refine_velocity=1.0;
        remeshing_params.size_max=trgt_len*3.0;
        remeshing_params.size_min=trgt_len*0.1;

    
    }

    // start = chrono::steady_clock::now();
    auto start_time_control= chrono::steady_clock::now();
    auto end_time_control= chrono::steady_clock::now();
    start = chrono::steady_clock::now();
    // Save_mesh(basic_name,12345);
    for(size_t current_t=0;current_t<=400000;current_t++ ){
        // for(size_t non_used_var=0;non_used_var<100;)
        // MemF.integrate(TS,sigma,kappa,H0,P,V0);
        if(arcsim){
            // std::cout<<"entering remesher\n";
            start_time_control=chrono::steady_clock::now();
            // std::cout<<"Translating mesh\n";
            

            arcsim::Mesh remesher_mesh2 = translate_to_arcsim(mesh,geometry);
            Cloth_1.mesh=remesher_mesh2;
            // if( current_t>48000 && current_t<49000 ){
            //     arcsim::save_obj(Cloth_1.mesh, basic_name +"Debugging_before_slot.obj");
            // }
            Cloth_1.remeshing=remeshing_params;
            // std::cout<<"remeshing\n";
            arcsim::dynamic_remesh(Cloth_1);
            
            Bead_1 = M3DG.Bead_1;
            prev_force = M3DG.Bead_1.Total_force.norm2();
            prev_strength = M3DG.Bead_1.strength;
            prev_E = M3DG.Bead_1.prev_E_stationary;
            Bead_force = M3DG.Bead_1.Total_force;
            // std::cout<<"accesing m3dg data\n";
            Save_ss = M3DG.Save_SS;
            Stop_increasing = M3DG.stop_increasing;
            Mem_Force = M3DG.Total_force;
            // std::cout<<"we accessed\n";
            delete mesh;
            delete geometry;
            // std::cout<<"translating back?\n";
            // if( current_t>48000 && current_t<49000 ){
            //     arcsim::save_obj(Cloth_1.mesh, basic_name + "Debugging_after.obj" );
            // }
            std::tie(mesh_uptr, geometry_uptr) = translate_to_geometry(Cloth_1.mesh);
            arcsim::delete_mesh(Cloth_1.mesh);

            // std::cout<<"defining pointers\n";
            mesh = mesh_uptr.release();
            geometry = geometry_uptr.release();
            Bead_1 = Bead(mesh,geometry,Bead_1.Pos,radius,Interaction_str,rc);
            Bead_1.prev_force = prev_force;
            Bead_1.interaction = "Shifted_LJ_Normal_nopush";
            Bead_1.strength = prev_strength;
            Bead_1.prev_E_stationary = prev_E;
            
            M3DG = Mem3DG(mesh,geometry,Bead_1);
            M3DG.Save_SS = Save_ss;
            M3DG.stop_increasing = Stop_increasing;
            end_time_control=chrono::steady_clock::now();
            remeshing_elapsed_time+=chrono::duration_cast<chrono::milliseconds>(end_time_control-start_time_control).count();
            // std::cout<<"leaving remesher\n";

        }
        else{
        prev_force = M3DG.Bead_1.Total_force.norm2();
        // Bead_1.prev_force=prev_force;
        M3DG.Bead_1.prev_force = prev_force;
        // if(current_t%10==0 ){
        start_time_control=chrono::steady_clock::now();
        n_vert_old=mesh->nVertices();
        n_vert_new=1;

        counter=0;
        // std::cout<<"Remeshing\n";
        while(n_vert_new!=n_vert_old && counter<10){

        
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
        if(n_vert_new>=15000){
            std::cout<<"Too many vertices, we are gonna crash \n";
            crash=true;
        }
        end_time_control=chrono::steady_clock::now();
        remeshing_elapsed_time+=chrono::duration_cast<chrono::milliseconds>(end_time_control-start_time_control).count();
        // std::cout<<"The number of vertices is "<<n_vert<< " and the timestep is "<< current_t << " (done remeshing)\n";
        }
        

        if( (M3DG.Save_SS && !Stop_increasing) || (Stop_increasing && current_t%500==0)){
            
            SS_index=SS_index+1;
            // Ok lets try to do this
            if(SS_index%2==1){
            Bead_data_SS = std::ofstream(filename4,std::ios_base::app);
            if(Stop_increasing && current_t&500==0 && flag_SS  ){
                flag_SS = false;
                Bead_data_SS << "# STARTING THE CONSTANT TUBE GROWTH THEREFORE THIS COUNTS AS A STATIONARY STATE # # # # # # # # # # # # # # # # # # # # #\n";
            }
            // Here we need to sabe the data
            Bead_data_SS  << Bead_1.Pos.x << " "<< Bead_1.Pos.y << " " << Bead_1.Pos.z << " " << Bead_force.x << " " << Bead_force.y << " " << Bead_force.z << " "<< Mem_Force.x << " "<< Mem_Force.y << " " << Mem_Force.z << " " << prev_strength << "\n";
            Bead_data_SS.close();

            Save_mesh(basic_name,"Steady_state",SS_index);
            }


            M3DG.Save_SS = false;
        }


        // psMesh->remove();
        
        // psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath("dodecahedra"), geometry->inputVertexPositions,
        //                                     mesh->getFaceVertexList(), polyscopePermutations(*mesh));

        // psMesh->setSurfaceColor({0.9607, 0.6627, 0.7215});
        // psMesh->setEdgeWidth(1.0);

        
        if(current_t%1000==0 || crash){
            // std::cout<<"THe bead position is "<< M3DG.Bead_1.Pos<<" \n";
            start_saving = chrono::steady_clock::now();
            std::cout<<"Saving mesh \n";

            Save_mesh(basic_name,current_t);
            Save_bead_data=true;
            Bead_data = std::ofstream(filename2,std::ios_base::app);
            Save_output_data=true;
            end_saving = chrono::steady_clock::now();
            saving_mesh_time+= chrono::duration_cast<chrono::milliseconds>(end_saving-start_saving).count();
            Sim_data = std::ofstream(filename3,std::ios_base::app);
        }
        if(current_t%100==0){
            Save_output_data=true;
        }
        if(current_t%1000==0) {

            end=chrono::steady_clock::now();
            start_report = chrono::steady_clock::now();
            std::cout<<"Started reporting process \n";
            n_vert=mesh->nVertices();
            std::cout<< "THe number of vertices is "<< n_vert <<"\n";    
            std::cout<< "the avg edge lenghth is " <<geometry->meanEdgeLength()<<"\n";
            std::cout << "The avg edge length is = " << std::fixed << std::setprecision(10) << geometry->meanEdgeLength() << std::endl;
            std::cout << "The interaction strength is "<< M3DG.Bead_1.strength <<" \n";
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
            end_report = chrono::steady_clock::now();
            std::cout<<"Reporting Process took "<<chrono::duration_cast<chrono::milliseconds>(end_report-start_report).count()   << " miliseconds\n";
            start = chrono::steady_clock::now();

        }
        
        nu_evol= time<50 ? nu_0 + (nu-nu_0)*time/50 : nu; 
        if(crash){
            break;
        }
        start_time_control = chrono::steady_clock::now();
        dt_sim=M3DG.integrate(TS,V_bar,nu_evol,c0,P0,KA,KB,Kd,Sim_data,time,Save_bead_data,Bead_data,Save_output_data,pulling);
        end_time_control = chrono::steady_clock::now();
        integrate_elapsed_time += chrono::duration_cast<chrono::milliseconds>(end_time_control-start_time_control).count(); 
        Save_output_data=false;
        Save_bead_data=false;
        Bead_data.close();
        Sim_data.close();
        if(current_t%1000==0){
            std::cout<<"THe position of the bead is "<< M3DG.Bead_1.Pos<<"\n";
            std::cout<< "Remeshing has taken a total of "<< remeshing_elapsed_time <<" milliseconds\n" << "Saving the mesh has taken a total of "<< saving_mesh_time<< "milliseconds \n Integrating the forces has taken a total of "<< integrate_elapsed_time <<" milliseconds \n\n"; 
        }
        if(dt_sim==-1){
            std::cout<<"Sim broke or timestep very small\n";
            break;
        }
        else{
            time+=dt_sim;
            
        }





    }
    // Sim_data.close();
    // Bead_data.close();

    Vector3 Pos;
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



