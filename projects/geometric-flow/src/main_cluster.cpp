// GEOMETRIC FLOW

// #include <stdlib.h>
#include <unistd.h>


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


#include <nlohmann/json.hpp>
using json = nlohmann::json;




using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace std;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;
// so we can more easily pass these to different classes
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;



std::vector<arcsim::Mesh> Saved_meshes(6);
std::vector<arcsim::Mesh> Saved_after_remesh(6);


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
E_Handler Sim_handler;
Bead Bead_1;
Bead Bead_2;

Eigen::SparseMatrix<double> Hessian;
Eigen::MatrixXd Hessian_matrix;
std::array<double, 3> BLUE = {0.11, 0.388, 0.89};
// glm::vec<3, float> ORANGE_VEC = {1, 0.65, 0};
std::array<double, 3> ORANGE = {1, 0.65, 0};


void showSelected() {
    // pass
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

void Save_dihedrals(std::string basic_name) {
    // Build member variables: mesh, geometry
    Vector3 Pos;
    std::ofstream o(basic_name+"dihedrals_evol.txt", std::ios::app);

    for(Edge e : mesh->edges()) {
        o << geometry->dihedralAngle(e.halfedge()) << " ";
        
    }
    o<<"\n";
    o.close();


    return ;
}

void Save_edgelengths(std::string basic_name){
    std::ofstream o(basic_name+"edgelengths.txt", std::ios::app);
    
    for(Edge e : mesh->edges()) {
        o << geometry->edgeLength(e) << " ";
        
    }
    o <<"\n";
    o.close();
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
    
    vector<int> flags(0);

    for(size_t v = 0 ; v < mesh.verts.size(); v++){
        arcsim::Vec3 pos_old = mesh.nodes[v]->x;
        
        v_pos.x=pos_old[0];
        v_pos.y=pos_old[1];
        v_pos.z=pos_old[2];
             
        
        simpleMesh.vertexCoordinates.push_back(v_pos);
    }

    // }
    int id1;
    int id2;
    int id3;


    bool non_manifold = false;
    for(size_t f = 0 ; f < mesh.faces.size();f++){
        
        std::vector<size_t> polygon(3);
        
        id1 = mesh.faces[f]->v[0]->index;
        id2 = mesh.faces[f]->v[1]->index;
        id3 = mesh.faces[f]->v[2]->index;

       
        polygon[0] = id1;
        polygon[1] = id2;
        polygon[2] = id3;
        
        simpleMesh.polygons.push_back(polygon);
        
    }
    
    // std::cout<<"Does this happen after loading the data to create the mesh?\n";
    // std::cout<<" THe information in the mesh is, "<< simpleMesh.vertexCoordinates.size()<<"number of vertices\n";
  auto lvals = makeManifoldSurfaceMeshAndGeometry(simpleMesh.polygons, simpleMesh.vertexCoordinates);
 
 return std::tuple<std::unique_ptr<ManifoldSurfaceMesh>,
                    std::unique_ptr<VertexPositionGeometry>>(std::move(std::get<0>(lvals)),  // mesh
                                                             std::move(std::get<1>(lvals))); // geometry
}
std::vector<std::string> split(std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

Vector3 Get_bead_pos(std::string filename, int step){
    std::cout<<"THe filename is" << filename <<" \n";
    // Ok i need to get the vector for the position;
    std::ifstream Bead_data(filename);

    if(!Bead_data.is_open()){
        cerr << "Error opening the file! "<< endl;
        return Vector3({0.0,0.0,0.0}); 
    }
    string line ;
    int counter = 0;
    Vector3 Bead_pos({0.0,0.0,0.0});
    // the frame i want is (500*nsim)/100 = 5*nsim

    while(std::getline(Bead_data,line)){
        std::vector<std::string> splitted = split(line, " ");
        if(line[0]=='#'){
            continue;
        }
        
        if(counter == int(step/50)){
            // We get the things we need
            // Bead_pos.x = std::stod(splitted[0]);
            Bead_pos = Vector3({std::stod(splitted[0]),std::stod(splitted[1]),std::stod(splitted[2])});
            break;
        }  
        counter+=1;            

    }
    if(counter<step/50.0){
        std::cout<<"The output file doesnt get that far\n";
    }

    std::cout<<"The initial bead position is "<< Bead_pos.x << " " << Bead_pos.y <<" " << Bead_pos.z <<" \n";


    return Bead_pos;
}



int main(int argc, char** argv) {

    

    // Here is where we start changing stufff
    // The parsing has to do many stuff before its completely functional
    std::fstream JsonFile;
    
    JsonFile.open(argv[1], std::ios::in);
    int Nsim = std::stoi(argv[2]);

    json Data = json::parse(JsonFile);

    // We loaded the json file 

    // std::cout << Data.dump(1);
    
    std::string Switch = "None";
    size_t Switch_t = 0;

    bool finish_sim = true;

    if(Data.contains("finish_sim")){
        finish_sim = Data["finish_sim"];
    }
    else{
        std::cout<<"The simulation will finish when the timestep decreases\n";
    }

    std::string filepath = Data["init_file"];
    int save_interval = Data["save_interval"];
    bool resize_vol = Data["resize_vol"];
    bool arcsim = Data["arcsim"];
    
    // Ok here 
    

    std::string Integration = "Gradient_descent";
    
    if(Data.contains("Integration")){
        Integration = Data["Integration"];
    }
    else{
        std::cout<<"The integration method is not defined, using Gradient descent\n";
    }

    if(Data.contains("Switch")){
        Switch = Data["Switch"];
        Switch_t = Data["Switch_t"];
    }
    else{
        std::cout<<"No switch in this run";
    }

    int remesh_every = 1;
    if( Data.contains("remesh_every")) remesh_every = Data["remesh_every"];
    std::cout<<"Remesh every is " << remesh_every << std::endl;

    Vector3 Recenter{0.0,0.0,0.0};
    if(Data.contains("Displacement")){
        Recenter.x = Data["Displacement"][0];
        Recenter.y = Data["Displacement"][1];
        Recenter.z = Data["Displacement"][2];
        std::cout<<"Displacing the membrane by " << Recenter <<" \n";
    }

    double scale_factor = 1.0;
    if(Data.contains("rescale")){
        scale_factor = Data["rescale"];
    }

    
    bool Saving_last_states = Data["saving_states"];
    size_t Final_t = Data["timesteps"];

    // Here i will load the geometry 
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    
    // Here i want to recenter and rescale.

    geometry->normalize(Recenter);
    geometry->rescale(scale_factor);
    geometry->refreshQuantities();

    V_bar = geometry->totalVolume();

    // We will deal with the energies now
    std::vector<std::string> Energies(0);
    std::vector<std::vector<double>> Energy_constants(0);
    std::vector<double> Constants(0);



    for( auto Energy : Data["Energies"]){
        Energies.push_back(Energy["Name"]);
        Constants = Energy["constants"].get<std::vector<double>>();
        std::cout<<"The constants for " << Energy["Name"]<< " are ";
        for(size_t z = 0; z < Constants.size(); z++) std::cout<<Constants[z] << " ";
        std::cout<<" \n";

        Energy_constants.push_back(Constants);
        Constants.resize(0);
    }
    
    std::vector<Bead> Beads;
    std::vector<std::string> bonds;
    std::vector<std::vector<double>> constants;

    Vector3 BPos;
    double radius;
    double interaction_str;
    std::string state;
    std::string interaction_mem;
    std::string Constraint;
    std::vector<double> Constraint_constants;
    Bead PBead;
    for( auto Bead_data : Data["Beads"]){
        std::cout<<"Adding a bead\n";
        if(Bead_data.contains("gradient_order")){
            Energies.push_back(Bead_data["gradient_order"]); 
        }
        else{
        Energies.push_back("Bead");
        }
        if(Bead_data.contains("Constraint")){
            Constraint = Bead_data["Constraint"];
            Constraint_constants = Bead_data["Constraint_constants"].get<std::vector<double>>();
        }
        else{
            Constraint = "None";
            Constraint_constants = {};
        }
        
        Energy_constants.push_back(Constants);

        BPos = Vector3({Bead_data["Pos"][0],Bead_data["Pos"][1],Bead_data["Pos"][2] });
        radius = Bead_data["radius"];
        state = Bead_data["state"];
        interaction_str = Bead_data["inter_str"];
        interaction_mem = Bead_data["mem_inter"];
        PBead = Bead(mesh, geometry, BPos, radius, interaction_str);
        PBead.interaction = interaction_mem;
        PBead.Bond_type = Bead_data["bonds"].get<vector<std::string>>();
        PBead.Interaction_constants_vector = Bead_data["bonds_constants"].get<std::vector<std::vector<double>>>();
        PBead.state = state;
        if(Bead_data["mem_inter"]=="Shifted_LJ"){
            PBead.rc = radius*pow(2,1.0/6.0);
        }
        else{
            PBead.rc = 2.0*radius;
        }
        PBead.Constraint = Constraint;
        PBead.Constraint_constants = Constraint_constants;
        std::cout<<"The bead has radius" << PBead.sigma <<" cutoff of " << PBead.rc <<" Interaction of " << PBead.interaction << " \n";
        Beads.push_back(PBead);

    }
    // Once all the beads have been added I need to make them point to each other
    std::vector<int> BeadBonds(0);
    size_t counter = 0;
    for(auto Bead_data: Data["Beads"]){
        // I am iterating again but i only care about the bonds
        BeadBonds = Bead_data["Beads"].get<std::vector<int>>();
        for(size_t i = 0; i < BeadBonds.size(); i++) 
        {Beads[counter].Beads.push_back(&Beads[BeadBonds[i]]);
        counter+=1;
        }
    }
    std::cout<<"There are " << counter/2 <<" bonds \n";

    if(counter>0) std::cout<< "The bead with radius " << Beads[0].sigma <<" is connected to the bead with radius " << Beads[0].Beads[0]->sigma << " \n";

    for(size_t i = 0; i < Beads.size(); i++) std::cout<<"The bead has " << Beads[i].Bond_type.size() << " bonds and interaction "<< Beads[i].interaction << " \n"; 
    // Now the beads point to each other ()

    // Lets define our integrator and all its values

    M3DG = Mem3DG(mesh,geometry);
    Sim_handler = E_Handler(mesh,geometry,Energies, Energy_constants);
    
    M3DG.recentering = Data["recentering"];
    M3DG.boundary = Data["boundary"];
    Sim_handler.boundary = Data["boundary"];

    for( size_t i = 0 ; i< Beads.size() ; i++) {
        M3DG.Add_bead(&Beads[i]);
        Sim_handler.Add_Bead(&Beads[i]);
    }
    
    M3DG.Sim_handler = &Sim_handler;

    if(Data.contains("Field")){
        M3DG.Field = Data["Field"];
        M3DG.Field_vals = Data["Field_vals"].get<std::vector<double>>();
        std::cout<<"The field is " << M3DG.Field << " and the values are ";
        for(size_t i = 0; i < M3DG.Field_vals.size(); i++) std::cout<< M3DG.Field_vals[i] << " ";
        std::cout<<"\n";

    }
    else{
        M3DG.Field = "None";
    }
    
    // Here i will do my alling
    // Face f = mesh->face(1);
    // M3DG.Face_sizing(f);
    // std::cout<<"Thats it\n";
    // return 1;



    arcsim::Cloth Cloth_1;
    arcsim::Cloth::Remeshing remeshing_params;

    if(arcsim){
        // You can make this only one function but you need to write the code for json helpers
        remeshing_params.aspect_min = Data["remesher"]["aspect_min"];
        remeshing_params.refine_angle = Data["remesher"]["refine_angle"];
        remeshing_params.refine_compression = Data["remesher"]["refine_compression"];
        remeshing_params.refine_velocity = Data["remesher"]["refine_velocity"];
        remeshing_params.size_max = Data["remesher"]["size_max"];
        remeshing_params.size_min = Data["remesher"]["size_min"];    
        remeshing_params.total_op = -1;
    }


    

    
    // return 1;

    
    int saved_mesh_idx = 0;
    std::vector<Vector3> Bead_pos_saved(6);

    

    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();

    auto start_full = chrono::steady_clock::now();
    auto end_full = chrono::steady_clock::now();
    
    auto start_time_control= chrono::steady_clock::now();
    auto end_time_control= chrono::steady_clock::now();
    
    double remeshing_elapsed_time=0;
    double integrate_elapsed_time=0;
    double saving_mesh_time=0;


    std::cout<< "Current path is " << argv[0];

    std::cout<<"The energy elements are \n";
    
    for(size_t z = 0 ; z < Energies.size(); z++){
        std::cout<<Energies[z]<<" ";
    }
    std::cout<<"\n";

    std::cout<<"Minimum edge length allowed is "<< remeshing_params.size_min<<" muak\n";


    double avg_dih = 0;
    double max_dih = 0;
    double min_dih = 0.1;
    double dih;
        

    for( Edge e : mesh->edges()){ 
    dih = fabs(geometry->dihedralAngle(e.halfedge()));
    avg_dih+=dih;
    if(dih > max_dih) max_dih = dih;
    if(dih < min_dih) min_dih = dih;
    }
    std::cout<<"First checkThe average dihedral is"<< avg_dih/mesh->nEdges()<<" \n";
    std::cout<<"The min dih is"<< min_dih << " and the max dih is " << max_dih <<" \n";
    avg_dih =0.0;

    
    FaceData<double> F_sizings = M3DG.Face_sizings();
    VertexData<double> Sizings = M3DG.Vert_sizing(F_sizings);
    double max_sizing = 0.0;
    double min_sizing = 1e4;
    double sizing;
    
    for(Vertex v : mesh->vertices()){
        // 
        sizing = Sizings[v];
        if(sizing > max_sizing ) max_sizing = sizing;
        if(sizing < min_sizing ) min_sizing = sizing;

    }
    std::cout<<"My algorithm says\n";
    std::cout<<"The max sizing is " << max_sizing << " and the min sizing is " << min_sizing <<" \n";
    arcsim::Mesh remesher_mesh = translate_to_arcsim(mesh,geometry);
    Cloth_1.mesh = remesher_mesh;
    Cloth_1.remeshing = remeshing_params; 
    arcsim::compute_masses(Cloth_1);
    arcsim::compute_ws_data(Cloth_1.mesh);

    std::cout<<"ARCSIM says \n";
    arcsim::dynamic_remesh(Cloth_1);
    // arcsim::static_remesh2(Cloth_1,true);

    std::cout<<"Done first remeshing\n";

    std::cout<<" \n\n";


    // return 1;

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
  

    ORIG_VPOS = geometry->inputVertexPositions;
    CoM = geometry->centerOfMass();
            





    // How do we create a name that makes sense 
    // I can do something like E _ param_ param_ E_ param_param _Nsim_
    std::string Directory = "";
    std::stringstream stream;
    std::string s;
    int bead_counter = 0;
    for(size_t z = 0; z < Energies.size(); z++){
        Directory = Directory + Energies[z]+"_";
        
        if(Energies[z] == "Bead" || Energies[z] == "H1_Bead" || Energies[z] == "H2_Bead" ){
            stream.str(std::string());
            stream << std::fixed << std::setprecision(4) << Beads[bead_counter].sigma;
            Directory = Directory +"radius_" +stream.str() +"_";

            stream.str(std::string());
            stream << std::fixed << std::setprecision(4) << Beads[bead_counter].interaction;
            Directory = Directory + stream.str()+"_";

            stream.str(std::string());
            stream << std::fixed << std::setprecision(4) << Beads[bead_counter].strength;
            Directory = Directory + "str_" +stream.str() +"_";

            if(Beads[bead_counter].Constraint == "Radial"){
                stream.str(std::string());
                stream << std::fixed << std::setprecision(4) << Beads[bead_counter].Constraint_constants[0];
                Directory = Directory + "theta_const_" +stream.str() +"_";
            }
            // I want to add the theta because if not there is no way of knowing the difference.
            // if(Beads
            

            bead_counter +=1 ;
        }
        for(size_t j = 0; j < Energy_constants[z].size(); j++) 
        {
            stream.str(std::string());
            stream << std::fixed << std::setprecision(4) << Energy_constants[z][j];
        
            Directory = Directory + stream.str() +"_";
            // stream << std::fixed << std::setprecision(2) << pi;
        }
        
    }
    
    // Lets add the fields

    if(Data.contains("Field")){
        
        Directory = Directory + M3DG.Field + "_";
        for(size_t i = 0; i < M3DG.Field_vals.size(); i++){
            stream.str(std::string());
            stream << std::fixed << std::setprecision(4) << M3DG.Field_vals[i];
            Directory = Directory + stream.str() + "_";
        }
    }

    // Lets add the bonds the types and the interaction strength

    bool bonds_exist = false;
    for(size_t z = 0; z< Beads.size(); z++){
        if(Beads[z].Bond_type.size()>0 && bonds_exist == false) {
            Directory = Directory + "Bonds_"; 
            bonds_exist = true;
        }
        for(size_t j = 0; j < Beads[z].Bond_type.size(); j ++){
            stream.str(std::string());
            stream << std::fixed << std::setprecision(4) << Beads[z].Interaction_constants_vector[j][0];
            Directory = Directory + Beads[z].Bond_type[j] + "_"+ stream.str() + "_";
        }
    }

    if(Switch !="None"){
        Directory = Directory + "Switch_" + Switch + "_";
        stream.str(std::string());
        stream << std::fixed << std::setprecision(1) << Switch_t;
        Directory = Directory + "Switch_t_" + stream.str() + "_";
    }


    
    Directory = Directory+ "Nsim_" + std::to_string(Nsim)+"/";
    std::cout<<"Directory is " << Directory << " \n";
        
    std::string first_dir = Data["first_dir"];

    int status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    

    std::string basic_name=first_dir+Directory;
    
    
    status = mkdir(basic_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    std::cout<<"\nIf this number is 0 the directory was created succesfully "<< status<<"\n" ;

    std::string filename = basic_name+"Output_data.txt";

    std::ofstream Sim_data(filename);
    Sim_data<<"T_Volume T_Area time Volume Area E_vol E_sur E_bend grad_norm backtrackstep\n";
    Sim_data.close();


    std::vector<std::string> Bead_filenames;
    std::ofstream Bead_datas;

    

    for(size_t i = 0 ; i < Beads.size(); i++){
        Bead_filenames.push_back(basic_name+ "Bead_"+std::to_string(i)+"_data.txt");
        Bead_datas = std::ofstream(Bead_filenames[i]);
        
        Bead_datas<<"####### This data is taken ever y" << save_interval <<" steps just like the mesh radius is " << radius<<" \n";
        Bead_datas.close();
    }

    // Here

    Bead_filenames.push_back(basic_name+ "Simulation_timings.txt");
    Bead_datas = std::ofstream(Bead_filenames[Beads.size()]);
    Bead_datas <<"Remeshing_time Gradients Backtracking  Construction compute Solve \n";
    Bead_datas.close();


    std::string filename2 = basic_name + "Bead_data.txt";
    
    std::ofstream Bead_data(filename2);

    bool Save_bead_data=false;
    bool Save_output_data=false;
    bool small_Ts;
    Bead_data<<"####### This data is taken every"<< save_interval <<" steps just like the mesh dump, radius is " << radius<<" \n";
    Bead_data.close();
    // Here i want to run my video
    size_t n_vert;
    size_t n_vert_old;
    size_t n_vert_new;
    double Volume;
    double Area;
    double nu_obs;
    // double nu_evol;
    double nu_0;
    double c_null;
    counter=0;
    double time=0.0;
    double dt_sim=0.0;
    int sys_time = 0;





    bool seam = false;
    Cloth_1.dump_info = false;
    start = chrono::steady_clock::now();

    // Save_mesh(basic_name,-1);



    // Save_mesh(basic_name,111);
    // M3DG.Smooth_vertices();
    // geometry->refreshQuantities()
    // Save_mesh(basic_name,112);

    
    // std::ofstream Dihedrals(basic_name+"dihedrals_evol.txt");
    // Dihedrals << "## THE EVOLUTION OF THE DIHEDRAL DISTRIBUTION\n";
    // Dihedrals.close();
    // std::ofstream EdgeLengths(basic_name+"edgelengths.txt");
    // EdgeLengths << "## THE EVOLUTION OF THE EDGE LENGTHS\n";
    // EdgeLengths.close();



    start_full = chrono::steady_clock::now();
    for(size_t current_t=0; current_t <= Final_t ;current_t++ ){

        // std::cout<<"Curren t t is " << current_t <<" \n";
        // if(current_t>400){
        //     save_interval = 1;
        // }
        if(Switch =="Newton" && current_t >= Switch_t){
            Integration = "Newton";
            remesh_every = -1;
            save_interval = 10;
        }
        
        if(Switch=="Free_beads" && current_t>= Switch_t){

            std::cout<<"Switching the beahaviour of the beads\n";
            // We activate the switch
            for(size_t i = 0; i < Beads.size(); i++){
                // We change their state to default
                Beads[i].state = "default";

            }

            std::vector<double> Interaction_const(0);
            Interaction_const.push_back(1.0);
            Beads[0].Add_bead( &Beads[1],"Shifted-LJ",Interaction_const);
            Beads[1].Add_bead( &Beads[0],"Shifted-LJ", Interaction_const);


            std::cout <<"The bonds of the first bead are "<<  Beads[0].Bond_type[0] << " with size" << Beads[0].Bond_type.size() <<" \n";
            
            // So now they are free to move (this still does not add the interaction with each other)

            Switch = false;
        }
        if(Switch=="No_remesh" && current_t>= Switch_t){

            std::cout<<"Switching off the remesher \n";
            

            arcsim = false;
            // So now they are free to move (this still does not add the interaction with each other)

            Switch = false;
        }

        start_time_control=chrono::steady_clock::now();

        if( arcsim  && (current_t%remesh_every==0 || dt_sim == 0.0)){
            // if(dt_sim==0.0){ std::cout<<"wE ARE REMESHING CAUSE THINGS DONT MAKE SENSE\n";}
            // std::cout<<"\t\t Remeshing\n";
            // Here we want to explore whats going on 
            


            int n_vert_old=0;
            int n_vert_new=0;

            n_vert_old  =  mesh->nVertices();
            

            double avg_dih1;
            double max_dih1 = 0.0;
            double min_dih1 = 1e2;
            double avg_dih2;
            double max_dih2 = 0.0;
            double min_dih2 = 1e2;
            // M3DG.Smooth_vertices();
            // avg_dih2 = 0.0;
            // avg_dih1 = 0.0;

            // for( Edge e : mesh->edges()){ 
            //     dih = fabs(geometry->dihedralAngle(e.halfedge()));
            //     avg_dih1+=dih;
            //     if(dih > max_dih1) max_dih1 = dih;
            //     if(dih < min_dih1) min_dih1 = dih;
            //     }
            // // std::cout<<"The average dihedral is"<< avg_dih/mesh->nEdges()<<" \n";
            // // std::cout<<"The min dih is"<< min_dih << " and the max dih is " << max_dih <<" \n";
            // // avg_dih =0.0;
            // avg_dih1 = avg_dih1/mesh->nEdges();


            arcsim::Mesh remesher_mesh2 = translate_to_arcsim(mesh,geometry);
            Cloth_1.mesh=remesher_mesh2;

            if(Saving_last_states){
                saved_mesh_idx = (saved_mesh_idx+1)%6;
                arcsim::delete_mesh(Saved_meshes[saved_mesh_idx]);
                if(Beads.size()>=01){
                Bead_pos_saved[saved_mesh_idx] = Beads[0].Pos;
                }
                Saved_meshes[saved_mesh_idx] = arcsim::deep_copy(remesher_mesh2);

            }

            
            Cloth_1.remeshing=remeshing_params;
            arcsim::compute_masses(Cloth_1);
            arcsim::compute_ws_data(Cloth_1.mesh);
            arcsim::dynamic_remesh(Cloth_1);

            // std::cout<<"The operation counter is "<< Cloth_1.remeshing.op_counter <<"\n";
            if(Cloth_1.remeshing.op_counter != 0){
                arcsim::compute_masses(Cloth_1);
                arcsim::compute_ws_data(Cloth_1.mesh);
                arcsim::dynamic_remesh(Cloth_1);
            }
            
            

            if(Saving_last_states){
                arcsim::delete_mesh(Saved_after_remesh[saved_mesh_idx]);
                Saved_after_remesh[saved_mesh_idx] = arcsim::deep_copy(Cloth_1.mesh);
   
            }
   
            small_Ts = M3DG.small_TS;
            sys_time = M3DG.system_time;

            delete mesh;
            delete geometry;
            // std::cout<<"translating back?\n";
            std::tie(mesh_uptr, geometry_uptr) = translate_to_geometry(Cloth_1.mesh);
            arcsim::delete_mesh(Cloth_1.mesh);

            mesh = mesh_uptr.release();
            geometry = geometry_uptr.release();
            
            
            M3DG.mesh = mesh;
            M3DG.geometry = geometry;
            Sim_handler.mesh = mesh;
            Sim_handler.geometry = geometry;

            // M3DG.Smooth_vertices();
            // geometry->refreshQuantities();


            for( Edge e : mesh->edges()){ 
                dih = fabs(geometry->dihedralAngle(e.halfedge()));
                avg_dih2+=dih;
                if(dih > max_dih2) max_dih2 = dih;
                if(dih < min_dih2) min_dih2 = dih;
                }
            
            avg_dih2 = avg_dih2/mesh->nEdges();
            
            n_vert_new = mesh->nVertices();


            for(size_t i = 0 ; i<Beads.size(); i++){
                Beads[i].Reasign_mesh(mesh,geometry);
            }

            if( abs(n_vert_new-n_vert_old)>= 300){
                std::cout<<"The change in the number of vertices is "<< n_vert_new-n_vert_old <<" \n";
                std::cout<<"Which is clearly too big :P \n";    
                break;
            }

            // Ok so here we have remeshed, and the idea is that everytime u remesh you have to recompute the lagrange multipliers
            if(Integration == "Newton" && current_t> Switch_t){

                std::cout<<"\t \t At this step\n";
            Sim_handler.Calculate_Jacobian();
            std::cout<<"\t \t The Jacobian is calculated\n";
            Sim_handler.Calculate_gradient();
            std::cout<<"\t \t The gradient is calculated\n";

            // std::vector<double> RHS(0);
            Eigen::VectorXd RHS(3* mesh->nVertices());
            int index;
            for(Vertex v: mesh->vertices()){
            
                index = v.getIndex();            
                RHS(3*index) = Sim_handler.Current_grad[v].x;
                RHS(3*index+1) = Sim_handler.Current_grad[v].y;
                RHS(3*index+2) = Sim_handler.Current_grad[v].z;
                
            }
            std::cout<<"Doing the solve now\n";

            std::cout<<"The dimension of the jacobian is " << Sim_handler.Jacobian_constraints.rows() << " x " << Sim_handler.Jacobian_constraints.cols() << "\n";
            std::cout<<"The dimension of the RHS is " << RHS.size() << "\n";
            // std::cout<<"The Jacobian is \n" << Sim_handler.Jacobian_constraints.transpose() << "\n";

            // std::cout<<"THe RHS is \n" << RHS.transpose() << "\n";
            Eigen::MatrixXd J_transpose = Sim_handler.Jacobian_constraints.transpose();
            Eigen::VectorXd  Lagrange_mul = J_transpose.colPivHouseholderQr().solve(RHS);
            // Eigen::VectorXd  Lagrange_mul = (J_transpose.transpose() * J_transpose).ldlt().solve(J_transpose.transpose()*RHS);
            // Sim_handler.Lagrange_multipliers = Lagrange_mul;
            std::cout<<"SOlve done\n";
            std::cout<<"THe previous Lagrange multipliers were " << Sim_handler.Lagrange_mult.transpose() <<" \n";
            Sim_handler.Lagrange_mult = Lagrange_mul;
            std::cout<<"THe new  Lagrange multipliers were " << Sim_handler.Lagrange_mult.transpose() <<" \n";
        }
             
        }
        end_time_control=chrono::steady_clock::now();
        remeshing_elapsed_time+=chrono::duration_cast<chrono::milliseconds>(end_time_control-start_time_control).count();
        Bead_datas = std::ofstream(Bead_filenames[Beads.size()],std::ios_base::app);

        Bead_datas << std::chrono::duration_cast<std::chrono::milliseconds>(end_time_control-start_time_control).count()  <<" ";
        Bead_datas.close();
                
        if(current_t%save_interval == 0 || (Integration == "Newton" && current_t%save_interval == save_interval-1 )){
            // Bead_data.close();
            // Sim_data.close();
            // std::cout<<"Saving\n";

            start_time_control=chrono::steady_clock::now();
            // if(current_t%100==0){
            Save_mesh(basic_name,current_t);
            // }
            end_time_control = chrono::steady_clock::now();
            saving_mesh_time+=chrono::duration_cast<chrono::milliseconds>(end_time_control-start_time_control).count();
            Save_bead_data=true;
            Bead_data = std::ofstream(filename2,std::ios_base::app);
            Save_output_data=true;
            Sim_data = std::ofstream(filename,std::ios_base::app);
            

        }
        
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
            
            for( Edge e : mesh->edges()){ 
                dih = fabs(geometry->dihedralAngle(e.halfedge()));
                avg_dih+=dih;
                if(dih > max_dih) max_dih = dih;
                if(dih < min_dih) min_dih = dih;
                }
            std::cout<<"The average dihedral is"<< avg_dih/mesh->nEdges()<<" \n";
            std::cout<<"The min dih is"<< min_dih << " and the max dih is " << max_dih <<" \n";
            avg_dih =0.0;
            

            std::cout<<"A thousand iterations took "<<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" miliseconds\n\n";

            
            // polyscope::screenshot(basic_name+std::to_string(current_t)+".jpg",true);
            // std::cout<<"SMALL TS?\n";
            start = chrono::steady_clock::now();
            if(M3DG.small_TS){
                    std::cout<<"Finishing\n";
            //     break;
            } 
            // std::cout<<"No\n";
    

        }
        // std::cout<<"Redeclaring M3DG\n";
        // std::cout<<"Bead_1 position changed? "<< Bead_1.Pos << " \n";
        
        
        start_time_control = chrono::steady_clock::now();
        // std::cout<<"3\n";
        // std::cout<<"Integrating\n";
        
        // dt_sim=M3DG.integrate(TS,V_bar,nu_evol,c0,P0,KA,KB,sigma,Sim_data, time,Save_bead_data,Bead_filenames,Save_output_data,pulling);
        geometry->refreshQuantities();
        mesh->compress();
        // std::cout<<"We integrate\n";
        if(Integration == "Gradient_descent"){
            dt_sim = M3DG.integrate(Sim_data, time, Bead_filenames, Save_output_data);
        }
        else if(Integration == "BFGS"){
            // std::cout<<"Integrating BFGS\n";
            if(current_t%remesh_every == 0){
                Hessian_matrix = Eigen::MatrixXd::Identity(mesh->nVertices()*3, mesh->nVertices()*3);
             }
                bool Save_Hessian = false;
                if(Save_Hessian){
                    // std::cout<<"Saving Hessian\n";
                    std::ofstream Hess_data(basic_name+"Hessian_"+std::to_string(current_t)+".txt");

                    // Now i just need to save the whole matrix
                    for(size_t i = 0; i < mesh->nVertices()*3; i++){
                        for(size_t j = 0; j < mesh->nVertices()*3; j++){
                            Hess_data << Hessian_matrix(i,j) << " ";
                        } 
                        Hess_data << "\n";
                    }

                    // for (int k=0; k<Hessian.outerSize(); ++k)
                    // for (SparseMatrix<double>::InnerIterator it(Hessian,k); it; ++it)
                    // {
                    //     Hess_data<< it.row() << " "<< it.col() << " " << it.value() << "\n";
                    // }
                    Hess_data.close();
                    // std::cout<<"Succesfully saved\n";
                }

                
                
            Hessian_matrix = M3DG.integrate_BFGS(Sim_data, time, Bead_filenames, Save_output_data, Hessian_matrix);
            // Hessian = M3DG.integrate_BFGS(Sim_data, time, Bead_filenames, Save_output_data, Hessian);
        
        }
        else if(Integration == "Newton"){
            if(current_t == 0 || current_t == Switch_t){
                std::cout<<"defining lagrange mults\n";
                Eigen::VectorXd Lagrange_mults(7);
                Lagrange_mults(0) = 0.0;
                Lagrange_mults(1) = 0.0;
                Lagrange_mults(2) = 0.0;
                Lagrange_mults(3) = 0.0;
                Lagrange_mults(4) = 0.0;
                Lagrange_mults(5) = 0.0;
                Lagrange_mults(6) = 0.0;
                Sim_handler.Lagrange_mult = Lagrange_mults;
                Sim_handler.Trgt_vol = geometry->totalVolume();
            }

            std::vector<std::string> Constraints(0);
            Constraints = std::vector<std::string>{"Volume","CMx","CMy","CMz","Rx","Ry","Rz"};
            Sim_handler.Constraints = Constraints;
            std::vector<std::string> Data_filenames(0);
            if(Integration == "Newton" && ( current_t%remesh_every == remesh_every-1  || current_t%remesh_every == 0) ){
                // This is the integration that happens just before the remeshing
                Save_output_data = true;
                Data_filenames.push_back(basic_name+"RHS_Norm"+std::to_string(current_t)+".txt"); 
            }
            M3DG.integrate_Newton(Sim_data, time, Energies, Save_output_data, Constraints, Data_filenames);
            // std::cout<<"INtegrated newton (: \n";
        }
        // if(dt_sim==0){
        //     Save_mesh(basic_name,-1);
        //     // std::cout<<"THe simulation went crazy i guess? " << dt_sim <<" \n"; 
        // }
        if (M3DG.small_TS && current_t>Final_t*0.2 && finish_sim ) {
            std::cout<<"The  current t is" << current_t << " and the condition is to be grater than " << Final_t*0.2 <<" \n";
            std::cout << "Ending sim due to small TS \n";
            break;
        }
        if(dt_sim<0){
            std::cout<<"Sim broke or timestep very small\n";
            std::cout<<"At timestep " << current_t << " \n";
            break;
        }
        else{
            // std::cout<<"Adding time\n";
            time+=dt_sim;
            M3DG.system_time+=1;
            // std::cout<<"SUccesfully\n";
        }
        if(time>10 && Beads.size()>1 ){
            Beads[1].state="froze";
        }


        // nanvertex = false;
        // for(Vertex v : mesh->vertices()) if(isnan(geometry->inputVertexPositions[v].x+ geometry->inputVertexPositions[v].y  + geometry->inputVertexPositions[v].z )) nanvertex = true;  

        // if(nanvertex) std::cout<< "After integrating one vertex is nan :( also the value of alpha is"<< dt_sim << " \n";




        Bead_data.close();
        Sim_data.close();
        

        // Then i need to multiply all the vertices by this value
        if(resize_vol){
            // std::cout<<"are we resizing?\n";
        double  k;
        
        double V;
        V = geometry->totalVolume();

        k = pow(V_bar/V,1.0/3.0);
        
        geometry->rescale(k);
        // geometry->inputVertexPositions *=k;
        geometry->refreshQuantities();
        
        
        }
        // std::cout<<"The current volume is " << geometry->totalVolume() << " \n";


        end_time_control = chrono::steady_clock::now();
        // std::cout<<"5\n";
        
        integrate_elapsed_time += chrono::duration_cast<chrono::milliseconds>(end_time_control-start_time_control).count(); 
        Save_output_data=false;
        Save_bead_data=false;


        if(current_t%1000==0){
            std::cout<< "Remeshing has taken a total of "<< remeshing_elapsed_time <<" milliseconds\n" << "Saving the mesh has taken a total of "<< saving_mesh_time<< "milliseconds \n Integrating the forces has taken a total of "<< integrate_elapsed_time <<" milliseconds \n\n"; 

        if(M3DG.small_TS && current_t>Final_t*0.2){ 
            std::cout<<"Ending sim due to small TS and long enough time \n";
            break;
        }
        }

        


    }
    end_full = chrono::steady_clock::now();
    Sim_data.close();
    Bead_data.close();

    Bead_data = std::ofstream(Bead_filenames[Beads.size()],std::ios_base::app);

    Bead_data << (chrono::duration_cast<chrono::milliseconds>(end_full-start_full).count()) << " 0 0 0 0 0 \n";
    Bead_data.close();

    if(Saving_last_states)
    {
        std::ofstream beads_saved(basic_name+"Saved_bead_info.txt");
        beads_saved << std::setprecision(std::numeric_limits<double>::max_digits10);

        std::cout<<"Printing bead info\n";
        for( int k = 0 ; k < 6 ; k++){
            std::cout<<Bead_pos_saved[k].x << " " << Bead_pos_saved[k].y << " " << Bead_pos_saved[k].z <<" \n";
            std::cout<<"Saved mesh idx is " << saved_mesh_idx<<" \n";
            arcsim::save_obj(Saved_meshes[saved_mesh_idx],basic_name+"Saved_final_frame_"+std::to_string(11-2*k)+".obj");
            arcsim::save_obj(Saved_after_remesh[saved_mesh_idx],basic_name+"Saved_final_frame_"+std::to_string(12-2*k)+".obj");

            beads_saved<< Bead_pos_saved[saved_mesh_idx].x<<" "<<Bead_pos_saved[saved_mesh_idx].y<<" "<<Bead_pos_saved[saved_mesh_idx].z<<" \n";

            
            saved_mesh_idx = (saved_mesh_idx -1 +6)%6;

            }   
    }
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



