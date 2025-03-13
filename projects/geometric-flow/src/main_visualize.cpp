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



arcsim::Cloth Cloth_1;
arcsim::Cloth::Remeshing remeshing_params;


SimplePolygonMesh Saved_mesh;

std::vector<arcsim::Mesh> Saved_meshes(6);
int counter_saved = 0;
std::vector<Vector3> Saved_beadpos(6);
VertexData<Vector3> Saved_vertex_positions;

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

std::vector<Bead> Beads;
std::vector<std::string> bonds;
std::vector<std::vector<double>> constants;

Bead Bead_1;
Bead Bead_2;


std::vector<std::string> Energies(0);
std::vector<std::vector<double>> Energy_constants(0);
std::vector<double> Constants(0);
std::ofstream Sim_data;
std::vector<std::string> Bead_filenames;
std::ofstream Bead_datas;



polyscope::SurfaceMesh* psMesh;
double TOTAL_ANGLE_DEFECT;
size_t EULER_CHARACTERISTIC;

std::map<int, std::vector<std::array<double, 3>>> sColors;
std::map<int, std::vector<std::array<double, 3>>> sColors2;

polyscope::SurfaceVertexColorQuantity* vertexColors;
polyscope::SurfaceFaceColorQuantity* faceColors;

enum shadeTpe{ sSHADED,sVertsizing, sFacesizing};
std::vector<std::string> sNames = {"Off", "VertSizing", "FaceSizing"};


VertexData<double> Vert_sizings;
FaceData<double> Face_sizings;
EdgeData<double> Edge_sizings;
std::vector<std::array<double,3>> Vert_sizing_C;
std::vector<std::array<double,3>> F_sizings_C;
std::vector<std::array<double,3>> shaded_C;

double scaling_factor = 1.0;



std::array<double, 3> BLUE = {0.11, 0.388, 0.89};
// glm::vec<3, float> ORANGE_VEC = {1, 0.65, 0};
std::array<double, 3> ORANGE = {1, 0.65, 0};


void showSelected() {
    // pass
}

void redraw() {
    psMesh->updateVertexPositions(geometry->inputVertexPositions);
    polyscope::requestRedraw();
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

SimplePolygonMesh save_geometry( ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry){

SimplePolygonMesh simpleMesh;

    // std::cout<<"This is being called\n";
//   processLoadedMesh(simpleMesh, loadType);
    Vector3 v_pos;
    
    vector<int> flags(0);
    for(Vertex v: mesh->vertices()){
        v_pos = geometry->inputVertexPositions[v];
        simpleMesh.vertexCoordinates.push_back(v_pos);
    }

    // }
    int id1;
    int id2;
    int id3;


    bool non_manifold = false;

    for(Face f: mesh->faces()){
        std::vector<size_t> polygon(3);

        Halfedge he = f.halfedge();
        id1 = he.vertex().getIndex();
        he = he.next();
        id2 = he.vertex().getIndex();
        he = he.next();
        id3 = he.vertex().getIndex();
        polygon[0] = id1;
        polygon[1] = id2;
        polygon[2] = id3;

        simpleMesh.polygons.push_back(polygon);
    }

    return simpleMesh;
    // std::cout<<"Does this happen after loading the data to create the mesh?\n";
    // std::cout<<" THe information in the mesh is, "<< simpleMesh.vertexCoordinates.size()<<"number of vertices\n";
//   auto lvals = makeManifoldSurfaceMeshAndGeometry(simpleMesh.polygons, simpleMesh.vertexCoordinates);
 
//  return std::tuple<std::unique_ptr<ManifoldSurfaceMesh>,
                    // std::unique_ptr<VertexPositionGeometry>>(std::move(std::get<0>(lvals)),  // mesh
                                                            //  std::move(std::get<1>(lvals))); // geometry
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


void Remesh(){
       int n_vert_old=0;
            int n_vert_new=0;

            n_vert_old  =  mesh->nVertices();
            double dih;
            double small_Ts;
            double sys_time;


            double avg_dih1;
            double max_dih1 = 0.0;
            double min_dih1 = 1e2;
            double avg_dih2;
            double max_dih2 = 0.0;
            double min_dih2 = 1e2;
            // M3DG.Smooth_vertices();
            avg_dih2 = 0.0;
            avg_dih1 = 0.0;

            for( Edge e : mesh->edges()){ 
                dih = fabs(geometry->dihedralAngle(e.halfedge()));
                avg_dih1+=dih;
                if(dih > max_dih1) max_dih1 = dih;
                if(dih < min_dih1) min_dih1 = dih;
                }
            // std::cout<<"The average dihedral is"<< avg_dih/mesh->nEdges()<<" \n";
            // std::cout<<"The min dih is"<< min_dih << " and the max dih is " << max_dih <<" \n";
            // avg_dih =0.0;
            avg_dih1 = avg_dih1/mesh->nEdges();


            arcsim::Mesh remesher_mesh2 = translate_to_arcsim(mesh,geometry);
            Cloth_1.mesh=remesher_mesh2;
            
            Cloth_1.remeshing=remeshing_params;
            arcsim::compute_masses(Cloth_1);
            arcsim::compute_ws_data(Cloth_1.mesh);
            arcsim::dynamic_remesh(Cloth_1);
            
               
            small_Ts = M3DG.small_TS;
            sys_time = M3DG.system_time;

            delete mesh;
            delete geometry;
            
            std::tie(mesh_uptr, geometry_uptr) = translate_to_geometry(Cloth_1.mesh);
            arcsim::delete_mesh(Cloth_1.mesh);

            mesh = mesh_uptr.release();
            geometry = geometry_uptr.release();
            
            
            M3DG.mesh = mesh;
            M3DG.geometry = geometry;
            

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
                // break;
            }

            geometry->requireVertexPositions();
            psMesh = polyscope::registerSurfaceMesh("MyMesh", geometry->vertexPositions,mesh->getFaceVertexList());

            return;

        }



void computeColors(){

    Face_sizings = M3DG.Face_sizings();
    Vert_sizings = M3DG.Vert_sizing(Face_sizings);

    double F_sizing_max = 0.0;
    double V_sizing_max = 0.0;

    double F_sizing_min = 1e3;
    double V_sizing_min = 1e3;

    for(Vertex v: mesh->vertices()){
        if(Vert_sizings[v] > V_sizing_max ) V_sizing_max = Vert_sizings[v];
        if(Vert_sizings[v] < V_sizing_min ) V_sizing_min = Vert_sizings[v]; 
    }
    for(Face f: mesh->faces()){
        if(Face_sizings[f] > F_sizing_max ) F_sizing_max = Face_sizings[f];
        if(Face_sizings[f] < F_sizing_min ) F_sizing_min = Face_sizings[f]; 
    }
    
    // We have the min and the max

    for(Vertex v: mesh->vertices()){

        shaded_C.push_back({1.0, 0.45, 0.0});

        // Vert_sizing_C.push_back(mapToColor(Vert_sizings[v],V_sizing_min,V_sizing_max,"viridis"));



    }


}


static const ImVec4 pressColor = ImColor::HSV(1. / 7.0f, 0.6f, 0.6f); // gold
static const ImVec4 releaseColor{0.35f, 0.61f, 0.49f, 0.62f};         // default green

static ImVec4 off_currColor = pressColor;
static ImVec4 V_sizing_currColor = releaseColor;
static ImVec4 F_sizing_currColor = releaseColor;

std::vector<ImVec4*> sState = {&off_currColor, &V_sizing_currColor, &F_sizing_currColor};


void functionCallback() {

    ImGui::Text("Total angle defect: %0.1fpi", TOTAL_ANGLE_DEFECT / M_PI);
    ImGui::Text("Euler characteristic: %zu", EULER_CHARACTERISTIC);
    ImGui::Text("");

    if(ImGui::Button("Remesh")){
        Remesh();
        redraw();
    }

    ImGui::PushStyleColor(ImGuiCol_Button, *sState[0]);
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, *sState[0]);
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, *sState[0]);
    
   if (ImGui::Button("Vertex Sizing")) {
    Face_sizings = M3DG.Face_sizings();
    Vert_sizings = M3DG.Vert_sizing(Face_sizings);
    psMesh->addVertexScalarQuantity("Vert sizing", Vert_sizings);
//         vertexColors = psMesh->addVertexColorQuantity("Plot", Vert_sizings);
//         // for (size_t j = 0; j < sColors.size(); j++) {
//             // *sState[j] = i == j ? pressColor : releaseColor;
        
        }
    if(ImGui::Button("Face Sizing")){
        Face_sizings = M3DG.Face_sizings();
        psMesh->addFaceScalarQuantity("Face sizing", Face_sizings);
    }
    if(ImGui::Button("Edge sizing")){
        Face_sizings = M3DG.Face_sizings();
        Vert_sizings = M3DG.Vert_sizing(Face_sizings);    
        Edge_sizings = M3DG.Edge_sizing(Vert_sizings);
        psMesh->addEdgeScalarQuantity("Edge sizing", Edge_sizings);

    }
    if(ImGui::Button("Save current state")){

        Saved_mesh = save_geometry(mesh,geometry);

        // Here i need the state and the vectors
        std::cout<<"THis button is being pressed\n";
        // arcsim::delete_mesh(Saved_meshes[0]);
        // Saved_vertex_positions = geometry->inputVertexPositions;
        // arcsim::Mesh remesher_mesh2 = translate_to_arcsim(mesh,geometry);
        // Saved_meshes[0] =  arcsim::deep_copy(remesher_mesh2);
        Saved_beadpos[0] = Beads[0].Pos;
        // std::cout<<"Saved mesh in a state\n";

    }
    if(ImGui::Button("Reload saved state")){

        std::cout<<"reloading\n";
        delete mesh;
        delete geometry;
        std::cout<<"deleted\n";
        
        auto lvals = makeManifoldSurfaceMeshAndGeometry(Saved_mesh.polygons, Saved_mesh.vertexCoordinates);

 
        std::tie(mesh_uptr, geometry_uptr) = std::tuple<std::unique_ptr<ManifoldSurfaceMesh>,
                    std::unique_ptr<VertexPositionGeometry>>(std::move(std::get<0>(lvals)),  // mesh
                                                             std::move(std::get<1>(lvals))); 

        std::cout<<"reassigned\n";
        
        mesh = mesh_uptr.release();
        geometry = geometry_uptr.release();
            
            
        M3DG.mesh = mesh;
        M3DG.geometry = geometry;
        std::cout<<"deleted\n";
        Beads[0].Pos = Saved_beadpos[0]; 
        geometry->requireVertexPositions();
        std::cout<<"relocated\n";
        psMesh = polyscope::registerSurfaceMesh("MyMesh", geometry->vertexPositions,mesh->getFaceVertexList());
        std::cout<<"re registered\n";
    }


    // scaling_factor = 1.0;
    ImGui::InputDouble("Scaling_factor", &scaling_factor);
    if(ImGui::Button("Rescale mesh ")){
        for(Vertex v: mesh->vertices()){
            geometry->inputVertexPositions[v]*=scaling_factor;
        }

        geometry->refreshQuantities();
        psMesh->updateVertexPositions(geometry->inputVertexPositions);
    }



    if(ImGui::Button("Step flow")){
        dt_sim = M3DG.integrate(Energies, Energy_constants , Sim_data, 0.0, Bead_filenames, Save_output_data);
        
        // M3DG.integrate()
    }

//     }
//     ImGui::PopStyleColor(3);
//     ImGui::SameLine();




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
    
    std::string filepath = Data["init_file"];
    int save_interval = Data["save_interval"];
    bool resize_vol = Data["resize_vol"];
    bool arcsim = Data["arcsim"];
    
    int remesh_every = 1;
    if( Data.contains("remesh_every")) remesh_every = Data["remesh_every"];
    
    
    std::cout<<"Remesh every is " << remesh_every << std::endl;
    bool Saving_last_states = Data["saving_states"];
    size_t Final_t = Data["timesteps"];

    // Here i will load the geometry 
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    

    V_bar = geometry->totalVolume();

    // We will deal with the energies now
    



    for( auto Energy : Data["Energies"]){
        Energies.push_back(Energy["Name"]);
        Constants = Energy["constants"].get<std::vector<double>>();
        std::cout<<"The constants for " << Energy["Name"]<< " are ";
        for(size_t z = 0; z < Constants.size(); z++) std::cout<<Constants[z] << " ";
        std::cout<<" \n";

        Energy_constants.push_back(Constants);
        Constants.resize(0);
    }
    
   

    Vector3 BPos;
    double radius;
    double interaction_str;
    std::string state;
    std::string interaction_mem;
    Bead PBead;
    for( auto Bead_data : Data["Beads"]){
        std::cout<<"Adding a bead\n";
        if(Bead_data.contains("gradient_order")){
            Energies.push_back(Bead_data["gradient_order"]); 
        }
        else{
        Energies.push_back("Bead");
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
        if(Bead_data["inter_str"]=="Shifted_LJ"){
            PBead.rc = radius*pow(2,1.0/6.0);
        }
        else{
            PBead.rc = 2.0*radius;
        }
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
    M3DG.recentering = Data["recentering"];
    M3DG.boundary = Data["boundary"]; 
    for( size_t i = 0 ; i< Beads.size() ; i++) M3DG.Add_bead(&Beads[i]);
    
    
    // Here i will do my alling
    // Face f = mesh->face(1);
    // M3DG.Face_sizing(f);
    // std::cout<<"Thats it\n";
    // return 1;



    // arcsim::Cloth Cloth_1;
    // arcsim::Cloth::Remeshing remeshing_params;

    if(arcsim){
        // You can make this only one function but you need to write the code for json helpers
        remeshing_params.aspect_min = Data["remesher"]["aspect_min"];
        remeshing_params.refine_angle = Data["remesher"]["refine_angle"];
        remeshing_params.refine_compression = Data["remesher"]["refine_compression"];
        remeshing_params.refine_velocity = Data["remesher"]["refine_velocity"];
        remeshing_params.size_max = Data["remesher"]["size_max"];
        remeshing_params.size_min = Data["remesher"]["size_min"];    

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
    
    
    // arcsim::Mesh remesher_mesh = translate_to_arcsim(mesh,geometry);
    // Cloth_1.mesh = remesher_mesh;
    // Cloth_1.remeshing = remeshing_params; 
    // arcsim::compute_masses(Cloth_1);
    // arcsim::compute_ws_data(Cloth_1.mesh);

    // std::cout<<"ARCSIM says \n";
    // arcsim::dynamic_remesh(Cloth_1);


    // std::cout<<" \n\n";


    // return 1;

    // delete mesh;
    // delete geometry;
    
    // std::tie(mesh_uptr, geometry_uptr) = translate_to_geometry(Cloth_1.mesh);
    // arcsim::delete_mesh(Cloth_1.mesh);
    // mesh = mesh_uptr.release();
    // geometry = geometry_uptr.release();

    M3DG.mesh = mesh;
    M3DG.geometry = geometry;

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
            stream << std::fixed << std::setprecision(4) << Beads[bead_counter].strength;
            Directory = Directory + "str_" +stream.str() +"_";

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
    
    // I need to add something that includes the bonds because then i will change that parameter(problem is)
    
    // Lets add the bonds 


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

    
    Directory = Directory+ "Nsim_" + std::to_string(Nsim)+"/";

    std::cout<<"Directory is " << Directory << " \n";
    
    
    
    
    
    std::string first_dir = Data["first_dir"];

    int status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    



    std::string basic_name=first_dir+Directory;
    
    
    status = mkdir(basic_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    std::cout<<"\nIf this number is 0 the directory was created succesfully "<< status<<"\n" ;

    std::string filename = basic_name+"Output_data.txt";


    Sim_data = std::ofstream(filename);
    Sim_data<<"T_Volume T_Area time Volume Area E_vol E_sur E_bend grad_norm backtrackstep\n";
    Sim_data.close();

    std::cout<<"Here1\n";

    // std::vector<std::string> Bead_filenames;
    // std::ofstream Bead_datas;
    

    std::cout<<"Here2\n";

    

    for(size_t i = 0 ; i < Beads.size(); i++){
        Bead_filenames.push_back(basic_name+ "Bead_"+std::to_string(i)+"_data.txt");
        Bead_datas = std::ofstream(Bead_filenames[i]);
        
        Bead_datas<<"####### This data is taken ever y" << save_interval <<" steps just like the mesh radius is " << radius<<" \n";
        Bead_datas.close();
    }

    
    std::cout<<"Here3\n";

    // Here

    Bead_filenames.push_back(basic_name+ "Simulation_timings.txt");
    Bead_datas = std::ofstream(Bead_filenames[Beads.size()]);
    Bead_datas <<"Remeshing_time Gradients Backtracking  Construction compute Solve \n";
    Bead_datas.close();


    std::cout<<"Here4\n";


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





    std::cout<<"Here5\n";

    bool seam = false;
    Cloth_1.dump_info = false;
    start = chrono::steady_clock::now();

    // Save_mesh(basic_name,-1);



    // Save_mesh(basic_name,111);
    // M3DG.Smooth_vertices();
    // geometry->refreshQuantities()
    // Save_mesh(basic_name,112);

    
    std::ofstream Dihedrals(basic_name+"dihedrals_evol.txt");
    Dihedrals << "## THE EVOLUTION OF THE DIHEDRAL DISTRIBUTION\n";
    Dihedrals.close();
    std::ofstream EdgeLengths(basic_name+"edgelengths.txt");
    EdgeLengths << "## THE EVOLUTION OF THE EDGE LENGTHS\n";
    EdgeLengths.close();


    

    // Here we start the polyscope thingy
    // We need a remesh function tho 


    std::cout<<"Here6\n";

    polyscope::init();
    TOTAL_ANGLE_DEFECT = geometry->totalAngleDefect();
    EULER_CHARACTERISTIC = geometry->eulerCharacteristic();

    polyscope::state::userCallback = functionCallback;
    geometry->requireVertexPositions();

    psMesh = polyscope::registerSurfaceMesh("MyMesh", geometry->vertexPositions,mesh->getFaceVertexList());


    std::cout<<"Here7\n";

    psMesh->setSmoothShade(true);
    psMesh->setSurfaceColor({1.0, 0.45, 0.0});
    
    polyscope::show();








    return 1;




    start_full = chrono::steady_clock::now();
    for(size_t current_t=0; current_t <= Final_t ;current_t++ ){

        // std::cout<<"Curren t t is " << current_t <<" \n";
        // if(current_t>400){
        //     save_interval = 1;
        // }

        Save_dihedrals(basic_name);
        Save_edgelengths(basic_name);
        

        start_time_control=chrono::steady_clock::now();

        if(arcsim && current_t%remesh_every==0){
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
            avg_dih2 = 0.0;
            avg_dih1 = 0.0;

            for( Edge e : mesh->edges()){ 
                dih = fabs(geometry->dihedralAngle(e.halfedge()));
                avg_dih1+=dih;
                if(dih > max_dih1) max_dih1 = dih;
                if(dih < min_dih1) min_dih1 = dih;
                }
            // std::cout<<"The average dihedral is"<< avg_dih/mesh->nEdges()<<" \n";
            // std::cout<<"The min dih is"<< min_dih << " and the max dih is " << max_dih <<" \n";
            // avg_dih =0.0;
            avg_dih1 = avg_dih1/mesh->nEdges();


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
            // std::cout<<"The change in the number of vertices is "<< n_vert_new-n_vert_old <<" \n";
            // std::cout<<"The average dihedral before was " << avg_dih1 << " and now it is " << avg_dih2 <<" \n";
            // std::cout<<"The previous min dih was" << min_dih1 << " and now it is " << min_dih2 <<" \n";
            // std::cout<<"The previous max dih was" << max_dih1 << " and now it is " << max_dih2 <<" \n\n";
             
        }
        end_time_control=chrono::steady_clock::now();
        remeshing_elapsed_time+=chrono::duration_cast<chrono::milliseconds>(end_time_control-start_time_control).count();
        Bead_datas = std::ofstream(Bead_filenames[Beads.size()],std::ios_base::app);

        Bead_datas << std::chrono::duration_cast<std::chrono::milliseconds>(end_time_control-start_time_control).count()  <<" ";
        Bead_datas.close();
                
        if(current_t%save_interval == 0){
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
                break;
            }
            // std::cout<<"No\n";
    

        }
        // std::cout<<"Redeclaring M3DG\n";
        // std::cout<<"Bead_1 position changed? "<< Bead_1.Pos << " \n";
        
        
        start_time_control = chrono::steady_clock::now();
        // std::cout<<"3\n";
        // std::cout<<"Integrating\n";
        
        // dt_sim=M3DG.integrate(TS,V_bar,nu_evol,c0,P0,KA,KB,sigma,Sim_data, time,Save_bead_data,Bead_filenames,Save_output_data,pulling);
        dt_sim = M3DG.integrate(Energies, Energy_constants , Sim_data, time, Bead_filenames, Save_output_data);
        
        if (M3DG.small_TS && current_t>Final_t*0.2) {
            std::cout << "Ending sim due to small TS \n";
            break;
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
        
        geometry->inputVertexPositions *=k;
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

        if(M3DG.small_TS) {
            std::cout<<"Ending sim due to small TS \n";
            break;
        }
        }

        if(dt_sim==-1){
            std::cout<<"Sim broke or timestep very small\n";
            std::cout<<"At timestep " << current_t << " \n";
            break;
        }
        else{
            // std::cout<<"Adding time\n";
            time+=dt_sim;
            // std::cout<<"SUccesfully\n";
        }
        if(time>10 && Beads.size()>1 ){
            Beads[1].state="froze";
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



