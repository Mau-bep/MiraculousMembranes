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

int Last_step(std::string base_dir){

    DIR *dh;
    struct dirent * contents;
    dh = opendir(base_dir.c_str());

    int index_underscore;
    int index_dot;
    string max_name;
    size_t max_index=0;
    size_t current_index=0;
    
    if ( !dh )
    {
        std::cout << "The given directory is not found";
        return -1;
    }
    while ( ( contents = readdir ( dh ) ) != NULL )
    {
        std::string name = contents->d_name;
        // std::cout<<name.size()<<endl;
        
        if(name.size()>7 && name.find('t')>10){
        index_underscore=name.find('_');
        index_dot=name.find('.');
        string str_index = name.substr(index_underscore+1,index_dot-index_underscore-1 );
        // std::cout << str_index << endl;
        // I need to avoid the final state
        current_index=stoi(str_index);
        // std::cout<<current_index<<endl;
        if(current_index>max_index){
            max_index=current_index;
            max_name= name;
        }

        // std::cout << name.substr(2) << endl;
        }
    }
    std::cout<< "The biggest index is "<< max_index << endl;
    closedir ( dh );
    



    return max_index;
}

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



int main(int argc, char** argv) {



    //

    double arr[] = { 1.0, 2.0, 3.0 }; 
    int n = sizeof(arr) / sizeof(arr[0]); 
  
    vector<double> radius(arr, arr + n); 

    double arr_2[] = {100.0,1000.0,10000.0,100000.0};
    n=sizeof(arr_2) / sizeof(arr_2[0]);
    Vector<double> KAs(arr_2,arr_2+n);
    
    KB=0.1;
    
    double arr_3[] = {0.1,0.5,1.0,5.0,10.0,15.0,15.0,20.0,40.0,50.0};
    n=sizeof(arr_3) / sizeof(arr_3[0]);
    
    Vector<double> strengths(arr_3,arr_3+n);
    
    nu=1.0;
    c0=0.0;

    int Init_cond=2;
    int Nsim=1;

    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    
    bool Save_data=false;
    double time;
    double dt_sim;
    double rad;
    double Interaction_str;
    int Last_ts;
    int num_steps;

    std::cout<< "Current path is " << argv[0]<<"\n";




    // Mesh related data 
    std::string filepath;
    filepath = "../../../input/sphere.obj";
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    std::cout<<"The mesh is correctly loaded\n";
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    trgt_len=geometry->meanEdgeLength();
    V_bar=geometry->totalVolume();
    ORIG_VPOS = geometry->inputVertexPositions;
    CoM = geometry->centerOfMass();
    
    std::stringstream nustream;
    std::stringstream c0stream;
    std::stringstream KAstream;
    std::stringstream KBstream;
    
    
    std::stringstream Curv_adapstream;
    std::stringstream Min_rel_lengthstream;
    std::stringstream radiusstream;
    std::stringstream Interactionstrstream;

    std::string basic_name;
    std::string first_dir;
    std::string filename;

    nustream << std::fixed << std::setprecision(3) << nu;
    c0stream << std::fixed << std::setprecision(3) << c0;
    KBstream << std::fixed << std::setprecision(6) << KB;
    Curv_adapstream << std::fixed << std::setprecision(2) << Curv_adap;
    Min_rel_lengthstream << std::fixed << std::setprecision(2) <<Min_rel_length;

    first_dir="../Results/Mem3DG_Bead_Reciprocal_finemesh/";
    std::string basic_name;
    for (int r_it = 0 ; r_it < radius.size(); r_it++){
        rad = radius[r_it];    
        radiusstream << std::fixed << std::setprecision(3) << rad;
    
        for(int KA_it = 0 ; KA_it<KAs.size(); KA_it++){
            KA = KAs[KA_it];
            KAstream << std::fixed << std::setprecision(3) << KA;
    
            for(int str_it = 0; str_it<strengths.size(); str_it++){
                Interaction_str=strengths[str_it];
                Interactionstrstream << std::fixed << std::setprecision(6) << Interaction_str;
                // Those are all the manual parameters
                
                // We now load the directory
                basic_name=first_dir+"nu_"+nustream.str()+"_radius_"+radiusstream.str()+"_curvadap_"+Curv_adapstream.str()+"_minrel_"+Min_rel_lengthstream.str()+"_KA_"+KAstream.str()+"_KB_"+KBstream.str()+"_strength_"+Interactionstrstream.str()+"_Init_cond_"+std::to_string(Init_cond)+"_Nsim_"+std::to_string(Nsim)+"/";

                // 
                // This is the last step to consider

                Last_ts = Last_step(basic_name);
                num_steps=Last_ts/500+1;
                Vector<Vector3> Bead_pos(num_steps);
                filename = basic_name+ "Bead_data.txt";
                // I now need to read the bead data 
                std::ifstream Bead_data(filename);
                if(!Bead_data.is_open()){
                    cerr<<"Error opening the file!"<< endl;
                    return 1;
                }
                string line;
                while (std::getline(Bead_data,line)){
                    std::vector<std::string> splitted = split(line," ");
                    // I want
                    
                }
                Bead_data.close()

                





            }
        }
    }

  
    


    std::string filename_basic = basic_name+"Output_data.txt";

    std::ofstream Sim_data;
    

    std::string filename2 = basic_name + "Bead_data.txt";
    std::ofstream Bead_data(filename2);
    bool Save_bead_data=false;
        


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



    std::string filename;
    std::ofstream Gradient_data_bead;

    std::ofstream Gradient_data_bead_dx;
    


    filename = basic_name+"Bead_Gradient_evaluation_dx_"+std::to_string(0) + ".txt";







    delete mesh;
    delete geometry;

    return EXIT_SUCCESS;
}



