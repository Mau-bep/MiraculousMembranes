// GEOMETRIC FLOW

// #include <stdlib.h>
#include <unistd.h>


#include <omp.h>


#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <string>

# include <dirent.h>

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

double Curv_adap=0.0;
double Min_rel_length=0.1;
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
    std::cout<<base_dir<<"\n";
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
        // std::cout<<"filename is "<<name<<" \n";
        // std::cout<<"Find t is"<<name.find('t')<<endl;
        
        if(name.size()>7 && name.find('t')>100){
        index_underscore=name.find('_');
        index_dot=name.find('.');
        string str_index = name.substr(index_underscore+1,index_dot-index_underscore-1 );
        // std::cout << str_index << endl;
        
        // I need to avoid the final state
        // std::cout<<"str_index is" << str_index <<"\n";
        if(str_index[0]=='f'){
            continue;
        }
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


double Interaction_E_from_Output(std::string  filename){

    
    std::ifstream Output_data(filename);
    if(!Output_data.is_open()){
        cerr<<"Error opening the file!"<< endl;
        return 1;
    }

    string line;
    int counter=0;
    double E_I = 0.0;

    while (std::getline(Output_data,line)){
        std::vector<std::string> splitted = split(line," ");
        if(line[0]=='#'|| line[0] =='V') {
            std::cout<<"First line\n";
            continue;
        }
        if(splitted.size()>10){
            // std::cout<<splitted[6] <<" \n";
        E_I  = std::stod(splitted[6]);    

        }
        
        // I want the first three
    }
    Output_data.close();





    return  E_I;
}


int main(int argc, char** argv) {



    //

    // double arr[] = { 1.0, 2.0, 3.0 };
    double arr[] = { 0.5}; 
    int n = sizeof(arr) / sizeof(arr[0]); 
  
    vector<double> radius(arr,arr+n);
    // radius.push_back(1.0);
    // radius.push_back(3.0);

    // double arr_2[] = {100.0};

    double arr_2[] = {45,50};
    n=sizeof(arr_2) / sizeof(arr_2[0]);
    vector<double> KAs(arr_2,arr_2+n);
    
    
    
    // double arr_3[] = {10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0 };
    double arr_3[] = {20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0,90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0};
    
    n=sizeof(arr_3) / sizeof(arr_3[0]);
    
    vector<double> strengths(arr_3,arr_3+n);
    
    KB=1.0;
    
    double arr_4[] = { 10.0};
    n = sizeof(arr_4) / sizeof(arr_4[0]);

    vector<double> KBs(arr_4,arr_4+n);

    nu=1.0;
    c0=0.0;

    double trgt_vol = 4*3.1415/3.0;
    double trgt_area = 4*3.1415;


    int Init_cond=2;
    int Nsim=2;

    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    
    bool Save_data=false;
    double time;
    double dt_sim;
    double rad;
    double Interaction_str;
    int Last_ts;
    int num_steps;
    int delta_steps=1000;
    Vector3 Bead_pos;
    std::cout<< "Current path is " << argv[0]<<"\n";
    bool check_coverage = true;
    bool check_forces = false;


    for ( size_t KB_it = 0; KB_it < KBs.size(); KB_it++){
        KB = KBs[KB_it];
    



    std::stringstream KAstream;
    std::stringstream KBstream;
    
    
    std::stringstream radiusstream;
    std::stringstream Interactionstrstream;

    std::string basic_name;
    std::string first_dir;
    std::string filename;

   
    KBstream << std::fixed << std::setprecision(4) << KB;

    first_dir="../Results/Particle_wrapping_on_plane_phase_space_sept/";
    
    filename = first_dir + "Coverage_final.txt" ;
    std::ofstream Coverage_final(filename,std::ios_base::app);
    if(KB_it == 0) Coverage_final<<"# # # Coverage data \n";
    Coverage_final.close(); 
    int counter=0;

    double avg_rmin =0.0;
    int avg_rmin_counter =0;
    for (int r_it = 0 ; r_it < radius.size(); r_it++){
        std::cout<<"Hre\n";
        rad = radius[r_it];    
        radiusstream.str(std::string());
        radiusstream.clear();

        radiusstream << std::fixed << std::setprecision(4) << rad;
    
        for(size_t KA_it = 0 ; KA_it<KAs.size(); KA_it++){
            KA = KAs[KA_it];
            KAstream.str(std::string());
            KAstream.clear();
            
            KAstream << std::fixed << std::setprecision(4) << KA;
    
            for(int str_it = 0; str_it<strengths.size(); str_it++){
                // std::cout<<"How many?\n";
                Interaction_str=strengths[str_it];
                Interactionstrstream.str(std::string());
                Interactionstrstream.clear();
                Interactionstrstream << std::fixed << std::setprecision(4) << Interaction_str;
                // Those are all the manual parameters
                
                // We now load the directory
                basic_name = first_dir +"Surface_tension_"+KAstream.str() + "_Bending_"+KBstream.str()+"_Edge_reg_1.0000_"+"Bead_radius_" + radiusstream.str() + "_Frenkel_Normal_nopush_str_" + Interactionstrstream.str()+ "_Switch_Newton_Nsim_" + std::to_string(Nsim) + "/";
            
                // 
                // This is the last step to consider

                Last_ts = Last_step(basic_name);
                
                num_steps=Last_ts/delta_steps+1;
                Vector<Vector3> Bead_pos(num_steps);
                filename = basic_name+ "Bead_0_data.txt";
                // I now need to read the bead data 
                std::ifstream Bead_data(filename);
                if(!Bead_data.is_open()){
                    cerr<<"Error opening the file!"<< endl;
                    return 1;
                }
                string line;
                counter=0;
                while (std::getline(Bead_data,line)){
                    if(counter>=num_steps) break;
                    std::vector<std::string> splitted = split(line," ");
                    // i NEED TO MAKE SURE THERE ARE NUMBERS HERE
                    if(line[0]=='#') {
                        std::cout<<"First line\n";
                        continue;
                    }
                    if(splitted.size()>5){
                    
                    Bead_pos[counter].x=std::stod(splitted[0]);
                    Bead_pos[counter].y=std::stod(splitted[1]);
                    Bead_pos[counter].z=std::stod(splitted[2]);
                    counter+=1;
                    }
                    
                    // I want the first three
                }
                Bead_data.close();

                // I have the bead positions and the directories of the objs, Lesgo

                std::cout<<"Bead data read\n";
                std::cout<<"Counter is"<< counter <<"\n";
                if(counter == 1) continue;
                // filename = basic_name + "Coverage_evol.txt";
                // std::ofstream Coverage(filename); 
                double covered_area=0;
                double relative_coverage=0;
                double rmin = 12;
                double x_max_mem = 0.0;
                bool first=true;


                double E_I = Interaction_E_from_Output(basic_name + "Output_data.txt");


                for(int step = counter-1 ; step<counter;step++){
                    
                // Mesh related data 
                std::string filepath;
                filepath = basic_name+"membrane_"+std::to_string(step*delta_steps)+".obj";

                // std::cout<<"Reading mesh\n";
                std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
                
                mesh = mesh_uptr.release();
                geometry = geometry_uptr.release();
                trgt_len=geometry->meanEdgeLength();
                V_bar=geometry->totalVolume();
                ORIG_VPOS = geometry->inputVertexPositions;
                CoM = geometry->centerOfMass();
                
                // The mesh and the bead position are loaded, time to measure

                Vector3 Vert_pos;
                Vector3 Bead_current;
                Vector3 Normal;
                Vector3 rij;
                double r_dist;

                double mag_Vol;
                double mag_Area;
                double mag_Bend;
                double mag_Bead;
                // I want to check
                
                
                // I can check the rdist
                Bead_current=Bead_pos[step];
                std::cout<<Bead_current<<"This is the bead position\n";
                // std::cout<<"The bead position is"<< Bead_current<<"\n";
                int touching_count=0;
                filename = first_dir + "Radius_distribution_strength_"+to_string(Interaction_str)+".txt";
                std::ofstream R_dist(filename);
                
                filename = first_dir + "Touching_strength_"+to_string(Interaction_str)+".txt";
                std::ofstream Touching_data(filename);

                covered_area=0.0;
                // std::cout<<"Iterating over vertices\n";

                for(int v = 0 ; v < mesh->nVertices() ; v++){
                    Vert_pos=geometry->inputVertexPositions[v];
                    // if(Vert_pos.x>x_max_mem) x_max_mem = Vert_pos.x;
                    rij = (Bead_current-Vert_pos);
                    r_dist = rij.norm();
                    if(r_dist<rmin) rmin = r_dist;

                }
                // x_max_mem+=1.0;
                // Vector3 Second_bead({x_max_mem,0.0,0.0});
                // double second_dist=0.0;
                for(int v =0; v<mesh->nVertices(); v++){
                    // So i have the radius and the 
                    Vert_pos=geometry->inputVertexPositions[v];
                    rij = (Bead_current-Vert_pos);
                    r_dist = rij.norm();
                    // second_dist = (Vert_pos-Second_bead).norm();

                    Normal = geometry->vertexNormalAreaWeighted(mesh->vertex(v));
                    Normal = Normal.unit();
                    Vector3 x_dir;
                    
                    if(first) {
                        R_dist<<r_dist;
                        first = false;
                        }
                    else R_dist<<" "<<r_dist;
                    rij = rij.unit();
                    // I want the distribution saved so 
                    
                
                    if(check_coverage){
                    // Now i need to do my part

                    if( (dot(rij,Normal)>0.0 && r_dist<rad*1.2 ) ){
                        

                        
                        // if(r_dist<rmin) rmin = r_dist;
                        Touching_data<<Vert_pos.x <<" "<< Vert_pos.y <<" "<<Vert_pos.z<<"\n";
                        touching_count+=1;
                        covered_area+=geometry->barycentricDualArea(mesh->vertex(v));
                        // std::cout<<"Dot"<< dot(rij,Normal)<<"\t";

                    }
                    

                    


                    }
                    // if(check_forces){
                    //     if(r_dist<rad*1.1 && dot(rij,Normal)>0){
                    //         // This are the vertices interacting with the bead
                            

                    //     }
                    // // There are multiple ways to go about this
                    // // I will check 




                    // }
                
                    // 
                
                
                }
                std::cout<<"\n";
                R_dist.close();
                Touching_data.close();
                // std::cout<<"The amount touching is"<< touching_count<<" \n";
                
                relative_coverage=covered_area/(4*PI*(rad*1.155)*(rad*1.155));
                avg_rmin += rmin;
                avg_rmin_counter+=1;
                // Coverage<<relative_coverage<<"\n";
                // R_dist.close();
                

                }
                // Coverage.close();
                filename = first_dir + "Coverage_final.txt" ;
                Coverage_final =std::ofstream(filename,std::ios_base::app); 
    
                Coverage_final<< rad<<" "<< KB << " "<< E_I<<" "<< relative_coverage<<" "<< Interaction_str <<" "<< KA <<"\n";
                Coverage_final.close();

                // R_dist.close()
            }
        }
    }
    std::cout<<"The average rmin overall is "<< avg_rmin/avg_rmin_counter <<"\n";

    }
    // Coverage_final.close();
    






    delete mesh;
    delete geometry;

    return EXIT_SUCCESS;
}



