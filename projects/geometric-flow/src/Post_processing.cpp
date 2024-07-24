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

    // double arr[] = { 1.0, 2.0, 3.0 };
    double arr[] = { 1.0}; 
    int n = sizeof(arr) / sizeof(arr[0]); 
  
    vector<double> radius(arr,arr+n);
    // radius.push_back(1.0);
    // radius.push_back(3.0);

    // double arr_2[] = {100.0};

    double arr_2[] = {100000.0};
    n=sizeof(arr_2) / sizeof(arr_2[0]);
    vector<double> KAs(arr_2,arr_2+n);
    
    KB=1.0;
    
    double arr_3[] = {1.0,5.0,10.0,20.0,30.0,40.0,50.0, 60.0, 70.0, 80.0, 90.0,100.0, 110.0, 120.0, 130.0, 140.0, 150.0};
    // double arr_3[] = {10.0};
    
    n=sizeof(arr_3) / sizeof(arr_3[0]);
    
    vector<double> strengths(arr_3,arr_3+n);
    
    nu=1.0;
    c0=0.0;

    double trgt_vol = 566.44;
    double trgt_area = 331.09;


    int Init_cond=1;
    int Nsim=777;

    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    
    bool Save_data=false;
    double time;
    double dt_sim;
    double rad;
    double Interaction_str;
    int Last_ts;
    int num_steps;
    Vector3 Bead_pos;
    std::cout<< "Current path is " << argv[0]<<"\n";
    bool check_coverage = true;
    bool check_forces = false;



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
    Min_rel_lengthstream << std::fixed << std::setprecision(4) <<Min_rel_length;

    first_dir="../Results/Mem3DG_Bead_Reciprocal_arcsim/";
    filename = first_dir + "Coverage_final.txt" ;
    std::ofstream Coverage_final(filename,std::ios_base::app);
    Coverage_final<<"# # # Coverage data \n";
    Coverage_final.close(); 
    int counter=0;
    for (int r_it = 0 ; r_it < radius.size(); r_it++){
        std::cout<<"Hre\n";
        rad = radius[r_it];    
        radiusstream.str(std::string());
        radiusstream.clear();

        radiusstream << std::fixed << std::setprecision(3) << rad;
    
        for(int KA_it = 0 ; KA_it<KAs.size(); KA_it++){
            KA = KAs[KA_it];
            KAstream.str(std::string());
            KAstream.clear();
            
            KAstream << std::fixed << std::setprecision(3) << KA;
    
            for(int str_it = 0; str_it<strengths.size(); str_it++){
                // std::cout<<"How many?\n";
                Interaction_str=strengths[str_it];
                Interactionstrstream.str(std::string());
                Interactionstrstream.clear();
                Interactionstrstream << std::fixed << std::setprecision(6) << Interaction_str;
                // Those are all the manual parameters
                
                // We now load the directory
                basic_name=first_dir+"nu_"+nustream.str()+"_radius_"+radiusstream.str()+"_KA_"+KAstream.str()+"_KB_"+KBstream.str()+"_strength_"+Interactionstrstream.str()+"_Init_cond_"+std::to_string(Init_cond)+"_Nsim_"+std::to_string(Nsim)+"/";

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

                // filename = basic_name + "Coverage_evol.txt";
                // std::ofstream Coverage(filename); 
                double covered_area=0;
                double relative_coverage=0;
                double rmin = 12;
                for(int step = counter-1 ; step<counter;step++){
                    
                // Mesh related data 
                std::string filepath;
                filepath = basic_name+"membrane_"+std::to_string(step*500)+".obj";

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
                // std::cout<<"The bead position is"<< Bead_current<<"\n";
                int touching_count=0;
                filename = basic_name + "Radius_distribution_step_"+to_string(str_it)+".txt";
                std::ofstream R_dist(filename);
                
                // filename = basic_name + "Touching_step_"+to_string(step)+".txt";
                // std::ofstream Touching_data(filename);
                covered_area=0.0;
                // std::cout<<"Iterating over vertices\n";

                for(int v = 0 ; v < mesh->nVertices() ; v++){
                    Vert_pos=geometry->inputVertexPositions[v];
                    rij = (Bead_current-Vert_pos);
                    r_dist = rij.norm();
                    if(r_dist<rmin) rmin = r_dist;
                }

                for(int v =0; v<mesh->nVertices(); v++){
                    // So i have the radius and the 
                    Vert_pos=geometry->inputVertexPositions[v];
                    rij = (Bead_current-Vert_pos);
                    r_dist = rij.norm();
                    Normal=geometry->vertexNormalAreaWeighted(mesh->vertex(v));

                    // I want the distribution saved so 
                    if(v==0) R_dist<<r_dist;
                    else R_dist<<" "<<r_dist;
                
                    if(check_coverage){
                    // Now i need to do my part

                    if(r_dist<rmin+0.3 && dot(rij,Normal)>0){
                        // if(r_dist<rmin) rmin = r_dist;
                        // Touching_data<<Vert_pos.x <<" "<< Vert_pos.y <<" "<<Vert_pos.z<<"\n";
                        // touching_count+=1;
                        covered_area+=geometry->barycentricDualArea(mesh->vertex(v));

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
                R_dist.close();
                // std::cout<<"The amount touching is"<< touching_count<<" \n";
                relative_coverage=covered_area/(4*PI*(rmin)*(rmin));
                // Coverage<<relative_coverage<<"\n";
                // R_dist.close();
                

                }
                // Coverage.close();
                filename = first_dir + "Coverage_final.txt" ;
                Coverage_final =std::ofstream(filename,std::ios_base::app); 
    
                Coverage_final<< rad<<" "<< KA << " "<< Interaction_str<<" "<< relative_coverage<<"\n";
                Coverage_final.close();

                // R_dist.close()
            }
        }
    }

    // Coverage_final.close();
    






    delete mesh;
    delete geometry;

    return EXIT_SUCCESS;
}



