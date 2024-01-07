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


#include <chrono>


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"


#include "Mem-3dg.h"

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
float KB=0.0001;
float Kd=0.0;
double TS=0.001;

size_t Area_evol_steps=100000;

double Curv_adap=0.1;
double Min_rel_length=0.5;
double trgt_len;
double avg_remeshing;
// bool edge_length_adj;
bool Save=false;
VertexData<Vector3> ORIG_VPOS; // original vertex positions
Vector3 CoM;                   // original center of mass

Mem3DG M3DG;


void showSelected() {
    // pass
}


void Save_mesh(std::string basic_name, size_t current_t) {
   // Build member variables: mesh, geometry
    Vector3 Pos;
    std::ofstream o(basic_name+"Membrane_"+std::to_string(current_t)+".obj");
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


int main(int argc, char** argv) {

    int Nsim;
    nu=std::stod(argv[1]);
    int Init_cond = std::stoi(argv[2]);
    Nsim = std::stoi(argv[3]);
    KB=std::stod(argv[4]);

    int init_step=0;
    c0=0.0;
    KA=10.0;
    // KB=0.001;

    // I will do it so i can give this values
 
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    


    TS=pow(10,-3);


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
    
    std::string first_dir="../Results/Mem3DG_Cell_Shape/";
    int status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"If this name is 0 the directory was created succesfully "<< status ;

    std::string basic_name=first_dir+"nu_"+nustream.str()+"_c0_"+c0stream.str()+"_KA_"+KAstream.str()+"_KB_"+KBstream.str()+"_init_cond_"+std::to_string(Init_cond)+"_Nsim_"+std::to_string(Nsim)+"/";
    status = mkdir(basic_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    std::cout<<"\nIf this number is 0 the directory was created succesfully "<< status<<"\n" ;


    // Here we will decide if we want to import the previous mesh or we want to use the default initial conditions



    bool Continue_sim=true;

    std::cout<< "Current path is " << argv[0];
    std::string filepath;
    

    if(Continue_sim){


    init_step = Last_step(basic_name);
    
    filepath= basic_name+"Membrane_"+std::to_string(init_step)+".obj";

    }
    else{

    // std::string filepath = "../../../input/sphere.obj";
    if(Init_cond==1){
        filepath = "../../../input/Pill_regular.obj";
    
    }
    if(Init_cond==2){
        filepath = "../../../input/bloodcell.obj";
    }
    if(Init_cond==3){
        filepath = "../../../input/Pushed_sphere_v_regular.obj";
    }
    }
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    
    
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    
    trgt_len=geometry->meanEdgeLength();
    V_bar=geometry->totalVolume();

    ORIG_VPOS = geometry->inputVertexPositions;
    CoM = geometry->centerOfMass();
    
    M3DG = Mem3DG(mesh,geometry);
    M3DG.pulling=false;
    
    std::string filename = basic_name+"Output_data.txt";


    std::ofstream Sim_data;
    if(Continue_sim){
        Sim_data=std::ofstream(filename,std::ios_base::app);
    }
    else{
        Sim_data=std::ofstream(filename);
        Sim_data<<"T_Volume T_Area time Volume Area E_vol E_sur E_bend grad_norm backtrackstep\n";
    
    }
        
    size_t n_vert;
    size_t n_vert_old;
    size_t n_vert_new;
    double Volume;
    double Area;
    double Area_bar;
    double nu_obs;
    double nu_evol;
    double nu_0;
    double c_null;
    size_t counter=0;
    double time=0.0;
    double dt_sim=0.0;
    // EdgeData<int> No_remesh(*mesh,0);

    start = chrono::steady_clock::now();
    for(size_t current_t=init_step;current_t<init_step + 8000000;current_t++ ){
        // for(size_t non_used_var=0;non_used_var<100;)
        // MemF.integrate(TS,sigma,kappa,H0,P,V0);
        if(true){
        n_vert_old=0;
        n_vert_new=mesh->nVertices();

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
        // 
        Options.remesh_list=false;
        // Options.No_remesh_list=No_remesh;
        MutationManager Mutation_manager(*mesh,*geometry);
        remesh(*mesh,*geometry,Mutation_manager,Options);
        n_vert_new=mesh->nVertices();
        counter=counter+1; 
        }
        }

        
        if(current_t%2500==0){
            Save_mesh(basic_name,current_t);

        }
        
        if(current_t%1000==0) {
            Save=true;
            end=chrono::steady_clock::now();
            n_vert=mesh->nVertices();
            std::cout<< "THe number of vertices is "<< n_vert <<"\n";    
            std::cout<< "the avg edge lenghth is " <<geometry->meanEdgeLength()<<"\n";
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

            // std::cout<< "The spontaneous curvature is " << H0<< "\n";
            std::cout << "The system time is " << M3DG.system_time <<"\n";
            
            std::cout<<"A thousand iterations took "<<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" miliseconds\n\n\n";

            start = chrono::steady_clock::now();

        }


        // I said i should add something that considers the change in area per timestep
        // Cause i dont want to add more area that i can handle per ts
        // I can estimate the target area

        Area_bar=4*PI*pow(3*V_bar/(4*PI*nu),2.0/3.0);
        
        // and the idea is that i can only add a certain percent of area , i think linearly is the best
        // I think a 1% pertime step should be alright
        



        nu_evol= (current_t-init_step) <Area_evol_steps ? nu_0 + (nu-nu_0)*(current_t-init_step)/Area_evol_steps : nu; 
        
        dt_sim=M3DG.integrate(TS,V_bar,nu_evol,c0,P0,KA,KB,Kd,Sim_data,time,Save);
        if(Save==true)
        {
            Save=false;
        }
        if(dt_sim==-1){
            std::cout<<"Sim broke or timestep very small\n";
            break;
        }
        else{
            time+=dt_sim;
            
        }





    }
    Sim_data.close();


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



