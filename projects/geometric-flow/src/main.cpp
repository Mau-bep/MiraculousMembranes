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


// #include "Mem-3dg.h"
#include "Mem-3dg_implicit.h"

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
polyscope::SurfaceMesh* psMesh;

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
double TS=0.0001;

double Curv_adap=0.1;
double Min_rel_length=0.5;
double trgt_len;
double avg_remeshing;
bool edge_length_adj;
VertexData<Vector3> ORIG_VPOS; // original vertex positions
Vector3 CoM;                   // original center of mass

IMem3DG IM3DG;

// RemeshBoundaryCondition defaultRemeshOptions;


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
    // c0=std::stod(argv[2]);
    // KA=std::stod(argv[3]);
    // KB=std::stod(argv[4]);
    // nu=1.0;
    KA=10.0;
    KB=0.005;
    double init_step=0.0;
    // I will do it so i can give this values
 
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    

    
    TS=pow(10,-4);


    std::cout<< "Current path is " << argv[0];

    // std::string filepath = "../../../input/deformed_sphere_2.obj";
    std::string filepath = "../../../input/Simple_cil_regular.obj";
    // Load mesh
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    
    trgt_len=geometry->meanEdgeLength();
    V_bar=geometry->totalVolume();
    
    // Initialize operators.
  
    

    ORIG_VPOS = geometry->inputVertexPositions;
    CoM = geometry->centerOfMass();
    
    // MCF = MeanCurvatureFlow(mesh, geometry);
    // ModMCF = ModifiedMeanCurvatureFlow(mesh, geometry);
    // NF =NormalFlow(mesh, geometry);
    // GCF = GaussCurvatureFlow(mesh, geometry);
    // WF = WillmoreFlow(mesh,geometry);
    // WF2 = WillmoreFlow2(mesh,geometry);
    // WFS = WillmoreFlowScho(mesh,geometry);
 
    IM3DG = IMem3DG(mesh,geometry);

    // Add visualization options.
    // polyscope::screenshot("./This_filename_is_nice_right.jpg",true);

    // size_t counter=0;
    
    std::stringstream nustream;
    std::stringstream c0stream;
    std::stringstream KAstream;
    std::stringstream KBstream;
    
    
    // std::stringstream H0stream;
    // std::stringstream kappastream;
    // std::stringstream sigmastream;
    std::stringstream Curv_adapstream;
    std::stringstream Min_rel_lengthstream;


    nustream << std::fixed << std::setprecision(3) << nu;
    c0stream << std::fixed << std::setprecision(3) << c0;
    KAstream << std::fixed << std::setprecision(3) << KA;
    KBstream << std::fixed << std::setprecision(6) << KB;



    Curv_adapstream << std::fixed << std::setprecision(2) << Curv_adap;
    Min_rel_lengthstream << std::fixed << std::setprecision(2) <<Min_rel_length;
    
    std::string first_dir="../Results/Implicit/";
    int status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"If this name is 0 the directory was created succesfully "<< status ;

    // first_dir="../Results/Mem3DG_IMG_serial/Curv_adap_"+Curv_adapstream.str()+"Min_rel_length_"+Min_rel_lengthstream.str();
    // status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"\nIf this name is 0 the directory was created succesfully "<< status ;
    
    std::string basic_name =first_dir+"nu_"+nustream.str()+"_c0_"+c0stream.str()+"_KA_"+KAstream.str()+"_KB_"+KBstream.str()+"/";
    status = mkdir(basic_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    std::cout<<"\nIf this number is 0 the directory was created succesfully "<< status<<"\n" ;

    std::string filename = basic_name+"Output_data.txt";

    std::ofstream Sim_data;

    // Sim_data=std::ofstream(filename);
    // Sim_data<<"T_Volume T_Area time Volume Area E_vol E_sur E_bend grad_norm backtrackstep\n";
    
    // Here i want to run my video
    size_t n_vert;
    size_t n_vert_old;
    size_t n_vert_new;
    double Volume;
    double Area;
    double Area_bar;
    double nu_obs;
    double nu_evol=0.0;
    double nu_0;
    double c_null;
    size_t counter=0;
    double time=0.0;
    double dt_sim=0.0;
    double nu_step;

    start = chrono::steady_clock::now();
    for(size_t current_t=0;current_t<10000;current_t++ ){
        // for(size_t non_used_var=0;non_used_var<100;)
        // MemF.integrate(TS,sigma,kappa,H0,P,V0);
        if(true){
        // if(current_t%10==0 ){
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
        Options.remesh_list=false;
        MutationManager Mutation_manager(*mesh,*geometry);
        remesh(*mesh,*geometry,Mutation_manager,Options);
        n_vert_new=mesh->nVertices();
        counter=counter+1; 
        }
        }

        
           if(current_t%100==0){
            Save_mesh(basic_name,current_t);

        }
        
        
        if(current_t%100==0) {
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
           
            // c_null=2*H0 *pow(Area/(4*PI),0.5);
            if(current_t==0){
                nu_0=nu_obs;
                // nu_= (nu-nu_0)/50;
            }
            std::cout<<" THe simulation time is "<< time <<"\n";
            // std::cout<< "The spontaneous curvature is " << H0<< "\n";
            // std::cout << "The system time is " << M3DG.system_time <<"\n\n";
            
            std::cout<<"Ten iterations took "<<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" miliseconds\n\n\n";


            // polyscope::screenshot(basic_name+std::to_string(current_t)+".jpg",true);
            start = chrono::steady_clock::now();

        }



        // I just want to change this part.
        double V = geometry->totalVolume();
        double D_P=-P0*(V-V_bar)/V_bar/V_bar;
        // VertexData<Vector3> Surf_grad=M3DG.SurfaceGrad();
        // VertexData<Vector3> Volume_grad=M3DG.OsmoticPressure(D_P);

        // Now i need to integrate 

        // geometry->inputVertexPositions+=1e-3*(KA*Surf_grad+Volume_grad);
        Area_bar=4*PI*pow(3*V_bar/(4*PI*nu),2.0/3.0);

        // nu_evol= (current_t-init_step) <Area_evol_steps ? nu_0 + (nu-nu_0)*(current_t-init_step)/Area_evol_steps : nu; 
        // if(current_t>=2*Area_evol_steps){
        // KB_evol= (current_t-init_step) <3*Area_evol_steps ? 0.1 + (KB-0.1)*(current_t-init_step-2*Area_evol_steps)/Area_evol_steps : KB; 
        // }
        // std::cout<<"Lets do the integration\n";
        IM3DG.integrate(TS,nu,V_bar,P0,KA,KB);
        time+=1e-3;
        // nu_evol= time<50 ? nu_0 + (nu-nu_0)*time/50 : nu; 
        // nu_evol=nu;
        // bool Save_output_data=false;
        // dt_sim=M3DG.integrate(TS,V_bar,nu_evol,c0,P0,KA,KB,Kd,Sim_data,time,Save_output_data);
        // if(dt_sim==-1){
            // std::cout<<"Sim broke\n";
            // break;
        // }
        // else{
            // time+=dt_sim;
            
        // }





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



