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


// #include "Mem-3dg.h"
#include "Mem-3dg_implicit.h"
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
float KA=10;
float KB=10.0;
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

int main(int argc, char** argv) {

    int Nsim;
    // nu=std::stod(argv[1]);
    // c0=std::stod(argv[2]);
    KA=std::stod(argv[1]);
    KB=std::stod(argv[2]);
    nu=0.650;
    // KA=50;
    // KB=10.0;
    Nsim = std::stoi(argv[3]);
    double init_step=0.0;
    // I will do it so i can give this values
 
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    

    bool Save_output_data=true;
    TS=pow(10,-4);



    arcsim::Cloth Cloth_1;
    arcsim::Cloth::Remeshing remeshing_params;

    remeshing_params.aspect_min = 0.4;
    remeshing_params.refine_angle = 0.6;
    remeshing_params.refine_compression = 1.0;
    remeshing_params.refine_velocity = 1.0;
    remeshing_params.size_max = 0.1;
    remeshing_params.size_min = 0.001;   


    std::cout<< "Current path is " << argv[0];

    // std::string filepath = "../../../input/Pushed_sphere_v_regular.obj";
    std::string filepath = "../../../input/Simple_cil_regular.obj";
    // std::string filepath = "../../../input/bloodcell.obj";

    // std::string filepath = "../Results/Mem3DG_Cell_Shape_KB_evol_flip/nu_0.625_c0_0.000_KA_10.000_KB_0.010000_init_cond_2_Nsim_11/Final_state.obj";
    // Load mesh
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    
    trgt_len=geometry->meanEdgeLength();
    V_bar=geometry->totalVolume();
    
    // Initialize operators.
  
    
    ORIG_VPOS = geometry->inputVertexPositions;
    CoM = geometry->centerOfMass();
    
    
    
    arcsim::Mesh remesher_mesh = translate_to_arcsim(mesh,geometry);
    Cloth_1.mesh = remesher_mesh;
    Cloth_1.remeshing = remeshing_params; 
    arcsim::compute_masses(Cloth_1);
    arcsim::compute_ws_data(Cloth_1.mesh);
    arcsim::dynamic_remesh(Cloth_1);
    delete mesh;
    delete geometry;

    std::tie(mesh_uptr, geometry_uptr) = translate_to_geometry(Cloth_1.mesh);
    arcsim::delete_mesh(Cloth_1.mesh);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();


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
    
    std::string first_dir="../Results/Sobolev/";
    int status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"If this name is 0 the directory was created succesfully "<< status ;

    // first_dir="../Results/Mem3DG_IMG_serial/Curv_adap_"+Curv_adapstream.str()+"Min_rel_length_"+Min_rel_lengthstream.str();
    // status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"\nIf this name is 0 the directory was created succesfully "<< status ;
    
    std::string basic_name =first_dir+"nu_"+nustream.str()+"_c0_"+c0stream.str()+"_KA_"+KAstream.str()+"_KB_"+KBstream.str()+"_Nsim"+std::to_string(Nsim)+"/";
    status = mkdir(basic_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    std::cout<<"\nIf this number is 0 the directory was created succesfully "<< status<<"\n" ;

    std::string filename = basic_name+"Output_data.txt";

    std::ofstream Sim_data(filename);

    // Sim_data=std::ofstream(filename);
    Sim_data<<"T_Volume T_Area time Volume Area E_bend grad_norm backtrackstep\n";
    Sim_data.close();
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

    geometry->normalize(Vector3({0.0,0.0,0.0}),false);
    geometry->refreshQuantities();

    std::cout<<"Sim already recentered \n";

    start = chrono::steady_clock::now();
    for(size_t current_t=0;current_t<1000;current_t++ ){
        // for(size_t non_used_var=0;non_used_var<100;)
        // MemF.integrate(TS,sigma,kappa,H0,P,V0);
        // std::cout<<current_t <<" \n";
        if(current_t%10 ==0 && current_t<200){
            arcsim::Mesh remesher_mesh2 = translate_to_arcsim(mesh,geometry);
            Cloth_1.mesh=remesher_mesh2;
            Cloth_1.remeshing=remeshing_params;
            arcsim::compute_masses(Cloth_1);
            arcsim::compute_ws_data(Cloth_1.mesh);
            arcsim::dynamic_remesh(Cloth_1);
            delete mesh;
            delete geometry;
            std::tie(mesh_uptr, geometry_uptr) = translate_to_geometry(Cloth_1.mesh);
            arcsim::delete_mesh(Cloth_1.mesh);
            mesh = mesh_uptr.release();
            geometry = geometry_uptr.release();
            IM3DG = IMem3DG(mesh,geometry);
                    



        }
        
        
        
        
        if(false){
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








        
           if(current_t%1==0){
            Save_mesh(basic_name,current_t);
            Save_output_data  = true;
            Sim_data = std::ofstream(filename, std::ios_base::app);


        }
        
        
        if(current_t%50==0) {
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
        TS = 1.0;
        // if(current_t >20) TS = 1e-3;/
        // if(current_t >25 && current_t <50) TS = 1e-2;
        // if(current_t>50) TS = 1e-2;
        // if(current_t>100) TS = 1e-1;
        // if(current_t>200) TS = 2e-1;
        // if(current_t>400) TS = e-1;
        
        TS = IM3DG.integrate_Sob(TS,KB,Sim_data,Save_output_data);
        if(TS <0 ) {
            std::cout<<"The simulation Ended \n";
            std::cout<<"At the step " << current_t <<"\n";
            break;
        }
        time+=TS;

        Sim_data.close();
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



