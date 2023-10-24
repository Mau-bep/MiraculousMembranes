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


#include <chrono>


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"


#include "Mem-3dg.h"
#include "mem-3dg_implicit.h"
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
polyscope::SurfaceMesh* psMesh;

// Some global variables
float TIMESTEP = -4;
float kappa = 1.0;


float H0 = 1.0;
float V_bar= (4/3)*PI*10; 



float nu;
float Interaction_str;
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

Mem3DG M3DG;
IMem3DG IM3DG;
Bead Bead_1;



std::array<double, 3> BLUE = {0.11, 0.388, 0.89};
glm::vec<3, float> ORANGE_VEC = {1, 0.65, 0};
std::array<double, 3> ORANGE = {1, 0.65, 0};

// RemeshBoundaryCondition defaultRemeshOptions;
    



void flipZ() {
    // Rotate mesh 180 deg about up-axis on startup
    glm::mat4x4 rot = glm::rotate(glm::mat4x4(1.0f), static_cast<float>(PI), glm::vec3(0, 1, 0));
    for (Vertex v : mesh->vertices()) {
        Vector3 vec = geometry->inputVertexPositions[v];
        glm::vec4 rvec = {vec[0], vec[1], vec[2], 1.0};
        rvec = rot * rvec;
        geometry->inputVertexPositions[v] = {rvec[0], rvec[1], rvec[2]};
    }
    psMesh->updateVertexPositions(geometry->inputVertexPositions);
}

void showSelected() {
    // pass
}

void redraw() {
    psMesh->updateVertexPositions(geometry->inputVertexPositions);
    polyscope::requestRedraw();
}




void Save_mesh() {
   // Build member variables: mesh, geometry
    Vector3 Pos;
    std::ofstream o("Mesh.obj");
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



/*
 * User-defined buttons
 */
void functionCallback() {

    // if (ImGui::Button("Mean curvature flow")) {
        
    //     std::cout << "The mean edge length is = " << std::fixed << std::setprecision(8) << geometry->meanEdgeLength() << std::endl;

    //     std::cout<< "The mean edge length is = "<< geometry->meanEdgeLength()<<"\n";
    //     std::cout<< "The number of vertices is = "<< mesh->nVertices()<<"\n";
        

    //     MCF = MeanCurvatureFlow(mesh, geometry);
        
    //     MCF.integrate(TS);
    //     geometry->normalize(CoM, false);
    //     redraw();
    // }
    // if (ImGui::Button("Gauss curvature flow")) {
    //     GCF.integrate(TS);
    //     geometry->normalize(CoM, false);
    //     redraw();
    // }
   
    // if (ImGui::Button("Willmore flow" )) {
    //     WF2.integrate(TS);
    //     geometry->normalize(CoM, false);
    //     redraw();
    // }

    
    
    if (ImGui::Button("Tangential Vertex Smoothing " )) {

        avg_remeshing = smoothByCircumcenter(*mesh,*geometry);
        std::cout<< "The average amount the vertex were moved is "<< avg_remeshing << "\n";
        geometry->normalize(CoM, false);
        psMesh->remove();
        psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath("dodecahedra"), geometry->inputVertexPositions,
                                            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

        psMesh->setSurfaceColor({0.9607, 0.6627, 0.7215});
        // psMesh->setSurfaceColor({1.0, 0.45, 0.0}); // orange        
        polyscope::requestRedraw();
        redraw();
    }
    if (ImGui::Button("Remesh " )) {

        std::cout<< "The mean edge length is = "<< geometry->meanEdgeLength()<<"\n";
        std::cout<< "The number of vertices is = "<< mesh->nVertices()<<"\n";
                
        RemeshOptions Options;
        Options.targetEdgeLength=0.21;
        Options.curvatureAdaptation=0.1;
        Options.maxIterations=10;
        Options.minRelativeLength=0.2;
        Options.smoothStyle=RemeshSmoothStyle::Circumcentric;
        Options.boundaryCondition=RemeshBoundaryCondition::Tangential;
        ManifoldSurfaceMesh& meshu= *mesh;
        VertexPositionGeometry& geometryu = *geometry;
        MutationManager Mutation_manager(meshu,*geometry);
        
        remesh(*mesh,*geometry,Mutation_manager,Options);

        geometry->normalize(CoM, false);
        std::cout<< "The mean edge length after the remeshing is = "<< geometry->meanEdgeLength()<<"\n";
        std::cout<< "The number of vertices afther the remeshing is  = "<< mesh->nVertices()<<"\n";
        psMesh->remove();
        psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath("dodecahedra"), geometry->inputVertexPositions,
                                            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

        psMesh->setSurfaceColor({0.9607, 0.6627, 0.7215}); // not orange        
        polyscope::requestRedraw();
        redraw();
        
        
        

    }

    // if (ImGui::Button("Normal Flow")) {
    //     NF.integrate(TS);
    //     geometry->normalize(CoM, false);
    //     redraw();
    // }
    // if (ImGui::Button("Modified mean curvature flow")) {
    //     ModMCF.integrate(TS);
    //     geometry->normalize(CoM, false);
    //     redraw();
    // }


    // if (ImGui::Button("Membrane flow")) {
    //     // double TS=TIMESTEP/1000;
    //     for(int i=0;i<100;i++){
    //         MemF.integrate(TS,sigma,kappa,H0,P,V0);
        
    //     }
        
    //     geometry->normalize(CoM, false);
    //     redraw();
    // }

    if (ImGui::Button("Reset")) {
        geometry->inputVertexPositions = ORIG_VPOS;
        psMesh->updateVertexPositions(ORIG_VPOS);
        polyscope::requestRedraw();
    }
    


    ImGui::SliderFloat("Timestep", &TIMESTEP, -10.0, 10.0);
    TS=pow(10,TIMESTEP);
    

}


int main(int argc, char** argv) {

    
    nu=std::stod(argv[1]);
    c0=std::stod(argv[2]);
    KA=std::stod(argv[3]);
    KB=std::stod(argv[4]);
    Interaction_str=std::stod(argv[5]);
    // I will do it so i can give this values
    std::cout<<Interaction_str<<"\n";
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    

    
    TS=pow(10,-3);


    std::cout<< "Current path is " << argv[0];

    // std::string filepath = "../../../input/deformed_sphere_2.obj";
    std::string filepath = "../../../input/sphere.obj";
    // Load mesh
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    


    trgt_len=geometry->meanEdgeLength();
    V_bar=geometry->totalVolume();
    polyscope::options::autocenterStructures = true;

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = functionCallback;

    // Add mesh to GUI
    psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(filepath), geometry->inputVertexPositions,
                                            mesh->getFaceVertexList(), polyscopePermutations(*mesh));
    psMesh->setSurfaceColor({0.9607, 0.6627, 0.7215});    
    
    // Initialize operators.
    flipZ();
    

    ORIG_VPOS = geometry->inputVertexPositions;
    CoM = geometry->centerOfMass();
    
    // MCF = MeanCurvatureFlow(mesh, geometry);
    // ModMCF = ModifiedMeanCurvatureFlow(mesh, geometry);
    // NF =NormalFlow(mesh, geometry);
    // GCF = GaussCurvatureFlow(mesh, geometry);
    // WF = WillmoreFlow(mesh,geometry);
    // WF2 = WillmoreFlow2(mesh,geometry);
    // WFS = WillmoreFlowScho(mesh,geometry);
    Bead_1 = Bead(mesh,geometry,Vector3({3.4,3.7,-5.6}),1.0,Interaction_str);
    M3DG = Mem3DG(mesh,geometry,Bead_1);
    IM3DG = IMem3DG(mesh,geometry);
    

    // Add visualization options.
    psMesh->setSmoothShade(false);
    
    psMesh->setSurfaceColor({0.9607, 0.6627, 0.7215});// not orange
    // polyscope::screenshot("./This_filename_is_nice_right.jpg",true);

    // size_t counter=0;
    
    std::stringstream nustream;
    std::stringstream c0stream;
    std::stringstream KAstream;
    std::stringstream KBstream;
    std::stringstream Interactionstrstream;
    
    
    // std::stringstream H0stream;
    // std::stringstream kappastream;
    // std::stringstream sigmastream;
    std::stringstream Curv_adapstream;
    std::stringstream Min_rel_lengthstream;


    nustream << std::fixed << std::setprecision(3) << nu;
    c0stream << std::fixed << std::setprecision(3) << c0;
    KAstream << std::fixed << std::setprecision(3) << KA;
    KBstream << std::fixed << std::setprecision(6) << KB;
    Interactionstrstream << std::fixed << std::setprecision(6) << Interaction_str;



    Curv_adapstream << std::fixed << std::setprecision(2) << Curv_adap;
    Min_rel_lengthstream << std::fixed << std::setprecision(2) <<Min_rel_length;
    
    std::string first_dir="./Mem3DG_Beads/";
    int status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"If this name is 0 the directory was created succesfully "<< status ;

    first_dir="./Mem3DG_Beads/Curv_adap_"+Curv_adapstream.str()+"Min_rel_length_"+Min_rel_lengthstream.str();
    status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"\nIf this name is 0 the directory was created succesfully "<< status ;
    
    std::string basic_name ="./Mem3DG_Beads/Curv_adap_"+Curv_adapstream.str()+"Min_rel_length_"+Min_rel_lengthstream.str()+ "/nu_"+nustream.str()+"_c0_"+c0stream.str()+"_KA_"+KAstream.str()+"_KB_"+KBstream.str()+"_Inter_"+ Interactionstrstream.str()+"/";
    status = mkdir(basic_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    std::cout<<"\nIf this number is 0 the directory was created succesfully "<< status<<"\n" ;

    std::string filename = basic_name+"Output_data.txt";
    std::string filename2 =basic_name + "Bead_data.txt";

    std::ofstream Sim_data(filename);
    std::ofstream Bead_data(filename2);

    bool Save_bead_data=false;
    Bead_data<<"####### This data is taken every 1000 steps just like the mesh\n";
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
    double nu_step;

    start = chrono::steady_clock::now();
    for(size_t current_t=0;current_t<100000;current_t++ ){ //I took two zeros
        // for(size_t non_used_var=0;non_used_var<100;)
        // MemF.integrate(TS,sigma,kappa,H0,P,V0);
        if(true){
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
        }

        psMesh->remove();
        
        psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath("dodecahedra"), geometry->inputVertexPositions,
                                            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

        psMesh->setSurfaceColor({0.9607, 0.6627, 0.7215});
        psMesh->setEdgeWidth(1.0);

        

        
        if(current_t%500==0) {
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
                // nu_= (nu-nu_0)/50;
            }

            std::cout<< "The spontaneous curvature is " << H0<< "\n";
            // std::cout << "The system time is " << M3DG.system_time <<"\n\n";
            
            std::cout<<"A thousand iterations took "<<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" miliseconds\n";


            polyscope::screenshot(basic_name+std::to_string(current_t)+".jpg",true);
            start = chrono::steady_clock::now();

        }
        nu_evol= time<50 ? nu_0 + (nu-nu_0)*time/50 : nu; 
        // nu_evol=nu;

        // dt_sim=M3DG.integrate(TS,V_bar,nu_evol,c0,P0,KA,KB,Kd,Sim_data,time);
        // dt_sim=M3DG.integrate(TS,V_bar,nu_evol,c0,P0,KA,KB,Kd,Sim_data,time,true);
        dt_sim=M3DG.integrate(TS,V_bar,nu_evol,c0,P0,KA,KB,Kd,Sim_data,time,Save_bead_data,Bead_data);
       
        // IM3DG.integrate(TS,nu_evol,V_bar,P0,KA,KB);
        // dt_sim=TS;
        Vector3 center =geometry->centerOfMass();
        
        geometry->normalize(Vector3({0.0,0.0,0.0}),false);
        geometry->refreshQuantities();
        
        if(dt_sim==-1){
            std::cout<<"Sim broke\n";
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



