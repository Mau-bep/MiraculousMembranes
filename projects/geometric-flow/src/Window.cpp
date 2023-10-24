// GEOMETRIC FLOW

// #include <stdlib.h>
#include <unistd.h>


#include <sys/stat.h>
#include <iostream>
#include <string>



#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/surface/remeshing.h"





#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "mean-curvature-flow.h"
// #include "modified-mean-curvature-flow.h"
#include "normal-flow.h"
// #include "gauss-curvature-flow.h"
// #include "willmore-flow.h"
// #include "willmore-flow-2.h"
// #include "willmore-flow-scho.h"
#include "Mem-3dg.h"
// #include "membrane-flow.h"
#include "mem-3dg_implicit.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


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
float V_bar; 
double init_vol;


float nu;
float c0;

float P0=10000.0;
float KA=10.0;
float KB=0.01;
// float Kd=1.0;

double TS=0.0001;

double Curv_adap=0.1;
double Min_rel_length=0.5;
double trgt_len;
double avg_remeshing;
bool edge_length_adj;
VertexData<Vector3> ORIG_VPOS; // original vertex positions
Vector3 CoM;                   // original center of mass
MeanCurvatureFlow MCF;
// ModifiedMeanCurvatureFlow ModMCF;
NormalFlow NF;
// GaussCurvatureFlow GCF;



// WillmoreFlow2 WF2;
// WillmoreFlowScho WFS;
Mem3DG M3DG;
IMem3DG IM3DG;


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


    if (ImGui::Button("Mean curvature flow")) {
        
        
        std::cout<< "The mean edge length is = "<< geometry->meanEdgeLength()<<"\n";
        std::cout<< "The number of vertices is = "<< mesh->nVertices()<<"\n";
        
        MCF.integrate(TS);
        geometry->normalize(CoM, false);
        geometry->refreshQuantities();
        redraw();
        // Save_mesh();
    }
    if (ImGui::Button("Implicit M3DG")) {
        
        
        // std::cout<< "The mean edge length is = "<< geometry->meanEdgeLength()<<"\n";
        // std::cout<< "The number of vertices is = "<< mesh->nVertices()<<"\n";
        
        IM3DG.integrate(TS,nu,V_bar,P0,KA,KB);
        geometry->normalize(CoM, false);
        geometry->refreshQuantities();
        redraw();
        std::cout<<"The new volume is"<<geometry->totalVolume()<<"\n";
        // Save_mesh();
    }
    if (ImGui::Button("Remesh " )) {

        std::cout<< "The mean edge length is = "<< geometry->meanEdgeLength()<<"\n";
        std::cout<< "The number of vertices is = "<< mesh->nVertices()<<"\n";
                
        RemeshOptions Options;
        Options.targetEdgeLength=trgt_len;
        Options.curvatureAdaptation=0.1;
        Options.maxIterations=10;
        Options.minRelativeLength=0.5;
        Options.smoothStyle=RemeshSmoothStyle::Circumcentric;
        Options.boundaryCondition=RemeshBoundaryCondition::Tangential;

        MutationManager Mutation_manager(*mesh,*geometry);
        
        remesh(*mesh,*geometry,Mutation_manager,Options);

        geometry->normalize(CoM, false);
        std::cout<< "The mean edge length after the remeshing is = "<< geometry->meanEdgeLength()<<"\n";
        std::cout<< "The number of vertices afther the remeshing is  = "<< mesh->nVertices()<<"\n";
        psMesh->remove();
        psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath("dodecahedra"), geometry->inputVertexPositions,
                                            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

        psMesh->setSurfaceColor({1.0, 0.45, 0.0}); // orange        
        polyscope::requestRedraw();
        redraw();
        
        
        

    }

    
    
    if (ImGui::Button("Normal Flow")) {
        NF.integrate(TS);
        geometry->normalize(CoM, false);
        geometry->refreshQuantities();
        redraw();
        std::cout<<"The new volume is"<<geometry->totalVolume()<<"\n";
    }
    if (ImGui::Button("Mem 3DG")) {
        // for(int m =0; m<100 ; m++){
        // M3DG.integrate(TS,V_bar,nu,c0,P0,KA,KB,Kd);
        
        geometry->normalize(CoM, false);
        redraw();
    // }
    }

    // if (ImGui::Button("Membrane flow")) {
    //     // double TS=TIMESTEP/1000;
    //     for(int i=0;i<1;i++){
    //         MemF.integrate(TS,sigma,kappa,H0,P,V_bar);
        
    //     }
        
    //     geometry->normalize(CoM, false);
    //     redraw();
    // }

    if (ImGui::Button("Reset")) {
        geometry->inputVertexPositions = ORIG_VPOS;
        psMesh->updateVertexPositions(ORIG_VPOS);
        polyscope::requestRedraw();
    }
    
    if (ImGui::Button("Savemesh")) {
        Save_mesh();
    
    }
    


    ImGui::SliderFloat("Timestep", &TIMESTEP, -10.0, 10.0);
    TS=pow(10,TIMESTEP);
    
        
    ImGui::SliderFloat("nu", &nu, 0.0, 1.0);
    ImGui::SliderFloat("c0", &c0, -4.0, 4.0);
    
    ImGui::SliderFloat("KA", &KA, 0.0, 1.0);
    ImGui::SliderFloat("KB", &KB, 0.0, 1.0);
    
    // ImGui::SliderFloat("P0", &P0, 0.0, 1.0);
    
    ImGui::SliderFloat("V_bar", &V_bar, 10.0, 100.0);
}



int main(int argc, char** argv) {

    char buffer[PATH_MAX];
   if (getcwd(buffer, sizeof(buffer)) != NULL) {
       printf("Current working directory : %s\n", buffer);
   } else {
       perror("getcwd() error");
   }
   
    // I will do it so i can give this values
    // H0=std::stod(argv[1]);
    // kappa=std::stod(argv[2]);
    // sigma=std::stod(argv[3]);
    // Min_rel_length = std::stod(argv[4]);
    // Curv_adap= std::stod(argv[5]);
    nu=std::stod(argv[1]);
    // c0=std::stod(argv[2]);
    c0=0.0;
    



    TS=pow(10,-6);
    // P=3.0;
    // 

    // Configure the argument parser
   
    // args::ArgumentParser parser("15-458 HW3");
    // args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");
    
    // // Parse args
    // try {
    //     parser.ParseCLI(argc, argv);
    // } catch (args::Help&) {
    //     std::cout << parser;
    //     return 0;
    // } catch (args::ParseError& e) {
    //     std::cerr << e.what() << std::endl;
    //     std::cerr << parser;
    //     return 1;
    // }

    // If a mesh name was not given, use sphere mesh.
    std::cout<< "Current path is " << argv[0];

    //this is for debugging   
    // std::string filepath = "./ddg-exercises/input/sphere.obj";
    
    // std::string filepath = "../../../input/bunny.obj";
    std::string filepath = "../../../input/Simple_cil_regular.obj";
    // std::string filepath = "../../../../Cluster_Folders/Results/projects/geometric-flow/build/Mem3DG_IMG_parallel/Curv_adap_0.10Min_rel_length_0.50/nu_0.800_c0_0.000_KA_11.000_KB_0.001000/150000.obj";

    // if (inputFilename) {
    //     filepath = args::get(inputFilename);
    // }
   
    // Load mesh
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    
    trgt_len=geometry->meanEdgeLength();
    V_bar=geometry->totalVolume();
    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = functionCallback;

    // Add mesh to GUI
    psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(filepath), geometry->inputVertexPositions,
                                            mesh->getFaceVertexList(), polyscopePermutations(*mesh));
    psMesh->setSurfaceColor({1.0, 0.45, 0.0}); // orange
    

    
    // Initialize operators.
    flipZ();
    

    ORIG_VPOS = geometry->inputVertexPositions;
    CoM = geometry->centerOfMass();
    MCF = MeanCurvatureFlow(mesh, geometry);
    // ModMCF = ModifiedMeanCurvatureFlow(mesh, geometry);
    NF =NormalFlow(mesh, geometry);
    // GCF = GaussCurvatureFlow(mesh, geometry);

    // WF2 = WillmoreFlow2(mesh,geometry);
    // WFS = WillmoreFlowScho(mesh,geometry);
    M3DG = Mem3DG(mesh,geometry);
    IM3DG = IMem3DG(mesh,geometry);

    // Add visualization options.
    psMesh->setSmoothShade(false);
    psMesh->setSurfaceColor({1.0, 0.45, 0.0}); // orange
 
    size_t counter=0;

    // std::stringstream H0stream;
    // std::stringstream kappastream;
    // std::stringstream sigmastream;
    // std::stringstream Curv_adapstream;
    // std::stringstream Min_rel_lengthstream;

    // H0stream << std::fixed << std::setprecision(2) << H0;
    // kappastream << std::fixed << std::setprecision(2) << kappa;
    // sigmastream << std::fixed << std::setprecision(2) << sigma;
    // Curv_adapstream << std::fixed << std::setprecision(2) << Curv_adap;
    // Min_rel_lengthstream << std::fixed << std::setprecision(2) <<Min_rel_length;
    
   






    std::cout<<"The initial volume is "<< geometry->totalVolume()<<"\n";

    // Give control to the polyscope gui
    polyscope::show();



    delete mesh;
    delete geometry;

    return EXIT_SUCCESS;
}


