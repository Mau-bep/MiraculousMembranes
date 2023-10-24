// GEOMETRIC FLOW

// #include <stdlib.h>
#include <unistd.h>


#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <string>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/surface/remeshing.h"





#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

// #include "mean-curvature-flow.h"
// #include "modified-mean-curvature-flow.h"
// #include "normal-flow.h"
// #include "gauss-curvature-flow.h"
// #include "willmore-flow.h"
// #include "willmore-flow-2.h"
// #include "willmore-flow-scho.h"
#include "Mem-3dg.h"

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
float V_bar= (4/3)*PI*10; 



float nu;

float c0;


float P0=10000.0;
float KA=1.0000;
float KB=0.0001;
float Kd=0.0;
double TS=0.0001;

double Curv_adap=0.1;
double Min_rel_length=0.5;
double trgt_len;
double avg_remeshing;
bool edge_length_adj;
VertexData<Vector3> ORIG_VPOS; // original vertex positions
Vector3 CoM;                   // original center of mass
// MeanCurvatureFlow MCF;
// ModifiedMeanCurvatureFlow ModMCF;
// NormalFlow NF;
// GaussCurvatureFlow GCF;


// WillmoreFlow WF;
// WillmoreFlow2 WF2;
// WillmoreFlowScho WFS;
Mem3DG M3DG;


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
  
    // I will do it so i can give this values
 
    
    TS=pow(10,-3);


    std::cout<< "Current path is " << argv[0];

    // std::string filepath = "../../../input/sphere.obj";
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

    size_t N_vert=mesh->nVertices();
    size_t index;
    Vector<double> Scalar_MC(N_vert);
    double Area_test=geometry->totalArea();
    double total_mean_curv=0;
    double mean_mean_curv=0;
    double r= pow((3.0/(4.0*PI))*V_bar,1.0/3.0); 
    double r_2=pow( Area_test/(4*PI),0.5 );

    for(Vertex v : mesh->vertices()) {
        index=v.getIndex();
        Scalar_MC.coeffRef(index)=geometry->scalarMeanCurvature(v)/geometry->barycentricDualArea(v);
        total_mean_curv+=geometry->scalarMeanCurvature(v);
        mean_mean_curv+=Scalar_MC[index];
        // Vertex_area.coeffRef(index)=geometry->barycentricDualArea(v);
        }   
    mean_mean_curv/=N_vert;
    std::cout << "The radius is"<< r<<"\n";
    std::cout << "The radius is"<< r_2<<"\n";
    
    std::cout<< "The theoretical H is "<<1/r<<"\n";
    std::cout<< "The measured H is "<< mean_mean_curv <<"\n";
    std::cout<<" The theoretical integrated H  is"<< 4*PI*r <<"\n";
    std::cout<< "The integrated H is "<< total_mean_curv<<"\n";



    // std::cout<<Scalar_MC<<"\n";





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
    M3DG = Mem3DG(mesh,geometry);

    // Add visualization options.
    psMesh->setSmoothShade(false);
    
    psMesh->setSurfaceColor({0.9607, 0.6627, 0.7215});// not orange
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
    
    std::string first_dir="./test_IMG/";
    int status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"If this name is 0 the directory was created succesfully "<< status ;

    first_dir="./test_IMG/Curv_adap_"+Curv_adapstream.str()+"Min_rel_length_"+Min_rel_lengthstream.str();
    status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"\nIf this name is 0 the directory was created succesfully "<< status ;
    
    std::string basic_name ="./test_IMG/Curv_adap_"+Curv_adapstream.str()+"Min_rel_length_"+Min_rel_lengthstream.str()+ "/nu_"+nustream.str()+"_c0_"+c0stream.str()+"_KA_"+KAstream.str()+"_KB_"+KBstream.str()+"/";
    status = mkdir(basic_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    std::cout<<"\nIf this number is 0 the directory was created succesfully "<< status<<"\n" ;

    std::string filename = basic_name+"Output_data.txt";

    std::ofstream Sim_data(filename);
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

    for(size_t current_t=0;current_t<0;current_t++ ){
        // for(size_t non_used_var=0;non_used_var<100;)
        // MemF.integrate(TS,sigma,kappa,H0,P,V0);
        if(true){
        // if(current_t%10==0 ){
        n_vert_old=mesh->nVertices();
        n_vert_new=1;

        counter=0;
        while(n_vert_new!=n_vert_old && counter<10){
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

        

        
        if(current_t%50==0) {
            n_vert=mesh->nVertices();
            std::cout<< "THe number of vertices is "<< n_vert <<"\n";    
            std::cout<< "the avg edge lenghth is " <<geometry->meanEdgeLength()<<"\n";
            
            Volume= geometry->totalVolume();
            Area=geometry->totalArea();
            nu_obs=3*Volume/(4*PI*pow( Area/(4*PI) ,1.5 ));
            H0=sqrt(4*PI/Area)*c0/2.0;



            // c_null=2*H0 *pow(Area/(4*PI),0.5);
            std::cout<< "The reduced volume is "<< nu_obs << "\n";
            if(current_t==0){
                nu_0=nu_obs;
            }

            std::cout<< "The spontaneous curvature is " << H0<< "\n\n";

            polyscope::screenshot(basic_name+std::to_string(current_t)+".jpg",true);
        }
        nu_evol= time<50 ? nu_0 + (nu-nu_0)*time/50 : nu; 



        time+=M3DG.integrate(TS,V_bar,nu_evol,c0,P0,KA,KB,Kd,Sim_data,time);
        // MCF.integrate(TS);


        


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



