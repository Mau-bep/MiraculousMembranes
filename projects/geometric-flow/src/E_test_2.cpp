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
double TS=0.0001;

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

// RemeshBoundaryCondition defaultRemeshOptions;
// void redraw() {
//     psMesh->updateVertexPositions(geometry->inputVertexPositions);
//     polyscope::requestRedraw();
// }

// void functionCallback() {

    
//     if (ImGui::Button("Tangential Vertex Smoothing " )) {

//         avg_remeshing = smoothByCircumcenter(*mesh,*geometry);
//         std::cout<< "The average amount the vertex were moved is "<< avg_remeshing << "\n";
//         geometry->normalize(CoM, false);
//         psMesh->remove();
//         psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath("dodecahedra"), geometry->inputVertexPositions,
//                                             mesh->getFaceVertexList(), polyscopePermutations(*mesh));

//         psMesh->setSurfaceColor({0.9607, 0.6627, 0.7215});
//         // psMesh->setSurfaceColor({1.0, 0.45, 0.0}); // orange        
//         polyscope::requestRedraw();
//         redraw();
//     }
//     if (ImGui::Button("Remesh " )) {

//         std::cout<< "The mean edge length is = "<< geometry->meanEdgeLength()<<"\n";
//         std::cout<< "The number of vertices is = "<< mesh->nVertices()<<"\n";
                
//         RemeshOptions Options;
//         Options.targetEdgeLength=0.21;
//         Options.curvatureAdaptation=0.1;
//         Options.maxIterations=10;
//         Options.minRelativeLength=0.2;
//         Options.smoothStyle=RemeshSmoothStyle::Circumcentric;
//         Options.boundaryCondition=RemeshBoundaryCondition::Tangential;
//         ManifoldSurfaceMesh& meshu= *mesh;
//         VertexPositionGeometry& geometryu = *geometry;
//         MutationManager Mutation_manager(meshu,*geometry);
        
//         remesh(*mesh,*geometry,Mutation_manager,Options);

//         geometry->normalize(CoM, false);
//         std::cout<< "The mean edge length after the remeshing is = "<< geometry->meanEdgeLength()<<"\n";
//         std::cout<< "The number of vertices afther the remeshing is  = "<< mesh->nVertices()<<"\n";
//         psMesh->remove();
//         psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath("dodecahedra"), geometry->inputVertexPositions,
//                                             mesh->getFaceVertexList(), polyscopePermutations(*mesh));

//         psMesh->setSurfaceColor({0.9607, 0.6627, 0.7215}); // not orange        
//         polyscope::requestRedraw();
//         redraw();
        
        
        

//     }



//     if (ImGui::Button("Reset")) {
//         geometry->inputVertexPositions = ORIG_VPOS;
//         psMesh->updateVertexPositions(ORIG_VPOS);
//         polyscope::requestRedraw();
//     }
    


//     ImGui::SliderFloat("Timestep", &TIMESTEP, -10.0, 10.0);
//     TS=pow(10,TIMESTEP);
    

// }

// void flipZ() {
//     // Rotate mesh 180 deg about up-axis on startup
//     glm::mat4x4 rot = glm::rotate(glm::mat4x4(1.0f), static_cast<float>(PI), glm::vec3(0, 1, 0));
//     for (Vertex v : mesh->vertices()) {
//         Vector3 vec = geometry->inputVertexPositions[v];
//         glm::vec4 rvec = {vec[0], vec[1], vec[2], 1.0};
//         rvec = rot * rvec;
//         geometry->inputVertexPositions[v] = {rvec[0], rvec[1], rvec[2]};
//     }
//     psMesh->updateVertexPositions(geometry->inputVertexPositions);
// }




void Save_mesh(std::string basic_name, size_t current_t) {
   // Build member variables: mesh, geometry
    Vector3 Pos;
    std::ofstream o(basic_name+std::to_string(current_t)+".obj");
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




int main(int argc, char** argv) {

    
    nu=std::stod(argv[1]);
    c0=std::stod(argv[2]);
    KA=std::stod(argv[3]);
    KB=std::stod(argv[4]);
    
    size_t target_index=120;
    // I will do it so i can give this values
 
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    
    bool with_bead=false;

    bool Save_data=false;
    TS=pow(10,-3);


    std::cout<< "Current path is " << argv[0]<<"\n";

    std::string filepath = "../../../input/sphere.obj";
    // std::string filepath = "../../../input/bloodcell_4k.obj";
    // std::string filepath = "../../../input/Simple_cil_regular.obj";
    
    // std::string filepath = "../../../input/bloodcell.obj";
    // std::string filepath = "../input/sphere.obj"; //this is for debug
    
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    std::cout<<"The mesh is correctly loaded\n";
    



    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    
    trgt_len=geometry->meanEdgeLength();
    V_bar=geometry->totalVolume();
    
    EdgeData<int> No_remesh(*mesh,0);
    No_remesh[0]=1;




    ORIG_VPOS = geometry->inputVertexPositions;
    CoM = geometry->centerOfMass();
    double radius=1.0;
    double Interaction_str=1.0;
    Bead_1 = Bead(mesh,geometry,Vector3({5.8,0.0,0.0}),radius,Interaction_str);
    
    M3DG = Mem3DG(mesh,geometry,Bead_1);



    // Add visualization options.
    // psMesh->setSmoothShade(false);
    
    // psMesh->setSurfaceColor({0.9607, 0.6627, 0.7215});// not orange
    
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
    std::string basic_name;
    std::string first_dir;

    if(with_bead){

    first_dir="../Results/Tests_bead/";
    int status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"If this name is 0 the directory was created succesfully "<< status ;

    // first_dir="./Tests/";
    // status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"\nIf this name is 0 the directory was created succesfully "<< status ;
    
    basic_name ="../Results/Tests_bead/nu_"+nustream.str()+"_c0_"+c0stream.str()+"_KA_"+KAstream.str()+"_KB_"+KBstream.str()+"/";
    status = mkdir(basic_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"\nIf this number is 0 the directory was created succesfully "<< status<<"\n" ;
    }
    else{

    first_dir="../Results/Tests_general/";
    int status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"If this name is 0 the directory was created succesfully "<< status ;

    // first_dir="./Tests/";
    // status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"\nIf this name is 0 the directory was created succesfully "<< status ;
    
    basic_name ="../Results/Tests_general/nu_"+nustream.str()+"_c0_"+c0stream.str()+"_KA_"+KAstream.str()+"_KB_"+KBstream.str()+"/";
    status = mkdir(basic_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        
    }

    std::string filename_basic = basic_name+"Output_data.txt";

    std::ofstream Sim_data;
    

    std::string filename2 = basic_name + "Bead_data.txt";
    std::ofstream Bead_data(filename2);
    bool Save_bead_data=false;
        
    if(with_bead){
        
        Bead_data<<"####### This data is taken every 250 steps just like the mesh radius is " << radius<<" \n";
    
    }

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
    double time=0.01;
    double dt_sim=0.0;


    // Here i have loaded the mesh and i am ready to start testing


    // I may need to run this and save things so i can se the evolution...
    std::ofstream Gradient_data;
    std::ofstream Gradient_data_area;
    std::ofstream Gradient_data_bending;
    std::ofstream Gradient_data_bending_norms;
    std::ofstream Gradient_data_tot_area;

    std::string filename;

    if(with_bead){
    std::ofstream Gradient_data_bead;

    std::ofstream Gradient_data_bead_dx;
    


    filename = basic_name+"Bead_Gradient_evaluation_dx_"+std::to_string(0) + ".txt";
    Gradient_data_bead_dx.open(filename);
    Gradient_data_bead_dx << "Bead_gradient dx evaluation\n";
    M3DG.Grad_Bead_dx(Gradient_data_bead_dx,true);

    Gradient_data_bead_dx.close();

    }


    // Since this is the test secion, i want to check something
    double Total_A=0;
    double Total_A_dual_bar=0;
    double Total_A_dual_circ=0;
    Total_A=geometry->totalArea();
    double max_x=0;
    double max_y=0;
    // double 
    for( Vertex v : mesh->vertices()){
        if(max_x<geometry->inputVertexPositions[v].x){
            max_x=geometry->inputVertexPositions[v].x;
        }    
        Total_A_dual_bar+=geometry->barycentricDualArea(v);
        Total_A_dual_circ+=geometry->circumcentricDualArea(v);
    }

    std::cout<<"THe areas a are: Total "<<Total_A <<" Barycentric "<< Total_A_dual_bar <<" and Circumcentric "<< Total_A_dual_circ <<" \n";
    
    Volume=geometry->totalVolume();
    std::cout<<"The radius (according to area) is " << sqrt(Total_A/(4*PI))<<"\n";
    std::cout<<"THe radius (according to volume) is" << pow(  3*Volume/(4*PI)   ,1.0/3.0)<<"\n"; 
    std::cout<<"THe maximum x coordinate is "<< max_x <<"\n";

    for(size_t current_t=0; current_t <1000;current_t++){
    if(current_t==0){
    Sim_data.open(filename_basic);
    }
    else{
    Sim_data.open(filename_basic,std::ios_base::app);
    }
   
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
        // Options.remesh_list=true;
        // Options.No_remesh_list=No_remesh;
        MutationManager Mutation_manager(*mesh,*geometry);
        remesh(*mesh,*geometry,Mutation_manager,Options);
        n_vert_new=mesh->nVertices();
        counter=counter+1; 
        }
        }

        n_vert=mesh->nVertices();
        std::cout<< "THe number of vertices is "<< n_vert <<"\n";    
        std::cout << "The avg edge length is = " << std::fixed << std::setprecision(10) << geometry->meanEdgeLength() << std::endl;

    if(current_t%100==0){

            Save_mesh(basic_name,current_t);
            if(with_bead)
            {
                Save_bead_data=false;
            }
        }
    
    if(current_t%200==0){
        filename = basic_name+"Vol_Gradient_evaluation_"+std::to_string(current_t) + ".txt";
        // std::ofstream Gradient_data(filename);
        Gradient_data.open(filename);
        // std::ofstream o(basic_name+std::to_string(current_t)+".obj");
        Gradient_data<< "Volume grad\n";
        M3DG.Grad_Vol(Gradient_data,P0,V_bar,true);

        Gradient_data.close();
        double A_bar=4*PI*pow(3*V_bar/(4*PI*nu_evol),2.0/3.0);

        filename = basic_name+"Area_Gradient_evaluation_"+std::to_string(current_t) + ".txt";
        Gradient_data_area.open(filename);
        Gradient_data_area<< "Area grad\n";
        M3DG.Grad_Area(Gradient_data_area,A_bar,KA,true);

        Gradient_data_area.close();

        filename = basic_name+"Bending_Gradient_evaluation_"+std::to_string(current_t) + ".txt";
        Gradient_data_bending.open(filename);
        Gradient_data_bending<< "Bending grad\n";
        double H_bar=sqrt(4*PI/A_bar)*c0/2.0;
        M3DG.Grad_Bending(Gradient_data_bending,H_bar,KB,true);
        Gradient_data_bending.close();


        filename = basic_name+"Bending_Gradient_evaluation_2_"+std::to_string(current_t) + ".txt";
        Gradient_data_bending_norms.open(filename);
        Gradient_data_bending_norms<< "Bending grad norm diff\n";
        // double H_bar=sqrt(4*PI/A_bar)*c0/2.0;
        M3DG.Grad_Bending_2(Gradient_data_bending_norms,H_bar,KB);
        Gradient_data_bending_norms.close();
        


        filename = basic_name+"Area_tot_Gradient_evaluation_"+std::to_string(current_t) + ".txt";
        Gradient_data_tot_area.open(filename);
        Gradient_data_tot_area<<"Area tot grad\n";
        M3DG.Grad_tot_Area(Gradient_data_tot_area,true);
        Gradient_data_tot_area.close();


        // filename = basic_name+"Bead_Gradient_evaluation_"+std::to_string(current_t) + ".txt";
        // Gradient_data_bead.open(filename);
        // Gradient_data_bead<< "Bead grad\n";
        // // double H_bar=sqrt(4*PI/A_bar)*c0/2.0;
        // M3DG.Grad_Bead(Gradient_data_bead,true,false);
        // Gradient_data_bead.close();


        Volume= geometry->totalVolume();
        Area=geometry->totalArea();
        nu_obs=3*Volume/(4*PI*pow( Area/(4*PI) ,1.5 ));
        H0=sqrt(4*PI/Area)*c0/2.0;
        
        if(current_t==0){
            nu_0=nu_obs;
        }
   

    }
     
    nu_evol= time<50 ? nu_0 + (nu-nu_0)*time/50 : nu; 
    if(with_bead){
        dt_sim=M3DG.integrate(TS,V_bar,nu_evol,c0,P0,KA,KB,Kd,Sim_data,time,Save_bead_data,Bead_data,Save_data);
        
    }
    else{
    dt_sim=M3DG.integrate(TS,V_bar,nu_evol,c0,P0,KA,KB,Kd,Sim_data,time,Save_data);
    }

    if(dt_sim==-1)
        {
        std::cout<<"Sim broke or timestep very small\n";
        break;
        }
    else{
        time+=dt_sim;
            
        }

    

    Sim_data.close();

    // filename = basic_name+"Bending_Gradient_analitical_terms.txt";
    // std::ofstream Gradient_data_bending_2(filename);
    // Gradient_data_bending_2<< "Bending gradients per component\n";
    // // double H_bar=sqrt(4*PI/A_bar)*c0/2.0;
    // M3DG.Bending_test(Gradient_data_bending_2,H_bar,KB);

    // Gradient_data_bending_2.close();

    
    
    
    }





    delete mesh;
    delete geometry;

    return EXIT_SUCCESS;
}



