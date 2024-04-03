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


// size_t get_init_state_index()

int main(int argc, char** argv) {




    nu=std::stod(argv[1]);
    c0=std::stod(argv[2]);
    KA=std::stod(argv[3]);
    KB=std::stod(argv[4]);
    


    bool DoTest = true;
    bool Test_Area = false;
    bool Test_Vol = false;
    bool Test_Bead = true;
    bool Test_bending = false;




    // std::string path = "/home/mrojasve/Documents/DDG/MiraculousMembranes/projects/geometric-flow/Results/Tests_cil_regular/Minus_sign_force";
  
    // DIR *dh;
    // struct dirent * contents; 
    // dh = opendir("/home/mrojasve/Documents/DDG/MiraculousMembranes/projects/geometric-flow/Results/Bunny/nu_0.000_c0_0.000_KA_10.000_KB_0.000000/");
    // int index_underscore;
    // int index_dot;
    // size_t max_index=0;
    // size_t current_index=0;
    
    // if ( !dh )
    // {
    //     std::cout << "The given directory is not found";
    //     return 1;
    // }
    // while ( ( contents = readdir ( dh ) ) != NULL )
    // {
    //     std::string name = contents->d_name;
    //     // std::cout<<name.size()<<endl;
        
    //     if(name.size()>7 && name.find('t')>10){
    //     index_underscore=name.find('_');
    //     index_dot=name.find('.');
    //     string str_index = name.substr(index_underscore+1,index_dot-index_underscore-1 );
    //     std::cout << str_index << endl;
    //     // I need to avoid the final state
    //     current_index=stoi(str_index);
    //     std::cout<<current_index<<endl;
    //     if(current_index>max_index){
    //         max_index=current_index;
    //     }

    //     // std::cout << name.substr(2) << endl;
    //     }
    // }
    // std::cout<< "The biggest index is "<< max_index << endl;
    // closedir ( dh );
    
    if(DoTest){


    size_t target_index=120;
    // I will do it so i can give this values
 
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    
    bool with_bead=true;

    
    TS=pow(10,-3);


    std::cout<< "Current path is " << argv[0]<<"\n";

    // std::string filepath = "../../../input/Simple_cil_regular.obj";
    // std::string filepath = "../../../input/sphere.obj";
    // std::string triangle_path ="../input/triangle_fan.obj"; 
    std::string triangle_path ="../../../input/triangle_fan.obj"; 
    
    std::ofstream Triangle_fan(triangle_path);

    std::uniform_real_distribution<double> unif(0.0,1.0);
    std::default_random_engine re;

    Triangle_fan<<"# This is a triangle fan developed only for testing \n";

    Triangle_fan<<"v "<<unif(re) << " " << unif(re)<<" " << unif(re)<<"\n";
    Triangle_fan<<"v "<<unif(re) << " " << unif(re)<<" " << unif(re)<<"\n";
    Triangle_fan<<"v "<<unif(re) << " " << unif(re)<<" " << unif(re)<<"\n";
    Triangle_fan<<"v "<<unif(re) << " " << unif(re)<<" " << unif(re)<<"\n";
    Triangle_fan<<"v "<<unif(re) << " " << unif(re)<<" " << unif(re)<<"\n";
    Triangle_fan<<"v "<<unif(re) << " " << unif(re)<<" " << unif(re)<<"\n";
    Triangle_fan<<"v "<<unif(re) << " " << unif(re)<<" " << unif(re)<<"\n";
    
    Triangle_fan<<"f 1 2 4 \n";
    Triangle_fan<<"f 2 5 4 \n";
    Triangle_fan<<"f 4 5 7 \n";
    Triangle_fan<<"f 4 7 6 \n";
    Triangle_fan<<"f 4 6 3 \n";
    Triangle_fan<<"f 4 3 1 \n";
    


    double a_random_double = unif(re);
    std::cout<<a_random_double<<"\n";

    Triangle_fan.close();
    
    std::string filepath = "../../../input/sphere.obj"; //this is for debug
    
    // std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);

    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(triangle_path);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    
    trgt_len=geometry->meanEdgeLength();
    V_bar=geometry->totalVolume();
  

    ORIG_VPOS = geometry->inputVertexPositions;
    CoM = geometry->centerOfMass();
    double radius=1.0;
    double Interaction_str=1.0;
    Bead_1 = Bead(mesh,geometry,Vector3({0.0,0.0,0.0}),radius,Interaction_str);
    Bead_1.interaction="pulling";
    // M3DG = Mem3DG(mesh,geometry);
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
    
    std::string first_dir="../Results/Tests_fan/";
    int status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"If this name is 0 the directory was created succesfully "<< status ;

    // first_dir="./Tests/";
    // status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"\nIf this name is 0 the directory was created succesfully "<< status ;
    
    std::string basic_name ="../Results/Tests_fan/nu_"+nustream.str()+"_c0_"+c0stream.str()+"_KA_"+KAstream.str()+"_KB_"+KBstream.str()+"/";
    status = mkdir(basic_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // std::cout<<"\nIf this number is 0 the directory was created succesfully "<< status<<"\n" ;

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
    std::ofstream Gradient_data_bead;
    std::ofstream Gradient_data_tot_area;

    std::ofstream Gradient_data_bead_dx;
    std::string filename;



    // OK first the sanity check
    std::cout<< "Are all areas positive? " << M3DG.Area_sanity_check()<< "\n ";



    // Since this is the test secion, i want to check something
    double Total_A=0;
    double Total_A_dual_bar=0;
    double Total_A_dual_circ=0;
    Total_A=geometry->totalArea();
    for( Vertex v : mesh->vertices()){
        Total_A_dual_bar+=geometry->barycentricDualArea(v);
        Total_A_dual_circ+=geometry->circumcentricDualArea(v);
    }

    std::cout<<"THe areas a are: Total "<<Total_A <<" Barycentric "<< Total_A_dual_bar <<" and Circumcentric "<< Total_A_dual_circ <<" \n";
    

    // Here the tests start

    // According to what i have written i need to create a function that receives the vertex where the angle is and the vertex where i want to calculate the gradient




    if(Test_Area){


    filename = basic_name+"Area_tot_Gradient_evaluation_"+std::to_string(0) + ".txt";
    Gradient_data_tot_area.open(filename);
    Gradient_data_tot_area<<"Area tot grad first finite then theory\n";
    VertexData<Vector3> Area_grad_finite = M3DG.Grad_tot_Area(Gradient_data_tot_area,false);
    VertexData<Vector3> Area_grad_th = M3DG.SurfaceGrad();
    Gradient_data_tot_area<<Area_grad_finite[3].x<<" "<<Area_grad_finite[3].y<<" "<<Area_grad_finite[3].z<<" \n";
    Gradient_data_tot_area<<Area_grad_th[3].x<<" "<<Area_grad_th[3].y<<" "<<Area_grad_th[3].z<<" \n";
    
    
    Gradient_data_tot_area.close();
    }

    if(Test_Bead){
    filename = basic_name+"Bead_Gradient_evaluation_"+std::to_string(0) + ".txt";
    Gradient_data_bead.open(filename);
    Gradient_data_bead<< "Bead grad first finite then theory\n";
    // double H_bar=sqrt(4*PI/A_bar)*c0/2.0;
 
    VertexData<Vector3> Bead_grad_finite = M3DG.Grad_Bead(Gradient_data_bead,false,false);
    VertexData<Vector3> Bead_grad_th = M3DG.Bead_1.Gradient();
 
    for(int i=0;i<7; i++){
    Gradient_data_bead<<Bead_grad_finite[i].x<<" "<<Bead_grad_finite[i].y<<" "<<Bead_grad_finite[i].z<<" \n";
    Gradient_data_bead<<Bead_grad_th[i].x<<" "<<Bead_grad_th[i].y<<" "<<Bead_grad_th[i].z<<" \n\n";
    }

    Gradient_data_bead.close();
    }






    // std::cout<<" This is where we will observe the problem in detail\n\n\n";


    // Vector3 r_ij = M3DG.Bead_1.Pos-geometry->inputVertexPositions[3];
    // double r_dist = r_ij.norm();
    // Vector3 unit_r=r_ij.unit();
    // std::cout<< "This two values should be equal "<< r_dist<< "and "<< dot(r_ij,unit_r)<<" which means the unit operation is correct\n";
    
    // double dual_area = geometry->circumcentricDualArea(mesh->vertex(3));
    // double strength= 1.0;
    // double sigma= 1.0;

    // double E_v = 4*strength*(pow(sigma/r_dist,12)-pow(sigma/r_dist,6));

    // Vector3 F_inter=(4*strength*( -12*pow(sigma/r_dist,12)+6*pow(sigma/r_dist,6))/r_dist)*unit_r;

    // Vector3 F_area_grad=-2*geometry->vertexNormalMeanCurvature(mesh->vertex(3));

    // Vector3 Gradient_finite_dif=Bead_grad_finite[3];


    // std::cout<<"The distance between them is "<<r_dist<<"\n";
    

    // std::cout<<"The Finite difference says "<< Gradient_finite_dif.x<<" "<< Gradient_finite_dif.y <<" "<< Gradient_finite_dif.z<<"\n";


    // std::cout<< "The interaction force says "<< F_inter.x<<" "<< F_inter.y<<" "<< F_inter.z <<"\n";

    // std::cout<< "The Area grad says "<< F_area_grad.x << " "<< F_area_grad.y <<" " << F_area_grad.z<<" \n";

    // std::cout<<"The dual area is "<< dual_area<< " and the interaction energy is "<< E_v<<" \n\n";

    // std::cout<< "The interaction force*dual_area says "<< dual_area*F_inter.x<<" "<< dual_area*F_inter.y<<" "<< dual_area*F_inter.z <<"\n";

    // std::cout<< "The Area grad*E_v says "<< E_v*F_area_grad.x << " "<< E_v*F_area_grad.y <<" " << E_v*F_area_grad.z<<" \n";
    

    // std::cout<<"The Finite difference says "<< Gradient_finite_dif.x<<" "<< Gradient_finite_dif.y <<" "<< Gradient_finite_dif.z<<"\n\n";

    // std::cout<<"If we add them we get "<< E_v*F_area_grad.x+dual_area*F_inter.x << " "<< E_v*F_area_grad.y+dual_area*F_inter.y<<" "<< E_v*F_area_grad.z+dual_area*F_inter.z<<" \n\n";




    // std::cout<<"If we assume that the contribution in the radial distance is correct we can get a new vector \n\n";

    // Vector3 Area_contribution= Gradient_finite_dif-dual_area*F_inter;



    // std::cout<<"If the E_v is correctly calculated then we can divide by E_v as long as it is not 0\n Then we will get the direction of the change\n This change is "<< Area_contribution.x/E_v<< " "<< Area_contribution.y/E_v<<" " << Area_contribution.z/E_v <<" \n";









    









    // for(size_t current_t=0; current_t <1000;current_t++){
    // if(current_t==0){
    // Sim_data.open(filename_basic);
    // }
    // else{
    // Sim_data.open(filename_basic,std::ios_base::app);
    // }

    // if(true){
    //     // if(current_t%10==0 ){
    //     n_vert_old=mesh->nVertices();
    //     n_vert_new=1;

    //     counter=0;
    //     while(n_vert_new!=n_vert_old && counter<10){
    //     // std::cout<<"Remeshing\n";
    //     n_vert_old=n_vert_new;
    //     RemeshOptions Options;
    //     Options.targetEdgeLength=trgt_len;
    //     Options.curvatureAdaptation=Curv_adap;
    //     Options.maxIterations=10;
    //     Options.minRelativeLength=Min_rel_length;
    //     Options.smoothStyle=RemeshSmoothStyle::Circumcentric;
    //     Options.boundaryCondition=RemeshBoundaryCondition::Tangential;
    //     MutationManager Mutation_manager(*mesh,*geometry);
    //     remesh(*mesh,*geometry,Mutation_manager,Options);
    //     n_vert_new=mesh->nVertices();
    //     counter=counter+1; 
    //     }
    //     }

    // if(current_t%100==0){
    //         Save_mesh(basic_name,current_t);
    //         Save_bead_data=true;
    //     }

    
    // if(current_t%200==0){
    //     // filename = basic_name+"Vol_Gradient_evaluation_"+std::to_string(current_t) + ".txt";
    //     // // std::ofstream Gradient_data(filename);
    //     // Gradient_data.open(filename);
    //     // // std::ofstream o(basic_name+std::to_string(current_t)+".obj");
    //     // Gradient_data<< "Volume grad\n";
    //     // M3DG.Grad_Vol(Gradient_data,P0,V_bar,true);

    //     // Gradient_data.close();
    //     // double A_bar=4*PI*pow(3*V_bar/(4*PI*nu_evol),2.0/3.0);

    //     // filename = basic_name+"Area_Gradient_evaluation_"+std::to_string(current_t) + ".txt";
    //     // Gradient_data_area.open(filename);
    //     // Gradient_data_area<< "Area grad\n";
    //     // M3DG.Grad_Area(Gradient_data_area,A_bar,KA,true);

    //     // Gradient_data_area.close();

    //     // filename = basic_name+"Bending_Gradient_evaluation_"+std::to_string(current_t) + ".txt";
    //     // Gradient_data_bending.open(filename);
    //     // Gradient_data_bending<< "Bending grad\n";
    //     // double H_bar=sqrt(4*PI/A_bar)*c0/2.0;
    //     // M3DG.Grad_Bending(Gradient_data_bending,H_bar,KB,true);
    //     // Gradient_data_bending.close();


    //     // filename = basic_name+"Bending_Gradient_evaluation_2_"+std::to_string(current_t) + ".txt";
    //     // Gradient_data_bending_norms.open(filename);
    //     // Gradient_data_bending_norms<< "Bending grad norm diff\n";
    //     // // double H_bar=sqrt(4*PI/A_bar)*c0/2.0;
    //     // M3DG.Grad_Bending_2(Gradient_data_bending_norms,H_bar,KB);
    //     // Gradient_data_bending_norms.close();
        


    //     filename = basic_name+"Area_tot_Gradient_evaluation_"+std::to_string(current_t) + ".txt";
    //     Gradient_data_tot_area.open(filename);
    //     Gradient_data_tot_area<<"Area tot grad\n";
    //     M3DG.Grad_tot_Area(Gradient_data_tot_area,true);
    //     Gradient_data_tot_area.close();


    //     filename = basic_name+"Bead_Gradient_evaluation_"+std::to_string(current_t) + ".txt";
    //     Gradient_data_bead.open(filename);
    //     Gradient_data_bead<< "Bead grad\n";
    //     // double H_bar=sqrt(4*PI/A_bar)*c0/2.0;
    //     M3DG.Grad_Bead(Gradient_data_bead,true,false);
    //     Gradient_data_bead.close();


    //     Volume= geometry->totalVolume();
    //     Area=geometry->totalArea();
    //     nu_obs=3*Volume/(4*PI*pow( Area/(4*PI) ,1.5 ));
    //     H0=sqrt(4*PI/Area)*c0/2.0;
        
    //     if(current_t==0){
    //         nu_0=nu_obs;
    //     }
   

    // }
     
    // nu_evol= time<50 ? nu_0 + (nu-nu_0)*time/50 : nu; 
    // if(with_bead){
    //     dt_sim=M3DG.integrate(TS,V_bar,nu_evol,c0,P0,KA,KB,Kd,Sim_data,time,Save_bead_data,Bead_data);
        
    // }
    // else{
    // dt_sim=M3DG.integrate(TS,V_bar,nu_evol,c0,P0,KA,KB,Kd,Sim_data,time);
    // }

    // if(dt_sim==-1)
    //     {
    //     std::cout<<"Sim broke or timestep very small\n";
    //     break;
    //     }
    // else{
    //     time+=dt_sim;
            
    //     }

    

    // Sim_data.close();

    // // filename = basic_name+"Bending_Gradient_analitical_terms.txt";
    // // std::ofstream Gradient_data_bending_2(filename);
    // // Gradient_data_bending_2<< "Bending gradients per component\n";
    // // // double H_bar=sqrt(4*PI/A_bar)*c0/2.0;
    // // M3DG.Bending_test(Gradient_data_bending_2,H_bar,KB);

    // // Gradient_data_bending_2.close();

    
    
    
    // }




    }

    delete mesh;
    delete geometry;
    // return 1;
    return EXIT_SUCCESS;
}



