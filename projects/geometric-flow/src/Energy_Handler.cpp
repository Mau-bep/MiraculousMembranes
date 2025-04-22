

#include "Energy_Handler.h"
#include "Beads.h"
#include <fstream>
#include <omp.h>

#include <chrono>
#include <Eigen/Core>


using namespace std;
using namespace geometrycentral;
using namespace geometrycentral::surface;
 

// THings to note, i should make this structure in main and then give a pointer of it to Mem3DG
// Because if i dont then it will be carried all the time by mem3dg and i dont think thats very smart.


E_Handler::E_Handler(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo){
    mesh = inputMesh;
    geometry = inputGeo;

}

E_Handler::E_Handler(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, std::vector<std::string> inputEnergyNames, std::vector<std::vector<double>> inputEnergyConstants ){
    mesh = inputMesh;
    geometry = inputGeo;
    Energies = inputEnergyNames;
    Energy_constants = inputEnergyConstants;
}

void E_Handler::Add_Bead(Bead *bead){

    Beads.push_back(bead);

}

void E_Handler::Add_Energy(std::string Energy_name, std::vector<double> Constants){
    Energies.push_back(Energy_name);
    Energy_constants.push_back(Constants);
}


// SparseMatrix<double> E_Handler::H1_operator(bool CM, bool Vol_c
// The H1 AND H2 OPERATORS WILL BE A TASK FOR LATER 

double E_Handler::E_Volume_constraint(std::vector<double> Constants) const {
double V = geometry->totalVolume();
double KV = Constants[0];
double V_bar = Constants[1];
return 0.5*KV*(V-V_bar)*(V-V_bar)/(V_bar*V_bar);

}

double E_Handler::E_Area_constraint(std::vector<double> Constants) const {
    double KA = Constants[0];
    double A_bar = Constants[1];
    double A = geometry->totalArea();
    // return 0.5*KA*A*A;
    return 0.5*KA*(A-A_bar)*(A-A_bar)/(A_bar*A_bar);
}

double E_Handler::E_SurfaceTension(std::vector<double> Constants) const{
    return geometry->totalArea()*Constants[0];
}

double E_Handler::E_Bending(std::vector<double> Constants) const {
    double KB = Constants[0];
    double H0 = 0.0;
    size_t index;
    double Eb=0;
    double H;
    double r_eff2;
    Vector3 Pos;
    for(Vertex v : mesh->vertices()) {
        //boundary_fix
        // std::cout<<"boundary \n";
        if(v.isBoundary()) continue;
        // std::cout<<" indxe\n";

        index=v.getIndex();
        // std::cout<<" pos\n";
        Pos = geometry->inputVertexPositions[v];
        // std::cout<<"reff \n";
        r_eff2 = Pos.z*Pos.z + Pos.y*Pos.y ;      
        // std::cout<<"boundary? \n";
        if(r_eff2 > 1.6 && boundary ) continue;
        // std::cout<<" MC and dual area\n";
        H=abs(geometry->scalarMeanCurvature(v)/
        geometry->barycentricDualArea(v));
        // std::cout<<" isnan\n";
        if(std::isnan(H)){
          
        //   std::cout<<"Dual area: "<< geometry->barycentricDualArea(v);
        //   std::cout<<"Scalar mean Curv"<< geometry->scalarMeanCurvature(v);
        //   std::cout<<"One of the H is not a number\n";
        continue;
        }        
        // std::cout<<" adding\n";
        Eb+=KB*H*H*geometry->barycentricDualArea(v);
        
        }   
    // std::cout<<" done with ?\n";
    // std::cout<<Eb<<"\n";
    return Eb;

}


VertexData<Vector3> E_Handler::F_Volume_constraint(std::vector<double> Constants) const{
    double KV = Constants[0];
    double V_bar = Constants[1];
    double V = geometry->totalVolume();
    
    size_t index;
    Vector3 Normal;
    size_t N_vert=mesh->nVertices();

    VertexData<Vector3> Force(*mesh);
    for(Vertex v : mesh->vertices()) {
     // do science here
        index=v.getIndex();
        Normal={0,0,0};
        for(Face f : v.adjacentFaces()) {
            Normal+=geometry->faceArea(f)*geometry->faceNormal(f);

        }
        // Force[v.getIndex()]=D_P*Normal/3.0;
        Force[v.getIndex()]=Normal/3.0;
    }




    return -1*KV*((V-V_bar)/(V_bar*V_bar))*Force;
}

VertexData<Vector3> E_Handler::F_SurfaceTension(std::vector<double> Constants) const{

    double sigma = Constants[0];

    size_t index;
    Vector3 Normal;
    size_t N_vert=mesh->nVertices();
    VertexData<Vector3> Force(*mesh);
    Vector3 u;
    Halfedge he_grad;
        
    for(Vertex v : mesh->vertices()) {
            Normal={0,0,0};
            Force[v]={0,0,0};
            Force[v]=-1*2*geometry->vertexNormalMeanCurvature(v);
            
        }

    // std::cout<< "THe surface tension force in magnitude is: "<< -1*lambda*sqrt(Force.transpose()*Force) <<"\n";
    return sigma*Force;

}

VertexData<Vector3> E_Handler::F_Area_constraint(std::vector<double> Constants) const{
    double KA = Constants[0];
    double A_bar = Constants[1];
    double A = geometry->totalArea();
    
    return ((A-A_bar)/(A_bar*A_bar))*F_SurfaceTension({KA});


}

VertexData<Vector3> E_Handler::F_Bending(std::vector<double> Constants) const{
    size_t neigh_index;
    size_t N_vert=mesh->nVertices();

    double KB = Constants[0];
    double H0 = 0.0;

    Vector3 Hij;
    Vector3 Kij;
    Vector3 Sij_1;
    Vector3 Sij_2;

    Vector3 F1={0,0,0};
    Vector3 F2={0,0,0};
    Vector3 F3={0,0,0};
    Vector3 F4={0,0,0};


    Vector3 Position_1;
    Vector3 Position_2;


    VertexData<Vector3> Force(*mesh);




    VertexData<double> Scalar_MC(*mesh,0.0);
    double factor;
    size_t index1;
    for(Vertex v1 : mesh->vertices()) {
        index1=v1.getIndex();
        Scalar_MC[index1]=geometry->scalarMeanCurvature(v1)/geometry->barycentricDualArea(v1);
        }   
    

    
    auto start = chrono::steady_clock::now();
    double H0i;
    double H0j;
    size_t index;
    for(Vertex v: mesh->vertices()){
        F1={0,0,0};
        F2={0,0,0};
        F3={0,0,0};
        F4={0,0,0};
        index=v.getIndex();
        // H0i= (system_time<50? (H_Vector_0[index]+H0)/2.0: H0);
        H0i=H0;
        Position_1=geometry->inputVertexPositions[v];
        for(Halfedge he: v.outgoingHalfedges()){

            
            neigh_index=he.tipVertex().getIndex();
            // H0j= (system_time<50? (H_Vector_0[neigh_index]+H0)/2: H0);
            H0j=H0;

            Position_2=geometry->inputVertexPositions[neigh_index];
      
            


            Kij=geometry->computeHalfedgeGaussianCurvatureVector(he);
            // factor=-1*(Scalar_MC[index]-(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0))-1*(Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
            factor=-1*(Scalar_MC[index]-H0i)-1*(Scalar_MC[neigh_index]-H0j);
            
            F1=F1+factor*Kij;


            Hij=2*geometry->computeHalfedgeMeanCurvatureVector(he);

            // factor=(1/3.0)*(Scalar_MC[index]*Scalar_MC[index] -(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0)*(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0))+(2.0/3.0)*(Scalar_MC[neigh_index]*Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0)*(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
            factor=(1/3.0)*(Scalar_MC[index]*Scalar_MC[index] -H0i*H0i)+(2.0/3.0)*(Scalar_MC[neigh_index]*Scalar_MC[neigh_index]-H0j*H0j);

            F2=F2+factor*Hij;


            Sij_1 =  geometry->edgeLength(he.edge()) * geometry->dihedralAngleGradient(he,he.vertex());

            // factor=-1*(Scalar_MC[index]-(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0));
            factor=-1*(Scalar_MC[index]-H0i);


            F3=F3+factor*Sij_1;

            Sij_2=(geometry->edgeLength(he.twin().edge())*geometry->dihedralAngleGradient(he.twin(),he.vertex())+ geometry->edgeLength(he.next().edge())*geometry->dihedralAngleGradient(he.next(),he.vertex()) + geometry->edgeLength(he.twin().next().next().edge())*geometry->dihedralAngleGradient(he.twin().next().next(), he.vertex()));
            
            // Sij_2=-1*( geometry->cotan(he.next().next())*geometry->faceNormal(he.face()) + geometry->cotan(he.twin())*geometry->faceNormal(he.twin().face()));

            // factor= -1*(Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
            factor= -1*(Scalar_MC[neigh_index]-H0j);
            
            F4=F4+factor*Sij_2;

    }
    



    Force[index]=F1+F2+F3+F4;

    }

  
    return KB*Force;



}


void E_Handler::Calculate_energies(double* E){

    // So this is the calculation of the energies
    // std::cout<<"reassinginn \n";
    *E = 0;

    // std::cout<<"The value of energy is " << *E <<"  \n";
    Energy_values.resize(Energies.size());

    int bead_count = 0;
    // std::cout<<" iterating\n";
    for(size_t i = 0; i < Energies.size(); i++){
        // std::cout<<" The energy is "<< Energies[i]<<"\n";
        if(Energies[i] == "Volume_constraint"){
            Energy_values[i] = E_Volume_constraint(Energy_constants[i]);
            *E += Energy_values[i];
            continue;
        
        }

        if(Energies[i]=="Area_constraint"){
            Energy_values[i] = E_Area_constraint(Energy_constants[i]);
            *E += Energy_values[i];
            continue;
        }

        if(Energies[i]=="Surface_tension" || Energies[i] == "H1_Surface_tension" || Energies[i] == "H2_Surface_tension"){
            Energy_values[i] = E_SurfaceTension(Energy_constants[i]);
            *E += Energy_values[i];
            // std::cout<<"The energy value is " << Energy_values[i]<<" \n";
            // std::cout<<"The value of E is" << *E <<" \n";
            continue;
        }

        if(Energies[i]=="Bending" || Energies[i] == "H1_Bending" || Energies[i]=="H2_Bending" ){
            // std::cout<<"Calculating bending energy \n";
            // std::cout<<"The energy constants are " << Energy_constants[i][0] << " " << Energy_constants[i][1] << "\n";
            // std::cout<<"The size of Energy values is" << Energy_values.size() << "\n";
            Energy_values[i] = E_Bending(Energy_constants[i]);
            // std::cout<<"succesfully bent \n";
            *E += Energy_values[i];
            // std::cout<<"The energy value is " << Energy_values[i]<<" \n";
            // std::cout<<"The value of E is" << *E <<" \n";
            continue;
        }
        if(Energies[i]=="Bead" || Energies[i]=="H1_Bead" || Energies[i]=="H2_Bead")
        {
            Energy_values[i] = Beads[bead_count]->Energy();
            *E += Energy_values[i];
            bead_count++;
            // std::cout<<"The energy value is " << Energy_values[i]<<" \n";
            // std::cout<<"The value of E is" << *E <<" \n";
            continue;
        }
        

    }




    return;
}

void E_Handler::Calculate_gradient(){
    // This function will calculate the gradient of the energy

    // std::cout<<"1 \n";
    Previous_grad = Current_grad;
    // std::cout<<"2 \n";
    Current_grad = VertexData<Vector3>(*mesh);
    // std::cout<<"3 \n";
    VertexData<Vector3> Force_temp;
    int bead_count = 0;
    double grad_norm = 0;
    // std::cout<<"4 \n";
    // std::cout<<"THe size of Energies is " << Energies.size() << "\n";
    // std::cout<<"4 1 \n";
    for(size_t i = 0; i < Energies.size(); i++){
        // std::cout<<"Energy is " << Energies[i]<<" \n";

        if(Energies[i] == "Volume_constraint"){
            Force_temp = F_Volume_constraint(Energy_constants[i]);
            grad_norm = 0.0;
            for(Vertex v : mesh->vertices()){
                grad_norm += Force_temp[v].norm2();
            }
            if(Gradient_norms.size() == i){
                Gradient_norms.push_back(grad_norm);
            
            }
            else{
                Gradient_norms[i] = grad_norm;
            }

            Current_grad+=Force_temp;


            continue;
        }

        if(Energies[i]== "Area_constraint"){
            Force_temp = F_Area_constraint(Energy_constants[i]);
            grad_norm = 0.0;
            for(Vertex v : mesh->vertices()){
                grad_norm += Force_temp[v].norm2();
            }
            if(Gradient_norms.size() == i){
                Gradient_norms.push_back(grad_norm);
            
            }
            else{
                Gradient_norms[i] = grad_norm;
            }

            Current_grad+=Force_temp;


            continue;
        }

        if(Energies[i]=="Surface_tension" ){
            Force_temp = F_SurfaceTension(Energy_constants[i]);
            grad_norm = 0.0;
            for(Vertex v : mesh->vertices()){
                grad_norm += Force_temp[v].norm2();
            }
            if(Gradient_norms.size() == i){
                Gradient_norms.push_back(grad_norm);
            
            }
            else{
                Gradient_norms[i] = grad_norm;
            }

            Current_grad+=Force_temp;

            continue;

            continue;
        }

        if(Energies[i]=="Bending"){
            Force_temp = F_Bending(Energy_constants[i]);
            grad_norm = 0.0;

            for(Vertex v :mesh->vertices()){
                grad_norm += Force_temp[v].norm2();
            }
            if(Gradient_norms.size() == i){
                Gradient_norms.push_back(grad_norm);
            
            }
            else{
                Gradient_norms[i] = grad_norm;
            }

            Current_grad+=Force_temp;


            continue;
        }
        // I need to add the beads 

        if(Energies[i]=="Bead"){
            // std::cout<<"Bead E \n";
            // std::cout<<"The number of beads is " << Beads.size() << "\n";
            Force_temp = Beads[bead_count]->Gradient();
            // std::cout<<"Interactions \n";
            Beads[bead_count]->Bead_interactions();
            // std::cout<<"? \n";
            grad_norm = 0;
            for(size_t j = 0; j < mesh->nVertices(); j++){
             grad_norm+= Force_temp[j].norm2();
            }
            // std::cout<<"not this \n";
            if(Gradient_norms.size() == i){
                Gradient_norms.push_back(grad_norm);
            
            }
            else{
                Gradient_norms[i] = grad_norm;
            }
            Current_grad+=Force_temp;
            bead_count +=1;

            continue;
        }

        }
    // std::cout<<"5 \n";
    grad_norm = 0.0;
    for(size_t i = 0; i < mesh->nVertices(); i++){
        grad_norm+=Current_grad[i].norm2();
    }
    // std::cout<<"6 \n";
    if(Gradient_norms.size() == Energies.size()){
            Gradient_norms.push_back(grad_norm);
        }
        else{
            Gradient_norms[Energies.size()] = grad_norm;
            }

    return;
}

void E_Handler::Do_nothing(){
    // This function does nothing, it is just a placeholder
    return;
}




//   for(size_t i = 0; i < Energies.size(); i++){
    
    
    // if(Energies[i]=="H1_Bead"){
    //   if(!H1_grad){
    //     construction_start = chrono::steady_clock::now();
    //     H1_grad = true;
        
    //     H1_mat = H1_operator(Energy_constants[i][1],Energy_constants[i][2],Energy_constants[i][3]);
    //     construction_end = chrono::steady_clock::now();
    //     time_construct += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();

    //     construction_start = chrono::steady_clock::now();
    //     solverH1.compute(H1_mat);
    //     construction_end = chrono::steady_clock::now();
    //     time_compute += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();
      
    //   }

    //   Energy_vals[i] = Beads[bead_count]->Energy(); 
    //   // std::cout<<"Bead count is" << bead_count <<"\n";
    //   Force_temp = Beads[bead_count]->Gradient();
    //   Beads[bead_count]->Bead_interactions();

    //   solve_start = chrono::steady_clock::now();
    //   Vector<double> RHS = Vector<double>(N_vert*3+N_constraints);
    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     RHS[3*index] = Force_temp[index].x;
    //     RHS[3*index+1] = Force_temp[index].y;
    //     RHS[3*index+2] = Force_temp[index].z;
    //   }
    //   for(int index = 0; index < N_constraints; index++){
    //     RHS[3*N_vert+index] = 0.0;
    //   }
      
    //   Vector<double> Solution = solverH1.solve(RHS);
    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     Force_temp[index] = Vector3({Solution.coeff(3*index),Solution.coeff(3*index+1),Solution.coeff(3*index+2)});
    //     // std::cout<< Force_temp[index].norm2() << " ";
    //   }
    //   solve_end = chrono::steady_clock::now();
    //   time_solve += std::chrono::duration_cast<std::chrono::milliseconds>(solve_end - solve_start).count();


    //   grad_value = 0;
    //   for(size_t j = 0; j < mesh->nVertices(); j++){
    //     grad_value+= Force_temp[j].norm2();
    //   }
    //   Gradient_norms.push_back(grad_value);
    //   Force+=Force_temp;
    //   bead_count +=1;
    //   continue;

    // }
    // if(Energies[i]=="H2_Bead"){
      
    //   if(!H2_grad){
    //     construction_start = chrono::steady_clock::now();
    //     H2_grad = true;
    //     N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];

    //     H2_mat = H2_operator(Energy_constants[i][1],Energy_constants[i][2],Energy_constants[i][3]);
    //     construction_end = chrono::steady_clock::now();
    //     time_construct += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();

    //     construction_start = chrono::steady_clock::now();
    //     solverH2.compute(H2_mat);
    //     construction_end = chrono::steady_clock::now();
    //     time_compute += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();
      
    //   }

    //   Energy_vals[i] = Beads[bead_count]->Energy(); 
    //   // std::cout<<"Bead count is" << bead_count <<"\n";
    //   Force_temp = Beads[bead_count]->Gradient();
    //   Beads[bead_count]->Bead_interactions();

    //   solve_start = chrono::steady_clock::now();
    //   Vector<double> RHS = Vector<double>(N_vert*3+N_constraints);
    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     RHS[3*index] = Force_temp[index].x;
    //     RHS[3*index+1] = Force_temp[index].y;
    //     RHS[3*index+2] = Force_temp[index].z;
    //   }
    //   for(int index = 0; index < N_constraints; index++){
    //     RHS[3*N_vert+index] = 0.0;
    //   }
      
    //   Vector<double> Solution = solverH2.solve(RHS);
    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     Force_temp[index] = Vector3({Solution.coeff(3*index),Solution.coeff(3*index+1),Solution.coeff(3*index+2)});
    //     // std::cout<< Force_temp[index].norm2() << " ";
    //   }
    //   solve_end = chrono::steady_clock::now();
    //   time_solve += std::chrono::duration_cast<std::chrono::milliseconds>(solve_end - solve_start).count();


    //   grad_value = 0;
    //   for(size_t j = 0; j < mesh->nVertices(); j++){
    //     grad_value+= Force_temp[j].norm2();
    //   }
    //   Gradient_norms.push_back(grad_value);
    //   Force+=Force_temp;
    //   bead_count +=1;
    //   // std::cout<<"Here\n";
    //   continue;


    // }

    // if(Energies[i]=="H2_Bending"){
    //   // std::cout<<"THe beding energy is hereee\n";
    //   if(!H2_grad){
    //     construction_start = chrono::steady_clock::now();
    //     H2_grad = true;
    //     N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];

    //     H2_mat = H2_operator(Energy_constants[i][1],Energy_constants[i][2],Energy_constants[i][3]);
    //     construction_end = chrono::steady_clock::now();
    //     time_construct += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();

    //     construction_start = chrono::steady_clock::now();
    //     solverH2.compute(H2_mat);
    //     construction_end = chrono::steady_clock::now();
    //     time_compute += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();
    //     // std::cout<<"Solving possible\n";
    //   }

    //   double KB = Energy_constants[i][0];
    //   Energy_vals[i] = E_Bending(0.0,KB);
    //   // Now i need to compute the original gradient
    //   Force_temp = KB*Bending(0.0);

    //   // I would say the solve considers the construction?

    //   N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];
    //   solve_start = chrono::steady_clock::now();
    //   Vector<double> RHS = Vector<double>(N_vert*3+N_constraints);
    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     RHS[3*index] = Force_temp[index].x;
    //     RHS[3*index+1] = Force_temp[index].y;
    //     RHS[3*index+2] = Force_temp[index].z;
    //     // std::cout<<"Force temp has magnitude" << Force_temp[index].norm2() <<"and barycentric " << Barycentric_area[index]<<"\n";
    //   }
    //   // std::cout<<"THe number of constraints is " << N_constraints<<" duh\n";
    //   for(int index = 0; index < N_constraints; index++){
    //     RHS[3*N_vert+index] = 0.0;
    //   }
    //   // std::cout<<"The RHS is readyy\n";
    //   // std::cout<<"Some components are " << RHS[0] <<" "<< RHS[1] <<" "<< RHS[2] <<" \n";
    //   // The RHS is readyy   
    //   // std::cout<<"Solving RHS\n";
      
    //   Vector<double> Solution = solverH2.solve(RHS);
    //   // std::cout<<"Solved\n";
    //   // std::cout<<"The solution is ready\n";
    //   // std::cout<<"Some components of the solution are " << Solution[0] <<" "<< Solution[1] <<" "<< Solution[2] <<" \n";
    //   // std::cout<<"The forces are ";
    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     Force_temp[index] = Vector3({Solution.coeff(3*index),Solution.coeff(3*index+1),Solution.coeff(3*index+2)});
    //     // std::cout<< Force_temp[index].norm2() << " ";
    //   }
    //   solve_end = chrono::steady_clock::now();
    //   time_solve += std::chrono::duration_cast<std::chrono::milliseconds>(solve_end - solve_start).count();
      
    //   // SO i have calculated the new force
    //   grad_value = 0;
    //   for(size_t j = 0; j < mesh->nVertices(); j++){
    //     grad_value+= Force_temp[j].norm2();
    //   }
    //   // std::cout<<"THe grad value is " << grad_value <<" \n";
    //   Gradient_norms.push_back(grad_value);
    //   Force+=Force_temp;

    //   continue;

    // }
    // if(Energies[i]=="H1_Bending"){
    //   // std::cout<<"THe beding energy is hereee\n";
    //   if(!H1_grad){
    //     construction_start = chrono::steady_clock::now();
    //     H1_grad = true;
    //     N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];

    //     // I need to build the sobolev operator
    //     // H2_mat = SparseMatrix<double>(N_vert*3+N_constraints,N_vert*3+N_constraints);
    //     H1_mat = H1_operator(Energy_constants[i][1],Energy_constants[i][2],Energy_constants[i][3]);
    //     construction_end = chrono::steady_clock::now();
    //     time_construct += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();

    //     // std::cout<<"THe solver is being solved\n"; 
    //     // std::cout<<"THe matrix has some nonzero components" << H2_mat.nonZeros() <<"\n";
    //     construction_start = chrono::steady_clock::now();
    //     solverH1.compute(H1_mat);
    //     construction_end = chrono::steady_clock::now();
    //     time_compute += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();
    //     // std::cout<<"Solving possible\n";
    //   }

    //   double KB = Energy_constants[i][0];
    //   Energy_vals[i] = E_Bending(0.0,KB);
    //   // Now i need to compute the original gradient
    //   Force_temp = KB*Bending(0.0);

    //   // double temp_grad = 0;
    //   // for(size_t j = 0; j < mesh->nVertices(); j++){
    //   //   temp_grad += Force_temp[j].norm2();
    //   // }
    //   // std::cout<<"Temp grad is " << temp_grad <<" \n";
    //   //Now i need to add the constraint

    //   // I would say the solve considers the construction?

    //   N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];
    //   solve_start = chrono::steady_clock::now();
    //   Vector<double> RHS = Vector<double>(N_vert*3+N_constraints);
    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     RHS[3*index] = Force_temp[index].x;
    //     RHS[3*index+1] = Force_temp[index].y;
    //     RHS[3*index+2] = Force_temp[index].z;
    //     // std::cout<<"Force temp has magnitude" << Force_temp[index].norm2() <<"and barycentric " << Barycentric_area[index]<<"\n";
    //   }
    //   // std::cout<<"THe number of constraints is " << N_constraints<<" duh\n";
    //   for(int index = 0; index < N_constraints; index++){
    //     RHS[3*N_vert+index] = 0.0;
    //   }
    //   // std::cout<<"The RHS is readyy\n";
    //   // std::cout<<"Some components are " << RHS[0] <<" "<< RHS[1] <<" "<< RHS[2] <<" \n";
    //   // The RHS is readyy   
    //   // std::cout<<"Solving RHS\n";
      
    //   Vector<double> Solution = solverH1.solve(RHS);
    //   // std::cout<<"Solved\n";
    //   // std::cout<<"The solution is ready\n";
    //   // std::cout<<"Some components of the solution are " << Solution[0] <<" "<< Solution[1] <<" "<< Solution[2] <<" \n";
    //   // std::cout<<"The forces are ";
    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     Force_temp[index] = Vector3({Solution.coeff(3*index),Solution.coeff(3*index+1),Solution.coeff(3*index+2)});
    //     // std::cout<< Force_temp[index].norm2() << " ";
    //   }
    //   solve_end = chrono::steady_clock::now();
    //   time_solve += std::chrono::duration_cast<std::chrono::milliseconds>(solve_end - solve_start).count();
      
    //   // SO i have calculated the new force
    //   grad_value = 0;
    //   for(size_t j = 0; j < mesh->nVertices(); j++){
    //     grad_value+= Force_temp[j].norm2();
    //   }
    //   // std::cout<<"THe grad value is " << grad_value <<" \n";
    //   Gradient_norms.push_back(grad_value);
    //   Force+=Force_temp;

    //   continue;

    // }
    // if(Energies[i]=="H2_Surface_tension"){
    //   if(!H2_grad){
    //     construction_start = chrono::steady_clock::now();
        
    //     H2_grad = true;
    //     N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];

    //     // I need to build the sobolev operator
    //     // H2_mat = SparseMatrix<double>(N_vert*3+N_constraints,N_vert*3+N_constraints);
    //     H2_mat = H2_operator(Energy_constants[i][1],Energy_constants[i][2],Energy_constants[i][3]);

    //     construction_end = chrono::steady_clock::now();
    //     time_construct += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();

    //     construction_start = chrono::steady_clock::now();
    //     solverH2.compute(H2_mat);
    //     construction_end = chrono::steady_clock::now();
    //     time_compute += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();
        
    //     // std::cout<<"Solving possible\n";
    //   }
    //   double sigma = Energy_constants[i][0];
    //   A = geometry->totalArea();
    //   Energy_vals[i] = A*sigma;
    //   Force_temp = sigma*SurfaceGrad();

    //   // I need to multiply by the mass dont I?


    //   N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];

    //   Vector<double> RHS = Vector<double>(N_vert*3+N_constraints);
    //   for (size_t i = 0; i < mesh->nVertices(); i++)
    //   {
    //     RHS[3*i] = Force_temp[i].x;
    //     RHS[3*i+1] = Force_temp[i].y;
    //     RHS[3*i+2] = Force_temp[i].z;
    //   }
    //   for(int index = 0; index < N_constraints; index++){
    //     RHS[3*N_vert+index] = 0.0;
    //   }

    //   solve_start = chrono::steady_clock::now();
    //   Vector<double> Solution = solverH2.solve(RHS);
    //   solve_end = chrono::steady_clock::now();
    //   time_solve += std::chrono::duration_cast<std::chrono::milliseconds>(solve_end - solve_start).count();

    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     Force_temp[index] = Vector3({Solution.coeff(3*index),Solution.coeff(3*index+1),Solution.coeff(3*index+2)});
    //   }
      
      
      
    //   grad_value  = 0.0;
    //   for(size_t j = 0; j < mesh->nVertices(); j++){
    //     grad_value+= Force_temp[j].norm2();
    //   }
    //   Gradient_norms.push_back(grad_value);
    //   Force+=Force_temp;
    //   continue;

    // }

    // if(Energies[i]=="H1_Surface_tension"){
    //   if(!H1_grad){
    //     construction_start = chrono::steady_clock::now();
        
    //     H1_grad = true;
    //     N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];

    //     // I need to build the sobolev operator
    //     // H2_mat = SparseMatrix<double>(N_vert*3+N_constraints,N_vert*3+N_constraints);
    //     H1_mat = H1_operator(Energy_constants[i][1],Energy_constants[i][2],Energy_constants[i][3]);

    //     construction_end = chrono::steady_clock::now();
    //     time_construct += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();

    //     construction_start = chrono::steady_clock::now();
    //     solverH1.compute(H1_mat);
    //     construction_end = chrono::steady_clock::now();
    //     time_compute += std::chrono::duration_cast<std::chrono::milliseconds>(construction_end - construction_start).count();
        
    //     // std::cout<<"Solving possible\n";
    //   }
    //   double sigma = Energy_constants[i][0];
    //   A = geometry->totalArea();
    //   Energy_vals[i] = A*sigma;
    //   Force_temp = sigma*SurfaceGrad();

    //   // I need to multiply by the mass dont I?


    //   N_constraints = 3*Energy_constants[i][1]+Energy_constants[i][2]+Energy_constants[i][3];

    //   Vector<double> RHS = Vector<double>(N_vert*3+N_constraints);
    //   for (size_t i = 0; i < mesh->nVertices(); i++)
    //   {
    //     RHS[3*i] = Barycentric_area[i]*Force_temp[i].x;
    //     RHS[3*i+1] = Barycentric_area[i]*Force_temp[i].y;
    //     RHS[3*i+2] = Barycentric_area[i]*Force_temp[i].z;
    //   }
    //   for(int index = 0; index < N_constraints; index++){
    //     RHS[3*N_vert+index] = 0.0;
    //   }

    //   solve_start = chrono::steady_clock::now();
    //   Vector<double> Solution = solverH1.solve(RHS);
    //   solve_end = chrono::steady_clock::now();
    //   time_solve += std::chrono::duration_cast<std::chrono::milliseconds>(solve_end - solve_start).count();

    //   for(size_t index = 0; index < mesh->nVertices(); index++){
    //     Force_temp[index] = Vector3({Solution.coeff(3*index),Solution.coeff(3*index+1),Solution.coeff(3*index+2)});
    //   }
      
      
      
    //   grad_value  = 0.0;
    //   for(size_t j = 0; j < mesh->nVertices(); j++){
    //     grad_value+= Force_temp[j].norm2();
    //   }
    //   Gradient_norms.push_back(grad_value);
    //   Force+=Force_temp;
    //   continue;

    // }

//   }