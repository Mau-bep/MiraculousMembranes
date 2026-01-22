

#include  "ShapeNLP.hpp"
#include <cassert>

#include <Eigen/Core>


#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

using namespace Ipopt;


ShapeNLP::ShapeNLP()
{ }

ShapeNLP::~ShapeNLP()
{ }

void ShapeNLP::update_positions(
    const Number* x
)
{
    // This function can only be called after get_starting_point has been called
    // std::cout<<"Updating positions\n";

    // Now the idea is to update inputvertexpositions
    Vector3 Pos;   
    for(Index i = 0; i < Index(M3DG->mesh->nVertices()); i++){
        // Pos = M3DG->geometry->inputVertexPositions[i];
        Pos.x = x[3*i];
        Pos.y = x[3*i+1];
        Pos.z = x[3*i+2];
        M3DG->geometry->inputVertexPositions[i] = Pos;
    }
    M3DG->geometry->refreshQuantities();
    
    size_t N_vert = M3DG->mesh->nVertices();

    for(Index i = 0; i < Index(M3DG->Beads.size()); i++){
        Pos.x = x[3*(N_vert + i)];
        Pos.y = x[3*(N_vert + i)+1];
        Pos.z = x[3*(N_vert + i)+2];
        M3DG->Beads[i]->Pos = Pos;
    }

    return;
}
bool ShapeNLP::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g,
   Index&          nnz_h_lag,
   IndexStyleEnum& index_style
)
{
    // std::cout<<"Getting NLP info \n";
    double N_vert = M3DG->mesh->nVertices();
    double N_beads = M3DG->Beads.size();

   // The number of variables is 3 times the number of vertices + the beads
   n = 3 * (N_vert + N_beads);
//    std::cout<<"The number of beads is "<< N_beads << "\n";

   // There are 2 constraints in this problem the volume and the area
   m = M3DG->Sim_handler->Constraints.size();

   // 2 constraints, so Jacobian has nonzero entries. For volume and area constraints, each depends on all vertices
   nnz_jac_g = m*(n-3*N_beads); //I need to calculate this
    // The Jacobian will have more terms when we constrain the beads

   
   // The Hessian has a lot of nonzero entries. 
   // There are the conectivity ones first, then we also have the contributions from the beads which i may need to make dense because even tho a bead wont interact with all the vertices it could interarct with many of them 0    
    //So here I will have to calculate the number of nonzeros in the hesssian.    
    Eigen::SparseMatrix<double> Hessian_temp = M3DG->Sim_handler->Calculate_Hessian_E() + M3DG->Sim_handler->H_Volume(std::vector<double>{1.0});
    for(int i = 0; i < M3DG->Beads.size(); i++) Hessian_temp += M3DG->Beads[i]->Bead_I->Hessian_IP();

    nnz_h_lag = (Hessian_temp.nonZeros()+n)/2 ; // Conectivity + Diagonal terms for regularization
    // nnz_h_lag = 144;
    // std::cout<<"The number of nonzeros in the hessian are are " << nnz_h_lag<<"\n";
   // We use the standard fortran index style for row/col entries
   index_style = FORTRAN_STYLE;
    // std::cout<<"NLP info obtained \n";
   return true;
}

bool ShapeNLP::get_bounds_info(
   Index   n,
   Number* x_l,
   Number* x_u,
   Index   m,
   Number* g_l,
   Number* g_u
)
{
   // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
   // If desired, we could assert to make sure they are what we think they are.
//    std::cout<<"\t\tGetting bounds info\n";
   assert(n == 3 * (M3DG->mesh->nVertices() + M3DG->Beads.size()));
   assert(m == M3DG->Sim_handler->Constraints.size());

    // Set bounds for vertex positions (In theory no bounds, but we set large bounds to help numerical stability)
    for (Index i = 0; i < n; i++) {
        x_l[i] = -1.0e19;
        x_u[i] = 1.0e19;
    }
    // Set bounds for the constraints (equality constraints for volume and area)
    g_l[0] = M3DG->Sim_handler->Trgt_vol;
    g_u[0] = M3DG->Sim_handler->Trgt_vol;

    if(m==2){
    g_l[1] = M3DG->Sim_handler->Trgt_area;
    g_u[1] = M3DG->Sim_handler->Trgt_area;
    }
    
    
    // std::cout<<"Bounds info obtained \n";
    return true;

}

bool ShapeNLP::get_starting_point(
   Index   n,
   bool    init_x,
   Number* x,
   bool    init_z, 
   Number* z_L,
   Number* z_U,
   Index   m,
   bool    init_lambda,
   Number* lambda
)
{
    // std::cout<<"Getting starting point\n";
   // Here, we assume we only have starting values for x, if you code
   // your own NLP, you can provide starting values for the others if
   // you wish.
   assert(init_x == true);
   assert(init_z == false);
   assert(init_lambda == false);

   // we initialize x with the current vertex positions and bead positions
   Index idx = 0;
   for (int v =0; v < M3DG->mesh->nVertices(); v++) {
       Vector3 pos = M3DG->geometry->inputVertexPositions[v];
       x[idx++] = pos[0];
       x[idx++] = pos[1];
       x[idx++] = pos[2];
   }
   for (auto bead : M3DG->Beads) {
       x[idx++] = bead->Pos[0];
       x[idx++] = bead->Pos[1];
       x[idx++] = bead->Pos[2];
   }

    // std::cout<<idx <<" variables well set \n";
    // std::cout<<"Starting point obtained \n";
   return true;
}

bool ShapeNLP::eval_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number&       obj_value
)
{
    // std::cout<<"Evaluating F\n";
    if(new_x) update_positions(x);
    
    double E = 0;

    M3DG->Sim_handler->Calculate_energies(&E);
    // std::cout<<"Current energy is "<< E << "\n";    

    obj_value = E;

    // std::cout<<"F evaluated \n";
    return true;
}

bool ShapeNLP::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
)
{
    if(new_x) update_positions(x);
    M3DG->Sim_handler->Calculate_gradient();
    
    VertexData<Vector3> Current_grad = M3DG->Sim_handler->Current_grad;

    Index idx = 0;
    for (int v = 0; v < M3DG->mesh->nVertices(); v++) {
       
        Vector3 grad = -1 * Current_grad[v];
        grad_f[idx++] = grad[0];
        grad_f[idx++] = grad[1];
        grad_f[idx++] = grad[2];
    }
    // Now the beads
    Vector3 beadGrad;
    for (size_t b = 0; b < M3DG->Beads.size(); b++) {
        beadGrad = -1 * M3DG->Beads[b]->Total_force;
        grad_f[idx++] = beadGrad[0];
        grad_f[idx++] = beadGrad[1];
        grad_f[idx++] = beadGrad[2];
    }  



    return true;
}

bool ShapeNLP::eval_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Number*       g
)
{
    // std::cout<<"Evaluating g\n";
    if(new_x) update_positions(x);
    
    // Volume constraint
    g[0] = M3DG->geometry->totalVolume();

    // Area constraint
    if(m==2){
    g[1] = M3DG->geometry->totalArea();
    }
    
    return true;
}

bool ShapeNLP::eval_jac_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Index         nele_jac,
   Index*        iRow,
   Index*        jCol,
   Number*       values
)
{
    // std::cout<<"Evaluating jac g\n";
    
    if(values == NULL)
    {
        int max_row = 0;
        int max_col = 0;
        // std::cout<<"Values are null\n";
        // return the structure of the jacobian of the constraints
        Index idx = 0;
        // iRow[0] = 1;
        // jCol[0] = 1;

        
        for (Index constr = 0; constr < m; constr++) {
            for (int v = 0; v < M3DG->mesh->nVertices(); v++) {
                // Each vertex contributes to each constraint
                for (int dim = 0; dim < 3; dim++) {
                    iRow[idx] = constr +1 ; // +1 for Fortran style
                    jCol[idx] = Index(3 * v + dim + 1); // +1 for Fortran style
                    idx++;
                }
            }

        }


    }
    else
    {
        if(new_x) update_positions(x);

        VertexData<Vector3> volGrad = M3DG->Sim_handler->F_Volume(std::vector<double>{1.0});
    
        Index idx = 0;
        // values[idx] = 1.0;
        for (size_t v = 0;  v < M3DG->mesh->nVertices(); v++) {
            Vector3 vVolGrad = -1 * volGrad[v];
            // Vector3 vVolGrad({0.0,0.0,0.0});
            // Volume constraint
            values[idx++] = Number(vVolGrad[0]);
            values[idx++] = Number(vVolGrad[1]);
            values[idx++] = Number(vVolGrad[2]);
            // Area constraint
        }
        if(m==2){
        VertexData<Vector3> areaGrad = M3DG->Sim_handler->F_SurfaceTension_2(std::vector<double>{1.0});

        for(size_t v = 0;  v < M3DG->mesh->nVertices(); v++){
            Vector3 vAreaGrad = -1 * areaGrad[v];
            values[idx++] = Number(vAreaGrad[0]);
            values[idx++] = Number(vAreaGrad[1]);
            values[idx++] = Number(vAreaGrad[2]);
        }
        }
        // std::cout<<"The value for index when evaluating is "<< idx<<"\n";
    }
    


    return true;

}

bool ShapeNLP::eval_h(
   Index         n,
   const Number* x,
   bool          new_x,
   Number        obj_factor,
   Index         m,
   const Number* lambda,
   bool          new_lambda,
   Index         nele_hess,
   Index*        iRow,
   Index*        jCol,
   Number*       values
)
{
    // std::cout<<"Evaluating h\n";
    // std::cout<<"The number of nonzero elements in the hessian are " << nele_hess << "\n";
    if( values == NULL )
    {   
        // std::cout<<"First eval\n";
        int max_row = 0;
        int max_col = 0;
        // return the structure. This is a symmetric matrix, fill the lower left
        // triangle only.
        // Ok so 
        // std::cout<<"Irow is " << iRow<<" \n";


        // Eigen::SparseMatrix<double> Hessian_temp = M3DG->Sim_handler->Calculate_Hessian_E();
        Eigen::SparseMatrix<double> Hessian_temp;
        
        if(m==1) Hessian_temp = M3DG->Sim_handler->Calculate_Hessian_E() + M3DG->Sim_handler->H_Volume(std::vector<double>{1.0}); 
        if(m==2) Hessian_temp = M3DG->Sim_handler->Calculate_Hessian_E() + M3DG->Sim_handler->H_Volume(std::vector<double>{1.0}) + M3DG->Sim_handler->H_SurfaceTension({1.0});
        
        for(int i = 0; i < M3DG->Beads.size(); i++) Hessian_temp += M3DG->Beads[i]->Bead_I->Hessian_IP();
        
        // Eigen::SparseMatrix<double> Hessian_temp = M3DG->Sim_handler->Calculate_Hessian_E() + M3DG->Sim_handler->H_Volume({1.0}) + M3DG->Sim_handler->H_SurfaceTension({1.0}); 
        Index idx = 0;
        // std::cout<<"Evaluating hessian of energy\n";
        for (int k=0; k<Hessian_temp.outerSize(); ++k){
            for (Eigen::SparseMatrix<double>::InnerIterator it(Hessian_temp,k); it; ++it){
                if(it.row() >= it.col()){ // Lower triangle
                    iRow[idx] = Index(it.row() +1); // Fortran style 
                    jCol[idx] = Index(it.col() +1); // Fortran style

                    idx++;
                }
            }
        }
       
        
    }
    else
    {
        if(new_x) update_positions(x);



        
        Eigen::SparseMatrix<double> Hessian_temp; 
        if(m == 1) Hessian_temp = obj_factor * M3DG->Sim_handler->Calculate_Hessian_E() + lambda[0]*M3DG->Sim_handler->H_Volume(std::vector<double>{1.0});
        if(m == 2) Hessian_temp = obj_factor * M3DG->Sim_handler->Calculate_Hessian_E() + lambda[0]*M3DG->Sim_handler->H_Volume(std::vector<double>{1.0}) + lambda[1]*M3DG->Sim_handler->H_SurfaceTension({1.0});
        for(int i = 0; i < M3DG->Beads.size(); i++) Hessian_temp += obj_factor * M3DG->Beads[i]->Bead_I->Hessian_IP();
        
        int idx = 0;
        for (int k=0; k<Hessian_temp.outerSize(); ++k){
            for (Eigen::SparseMatrix<double>::InnerIterator it(Hessian_temp,k); it; ++it){
                if(it.row() >= it.col()){ // Lower triangle
                    values[idx] = Number(it.value());
                    idx++;
                    // std::cout<<"The row is "<< it.row()+1 << " and the col is "<< it.col()+1 << "\n";
                }
            }
        }
              

    
    }
    
    return true; 
}

void ShapeNLP::finalize_solution(
   SolverReturn               status,
   Index                      n,
   const Number*              x,
   const Number*              z_L,
   const Number*              z_U,
   Index                      m,
   const Number*              g,
   const Number*              lambda,
   Number                     obj_value,
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq
)
{
    // Here is where we would store the solution to variables, or write to a file, etc
    // For this example, we just print the solution to the console.
    std::cout <<"\n\n************************************* \n\n";
    std::cout << "        YEIII WE MADE IT\n\n";
    std::cout <<"************************************* \n";
    
    update_positions(x);
    // std::cout << "\n\nObjective value\n";
    // std::cout << "f(x*) = " << obj_value << "\n";
    // std::cout<<" The vertex positions are :\n";
    // for (int v =0; v < M3DG->mesh->nVertices(); v++) {
    //     Vector3 pos = M3DG->geometry->inputVertexPositions[v];
    //     std::cout << "Vertex " << v << ": (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
    // }

    
    std::cout<<"The target volume is " << M3DG->Sim_handler->Trgt_vol << " and the final volume is "<< M3DG->geometry->totalVolume() << "\n";
    std::cout<<"The target area is " << M3DG->Sim_handler->Trgt_area << " and the final area is "<< M3DG->geometry->totalArea() << "\n";
}


bool ShapeNLP::intermediate_callback(
      AlgorithmMode              mode,
      Index                      iter,
      Number                     obj_value,
      Number                     inf_pr,
      Number                     inf_du,
      Number                     mu,
      Number                     d_norm,
      Number                     regularization_size,
      Number                     alpha_du,
      Number                     alpha_pr,
      Index                      ls_trials,
      const IpoptData*           ip_data,
      IpoptCalculatedQuantities* ip_cq
   ){

    M3DG->Save_mesh(iter); 
    // 
    // if(iter>1){
    // double T_A;
    // for(Face f: M3DG->mesh->faces()){
    //     T_A = M3DG->geometry->faceArea(f);
    //     if(T_A<1e-4){
    //         std::cout<<"Face "<< f.getIndex() << " has area "<< T_A << " at iteration "<< iter << "\n";
    //         return false;
    //     }
    // }
    // }


    return true;
   }