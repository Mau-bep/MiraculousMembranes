

#include "NormalNLP.hpp"
#include <cassert>
#include <fstream>

#include <Eigen/Core>

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

using namespace Ipopt;

NormalNLP::NormalNLP()
{
}

NormalNLP::~NormalNLP()
{
}

void NormalNLP::update_positions(
    const Number *x)
{
    // This function can only be called after get_starting_point has been called
    // std::cout << "Updating positions\n";

    std::ofstream file("T_evol.txt", std::ios::app);
    if (file.is_open())
    {
        Index p = M3DG->mesh->nVertices();
        for (Index i = 0; i < p; i++)
        {
            file << x[i] << " ";
        }
        file << "\n";
        file.close();
    }

    VertexData<double> ts(*M3DG->mesh);

    for (Vertex v : M3DG->mesh->vertices())
    {
        ts[v] = x[v.getIndex()];
    }

    // geometry->inputVertexPositions = ORIG_VPOS + Normals*ts;
    M3DG->geometry->inputVertexPositions = ORIG_VPOS + Normals * ts;
    // geometry->refreshQuantities();
    M3DG->geometry->refreshQuantities();
    // return false;
    // std::cout << "The bead should go to " << x[M3DG->mesh->nVertices()] << " " << x[M3DG->mesh->nVertices() + 1] << " " << x[M3DG->mesh->nVertices() + 2] << " \n";

    return;
}
bool NormalNLP::get_nlp_info(
    Index &n,
    Index &m,
    Index &nnz_jac_g,
    Index &nnz_h_lag,
    IndexStyleEnum &index_style)
{
    // std::cout<<"Getting NLP info \n";
    double N_vert = M3DG->mesh->nVertices();
    // double N_beads = M3DG->Beads.size();
    // std::cout << "The number of beads is " << N_beads << "\n";
    // The number of variables is 3 times the number of vertices + the beads
    n = N_vert;
    // std::cout << "The number of vertices is " << N_vert << "\n";

    // There are 2 constraints in this problem the volume and the area

    // std::cout << "The number of constraints is " << M3DG->Sim_handler->Constraints.size() << "\n";
    m = M3DG->Sim_handler->Constraints.size();

    // 2 constraints, so Jacobian has nonzero entries. For volume and area constraints, each depends on all vertices
    nnz_jac_g = m * (N_vert); // I need to calculate this
                              //  The Jacobian will have more terms when we constrain the beads

    // std::cout << "CALCULATING HESSIAN\n";

    Eigen::SparseMatrix<double> Hessian_temp = M3DG->Sim_handler->Calculate_Hessian_E_Normal() + M3DG->Sim_handler->Calculate_Hessian_Constraints_Normal();

    // Here i want to iterate over this
    std::ofstream file("../Results/Nonzeroelements_init.txt", std::ios::app);
    Index idex = 0;
    if (file.is_open())
    {
        // Here i want to iterate over the nonzeros.
        for (int k = 0; k < Hessian_temp.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(Hessian_temp, k); it; ++it)
            {
                // file.write( std::to_string(it.row())+ ","+ std::to_string(it.col())+" ");
                file << it.row() << "," << it.col() << " ";
                if (it.row() <= it.col())
                    idex++;
            }
        }
    }
    file << "\n";

    file.close();

    nnz_h_lag = (Hessian_temp.nonZeros() + n) / 2; // Conectivity + Diagonal terms for regularization

    // nnz_h_lag = 144;
    std::cout << "The number of nonzeros in the hessian are are " << nnz_h_lag << "\n";
    // We use the standard fortran index style for row/col entries
    index_style = FORTRAN_STYLE;
    std::cout << "NLP info obtained \n";
    return true;
}

bool NormalNLP::get_bounds_info(
    Index n,
    Number *x_l,
    Number *x_u,
    Index m,
    Number *g_l,
    Number *g_u)
{
    // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
    // If desired, we could assert to make sure they are what we think they are.
    std::cout << "\t\tGetting bounds info\n";
    assert(n == (M3DG->mesh->nVertices()));
    assert(m == M3DG->Sim_handler->Constraints.size());

    // Set bounds for vertex positions (In theory no bounds, but we set large bounds to help numerical stability)
    for (Index i = 0; i < n; i++)
    {
        x_l[i] = -1.0e19;
        x_u[i] = 1.0e19;
    }
    // Set bounds for the constraints (equality constraints for volume and area)
    g_l[0] = M3DG->Sim_handler->Trgt_vol;
    g_u[0] = M3DG->Sim_handler->Trgt_vol;

    if (m == 2)
    {
        g_l[1] = M3DG->Sim_handler->Trgt_area;
        g_u[1] = M3DG->Sim_handler->Trgt_area;
    }

    std::cout << "Bounds info obtained \n";
    return true;
}

bool NormalNLP::get_starting_point(
    Index n,
    bool init_x,
    Number *x,
    bool init_z,
    Number *z_L,
    Number *z_U,
    Index m,
    bool init_lambda,
    Number *lambda)
{
    std::cout << "Getting starting point\n";
    // Here, we assume we only have starting values for x, if you code
    // your own NLP, you can provide starting values for the others if
    // you wish.
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    // we initialize x with the current vertex positions and bead positions
    Index idx = 0;

    for (int v = 0; v < M3DG->mesh->nVertices(); v++)
    {
        x[idx++] = 0;
    }
    // for (auto bead : M3DG->Beads)
    // {
    //     x[idx++] = bead->Pos[0];
    //     x[idx++] = bead->Pos[1];
    //     x[idx++] = bead->Pos[2];
    // }

    return true;
}

bool NormalNLP::eval_f(
    Index n,
    const Number *x,
    bool new_x,
    Number &obj_value)
{
    // std::cout << "Evaluating F\n";
    // if (new_x)
    update_positions(x);

    double E = 0;

    M3DG->Sim_handler->Calculate_energies(&E);
    // std::cout<<"Current energy is "<< E << "\n";

    obj_value = E;

    // std::cout << "F evaluated \n";
    return true;
}

bool NormalNLP::eval_grad_f(
    Index n,
    const Number *x,
    bool new_x,
    Number *grad_f)
{
    // std::cout << " Eval grad f\n";
    // if (new_x)
    update_positions(x);
    M3DG->Sim_handler->Calculate_gradient();

    VertexData<Vector3> Current_grad = M3DG->Sim_handler->Current_grad;

    Index idx = 0;
    Vector3 beadGrad;
    for (int i = 0; i < M3DG->mesh->nVertices(); i++)
    {
        grad_f[i] = -1 * dot(Current_grad[i], Normals[i]); // Start with zero gradient
        idx++;
    }
    // for (size_t b = 0; b < M3DG->Beads.size(); b++)
    // {
    //     std::cout << "Evaluating gradient for bead " << b << "\n";
    //     beadGrad = -1 * M3DG->Beads[b]->Total_force;
    //     grad_f[idx++] = beadGrad[0];
    //     grad_f[idx++] = beadGrad[1];
    //     grad_f[idx++] = beadGrad[2];
    // }
    // std::cout << "Gradient evaluated \n";
    return true;
}

bool NormalNLP::eval_g(
    Index n,
    const Number *x,
    bool new_x,
    Index m,
    Number *g)
{
    // std::cout << "Evaluating g\n";
    // if (new_x)
    update_positions(x);

    // Volume constraint
    g[0] = M3DG->geometry->totalVolume();

    // Area constraint
    if (m == 2)
    {
        g[1] = M3DG->geometry->totalArea();
    }

    return true;
}

bool NormalNLP::eval_jac_g(
    Index n,
    const Number *x,
    bool new_x,
    Index m,
    Index nele_jac,
    Index *iRow,
    Index *jCol,
    Number *values)
{
    // std::cout << "Evaluating jac g\n";

    if (values == NULL)
    {
        int max_row = 0;
        int max_col = 0;
        // std::cout<<"Values are null\n";
        // return the structure of the jacobian of the constraints
        Index idx = 0;
        // iRow[0] = 1;
        // jCol[0] = 1;

        for (Index constr = 0; constr < m; constr++)
        {
            for (int v = 0; v < M3DG->mesh->nVertices(); v++)
            {
                // Each vertex contributes to each constraint
                // for (int dim = 0; dim < 3; dim++) {
                iRow[idx] = constr + 1;   // +1 for Fortran style
                jCol[idx] = Index(v + 1); // +1 for Fortran style
                idx++;
                // }
            }
        }
    }
    else
    {
        // if (new_x)
        update_positions(x);

        VertexData<Vector3> volGrad = M3DG->Sim_handler->F_Volume(std::vector<double>{1.0});

        Index idx = 0;

        for (size_t i = 0; i < M3DG->mesh->nVertices(); i++)
        {
            values[idx++] = -1 * dot(volGrad[i], Normals[i]);
        }

        if (m == 2)
        {
            VertexData<Vector3> areaGrad = M3DG->Sim_handler->F_SurfaceTension_2(std::vector<double>{1.0});

            for (size_t i = 0; i < M3DG->mesh->nVertices(); i++)
            {
                values[idx++] = -1 * dot(areaGrad[i], Normals[i]);
            }
        }
        // std::cout<<"The value for index when evaluating is "<< idx<<"\n";
    }

    return true;
}

bool NormalNLP::eval_h(
    Index n,
    const Number *x,
    bool new_x,
    Number obj_factor,
    Index m,
    const Number *lambda,
    bool new_lambda,
    Index nele_hess,
    Index *iRow,
    Index *jCol,
    Number *values)
{
    // std::cout << "Evaluating h\n";
    // std::cout << "The number of nonzero elements in the hessian are " << nele_hess << "\n";
    if (values == NULL)
    {
        Eigen::SparseMatrix<double> Hessian_temp;
        Hessian_temp = M3DG->Sim_handler->Calculate_Hessian_E_Normal() + M3DG->Sim_handler->Calculate_Hessian_Constraints_Normal();
        // std::cout << "The number of nonzero elements in the hessian are " << Hessian_temp.nonZeros() << "\n";
        Index idx = 0;

        std::ofstream file("../Results/Nonzeroelements_first_eval.txt", std::ios::app);

        for (int k = 0; k < Hessian_temp.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(Hessian_temp, k); it; ++it)
            {
                file << it.row() << "," << it.col() << " ";
                if (it.row() >= it.col())
                { // Lower triangle
                    // std::cout << "Adding nonzero entry to Hessian structure at row " << it.row() + 1 << " and column " << it.col() + 1 << "\n";
                    iRow[idx] = Index(it.row() + 1);
                    jCol[idx] = Index(it.col() + 1);
                    idx++;
                }
            }
        }

        file << "\n";
        file.close();
        // std::cout << "Hessian structure obtained with " << idx << " nonzero entries\n";
    }
    else
    {
        // if (new_x)
        update_positions(x);

        Eigen::SparseMatrix<double> Hessian_temp;
        if (m == 1)
            Hessian_temp = obj_factor * M3DG->Sim_handler->Calculate_Hessian_E_Normal() + lambda[0] * M3DG->Sim_handler->H_Volume_Normal(std::vector<double>{1.0});
        if (m == 2)
            Hessian_temp = obj_factor * M3DG->Sim_handler->Calculate_Hessian_E_Normal() + lambda[0] * M3DG->Sim_handler->H_Volume_Normal(std::vector<double>{1.0}) + lambda[1] * M3DG->Sim_handler->H_SurfaceTension_Normal({1.0});

        // std::cout << "The number of nonzero elements in the hessian are " << Hessian_temp.nonZeros() << "\n";
        int idx = 0;
        std::ofstream file("../Results/Nonzeroelements_real_eval.txt", std::ios::app);

        for (int k = 0; k < Hessian_temp.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(Hessian_temp, k); it; ++it)
            {
                file << it.row() << "," << it.col() << " ";
                if (it.row() >= it.col())
                { // Lower triangle
                    // std::cout << "Adding nonzero entry to Hessian structure at row " << it.row() + 1 << " and column " << it.col() + 1 << "\n";
                    values[idx] = Number(it.value());
                    idx++;
                }
            }
        }
        file << "\n";
        file.close();
        // std::cout << "The number of nonzero entries added to the hessian values is " << idx << "\n";
        // if(idx != nele_hess) std::cout<<"There is a difference"
        if (idx < nele_hess)
        {
            std::cout << "Warning: The number of nonzero entries in the Hessian (" << idx << ") is less than the expected number (" << nele_hess << "). This may indicate an issue with the Hessian structure or values.\n";
            // std::cout << "Filling remaining entries with zero.\n";
            for (int val = idx; val < nele_hess; val++)
            {
                // std::cout << "Filling entry " << val << " with zero.\n";
                // std::cout << "The size of rows is " << sizeof(iRow) / sizeof(iRow[0]) << " and the size of cols is " << sizeof(jCol) / sizeof(jCol[0]) << "\n";
                values[val] = 0.0; // Fill remaining entries with zero
                // std::cout << "Succesfully done\n";
                // std::cout << "THis value corresponsds to row " << iRow[val - 1] << " and column " << jCol[val - 1] << "\n";
            }
        }
    }
    // std::cout << "Done with this particular op\n";
    return true;
}

void NormalNLP::finalize_solution(
    SolverReturn status,
    Index n,
    const Number *x,
    const Number *z_L,
    const Number *z_U,
    Index m,
    const Number *g,
    const Number *lambda,
    Number obj_value,
    const IpoptData *ip_data,
    IpoptCalculatedQuantities *ip_cq)
{
    // Here is where we would store the solution to variables, or write to a file, etc
    // For this example, we just print the solution to the console.
    std::cout << "\n\n"
              << "************************************* \n\n";
    std::cout << "        YEIII WE MADE IT\n\n";
    std::cout << "************************************* \n";

    update_positions(x);
    std::cout << "\n\nObjective value\n";
    std::cout << "f(x*) = " << obj_value << "\n";
    // std::cout<<" The vertex positions are :\n";
    // for (int v =0; v < M3DG->mesh->nVertices(); v++) {
    //     Vector3 pos = M3DG->geometry->inputVertexPositions[v];
    //     std::cout << "Vertex " << v << ": (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
    // }

    std::cout << "The target volume is " << M3DG->Sim_handler->Trgt_vol << " and the final volume is " << M3DG->geometry->totalVolume() << "\n";
    std::cout << "The target area is " << M3DG->Sim_handler->Trgt_area << " and the final area is " << M3DG->geometry->totalArea() << "\n";
}

bool NormalNLP::intermediate_callback(
    AlgorithmMode mode,
    Index iter,
    Number obj_value,
    Number inf_pr,
    Number inf_du,
    Number mu,
    Number d_norm,
    Number regularization_size,
    Number alpha_du,
    Number alpha_pr,
    Index ls_trials,
    const IpoptData *ip_data,
    IpoptCalculatedQuantities *ip_cq)
{

    M3DG->Save_mesh(iter);

    // M3DG->geometry->totalArea();
    // M3DG->geometry->totalVolume();
    // std::cout << "The total area is " << M3DG->geometry->totalArea() << " and the total volume is " << M3DG->geometry->totalVolume() << "\n";

    return true;
}