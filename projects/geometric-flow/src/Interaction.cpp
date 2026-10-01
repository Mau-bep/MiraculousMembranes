

#include <Eigen/Core>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include <fstream>
#include <omp.h>
#include "Interaction.h"
#include "GeometryHelpers.h"
#include "BeadGeometry.h"
#include "Beads.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
typedef Eigen::Triplet<double> T;

#include <stdexcept>


void Interaction::validateSetup() const
{
    if (mesh == nullptr)
    {
        throw std::runtime_error(
            "Interaction: mesh has not been assigned.");
    }

    if (geometry == nullptr)
    {
        throw std::runtime_error(
            "Interaction: geometry has not been assigned.");
    }

    if (Bead_1 == nullptr)
    {
        throw std::runtime_error(
            "Interaction: Bead_1 has not been assigned.");
    }
}

double Interaction::Bond_energy()
{
    double E_bond = 0.0;
    Vector3 vec_r;
    double r;
    std::vector<double> params;

    // if (Bead_1->state == "manual")
    //     return E_bond;

    // std::cout<<"Calculating bond energy\n";

    for (size_t i = 0; i < Bead_1->Beads.size(); i++)
    {
        params = Bead_1->Interaction_constants_vector[i];
        if (Bead_1->Bond_type[i] == "Harmonic")
        {
            vec_r = Bead_1->Pos - Bead_1->Beads[i]->Pos;
            r = vec_r.norm();
            E_bond += 0.25 * params[0] * (r - params[1]) * (r - params[1]);
        }
        if (Bead_1->Bond_type[i] == "Lineal")
        {
            vec_r = Bead_1->Pos - Bead_1->Beads[i]->Pos;
            r = vec_r.norm();
            E_bond += 0.5 * params[0] * r;
        }
    }

    return E_bond;
}

std::vector<T> Interaction::Hessian_bonds_triplet()
{
    int B_1 = Bead_1->Bead_id;
    int B_2;
    int N_vert = mesh->nVertices();
    std::vector<T> tripletList;
    Eigen::Matrix<double, 3, 3> H_bond;

    // We decided that frozen beads will contribute to the Hessian
    // if (Bead_1->state == "manual")
    // return tripletList;

    for (int i = 0; i < Bead_1->Beads.size(); i++)
    {
        B_2 = Bead_1->Beads[i]->Bead_id;

        if (Bead_1->Bond_type[i] == "Harmonic")
        {
            //
            H_bond = 0.5 * Eigen::Matrix3d::Identity() * Bead_1->Interaction_constants_vector[i][0];
            for (int diag = 0; diag < 3; diag++)
            {
                // I will be adding this onto the corresponding beads
                tripletList.push_back(T(3 * (N_vert + B_1) + diag, 3 * (N_vert + B_1) + diag, H_bond(diag, diag)));
                tripletList.push_back(T(3 * (N_vert + B_2) + diag, 3 * (N_vert + B_2) + diag, H_bond(diag, diag)));
                tripletList.push_back(T(3 * (N_vert + B_1) + diag, 3 * (N_vert + B_2) + diag, -1 * H_bond(diag, diag)));
                tripletList.push_back(T(3 * (N_vert + B_2) + diag, 3 * (N_vert + B_1) + diag, -1 * H_bond(diag, diag)));
            }
        }

        if (Bead_1->Bond_type[i] == "Lineal")
            continue;
    }
    return tripletList;
}

VertexData<double> Normal_dot_Interaction::V_Tot_Energy()
{

    Vector3 vec_r;
    Vector3 Normal;
    double Q;
    int outside = 1;
    if (Energy_constants.size() > 3)
        outside = static_cast<int>(Energy_constants[3]);

    VertexData<double> VertexE(*mesh, 0.0);

    for (Face f : mesh->faces())
    {

        Normal = geometry->faceNormal_vec(f);
        vec_r = Bead_1->Pos - geometry->inputVertexPositions[f.halfedge().vertex()];
        Q = outside * dot(Normal, vec_r);
        if (Q < 0)
            continue;
        for (Vertex v : f.adjacentVertices())
        {

            vec_r = Bead_1->Pos - geometry->inputVertexPositions[v];
            Q = outside * dot(Normal, vec_r); // I will need to check this
            double r = vec_r.norm();
            VertexE[v] += (1.0 / 6.0) * E_r(r, Energy_constants) * Q / r;
        }
    }

    return VertexE;
    //
}

double Normal_dot_Interaction::V_Energy(Vertex v)
{
    Vector3 vec_r;
    Vector3 Normal;
    double Q;
    int outside = 1;
    if (Energy_constants.size() > 3)
        outside = static_cast<int>(Energy_constants[3]);
    double Energy = 0.0;
    for (Face f : v.adjacentFaces())
    {
        Normal = geometry->faceNormal_vec(f);
        vec_r = Bead_1->Pos - geometry->inputVertexPositions[f.halfedge().vertex()];
        Q = outside * dot(Normal, vec_r);
        if (Q < 0)
            continue;
        double r = vec_r.norm();
        Energy += (1.0 / 6.0) * E_r(r, Energy_constants) * Q / r;
    }
    return Energy;
}
double Normal_dot_Interaction::Tot_Energy()
{
    // std::cout<<"THe total energy is being calculated\n";
    double Total_E = 0.0;
    double rc = Energy_constants[2];
    // std::cout<<"THe cutoff is" <<rc <<" \n";
    Vector3 vec_r;
    Vector3 Normal;
    double Q;
    int outside = 1;
    if (Energy_constants.size() > 3)
        outside = static_cast<int>(Energy_constants[3]);

    for (Face f : mesh->faces())
    {

        Normal = geometry->faceNormal_vec(f);
        vec_r = Bead_1->Pos - geometry->inputVertexPositions[f.halfedge().vertex()];
        Q = outside * dot(Normal, vec_r);
        if (Q < 0)
            continue;
        for (Vertex v : f.adjacentVertices())
        {
            vec_r = Bead_1->Pos - geometry->inputVertexPositions[v];
            // Q should be recalculated here
            double r = vec_r.norm();
            Total_E += (1.0 / 6.0) * E_r(r, Energy_constants) * Q / r;
        }
    }

    Total_E += Bond_energy();
    return Total_E;
}

Vector3 Interaction::Bond_force()
{

    Vector3 r_vec;
    Vector3 Force{0.0, 0.0, 0.0};
    std::vector<double> params;
    double r;
    // if (Bead_1->state == "manual")
    // return Force;
    for (size_t i = 0; i < Bead_1->Beads.size(); i++)
    {
        params = Bead_1->Interaction_constants_vector[i];
        r_vec = Bead_1->Pos - Bead_1->Beads[i]->Pos;
        r = r_vec.norm();

        if (Bead_1->Bond_type[i] == "Harmonic")
        {
            Force += -1 * params[0] * (r - params[1]) * r_vec / r;
        }
        if (Bead_1->Bond_type[i] == "Lineal")
        {
            Force += -1 * params[0] * r_vec.unit();
        }
    }

    return Force;
}

VertexData<Vector3> Normal_dot_Interaction::Gradient()
{

    Vertex v;

    VertexData<Vector3> Force(*mesh, {0.0, 0.0, 0.0});
    double rc = Energy_constants[2];
    int outside = 1;
    if (Energy_constants.size() > 3)
        outside = static_cast<int>(Energy_constants[3]);

    // Ok so what is going on here

    Halfedge he;
    // I need to check that there is a nonzero total force
    Bead_1->Prev_Total_force = Bead_1->Total_force;
    Bead_1->Total_force = Vector3{0.0, 0.0, 0.0};
    Bead_1->Total_force += Bond_force();
    Bead_1->Total_force += Bead_1->CoverageForce;

    Eigen::Vector<double, 13> Positions_triangle;
    Eigen::Vector<double, 6> Positions_r;
    Eigen::Vector<double, 12> Grad_Q;
    Eigen::Vector<double, 6> Grad_r;
    std::array<int, 4> Vertices_triangle;

    Vector3 Force_vector;
    Vector3 Positions_r_vec;

    VertexData<Vector3> vector_bead(*mesh);
    VertexData<double> distance_bead(*mesh, 0.0);
    VertexData<double> Energy_contribution(*mesh, 0.0);
    VertexData<double> d_Energy_contribution(*mesh, 0.0);
    Vector3 Bead_pos = Bead_1->Pos;

    for (Vertex v : mesh->vertices())
    {
        vector_bead[v] = Bead_pos - geometry->inputVertexPositions[v];
        distance_bead[v] = vector_bead[v].norm();
        Energy_contribution[v] = E_r(distance_bead[v], Energy_constants);
        d_Energy_contribution[v] = dE_r(distance_bead[v], Energy_constants);
    }

    Vector3 Normal;
    double Q;
    double C1;
    double C2;

    for (Face f : mesh->faces())
    {

        he = f.halfedge();
        Vertices_triangle[0] = he.vertex().getIndex();
        Vertices_triangle[1] = he.next().vertex().getIndex();
        Vertices_triangle[2] = he.next().next().vertex().getIndex();
        Vertices_triangle[3] = mesh->nVertices();
        Positions_triangle << geometry->inputVertexPositions[Vertices_triangle[0]].x, geometry->inputVertexPositions[Vertices_triangle[0]].y, geometry->inputVertexPositions[Vertices_triangle[0]].z,
            geometry->inputVertexPositions[Vertices_triangle[1]].x, geometry->inputVertexPositions[Vertices_triangle[1]].y, geometry->inputVertexPositions[Vertices_triangle[1]].z,
            geometry->inputVertexPositions[Vertices_triangle[2]].x, geometry->inputVertexPositions[Vertices_triangle[2]].y, geometry->inputVertexPositions[Vertices_triangle[2]].z,
            Bead_pos.x, Bead_pos.y, Bead_pos.z, 1;
        Normal = geometry->faceNormal_vec(f);

        Q = outside * dot(vector_bead[Vertices_triangle[0]], Normal);

        if ((Q < 0))
        {
            continue;
        }

        Grad_Q = outside * geometry->gradient_triple_product(Positions_triangle);
        C2 = 0;

        C2 = C2 + Energy_contribution[Vertices_triangle[0]] / distance_bead[Vertices_triangle[0]];
        C2 = C2 + Energy_contribution[Vertices_triangle[1]] / distance_bead[Vertices_triangle[1]];
        C2 = C2 + Energy_contribution[Vertices_triangle[2]] / distance_bead[Vertices_triangle[2]];

        C2 = C2 / 6.0;

        for (size_t i = 0; i < 3; i++)
        {
            if (Q < 0 && (C2 > 1e-7 || C2 < -1e-7))
            {
                std::cout << " THe dot product is negative but C2 is not? \n";
                std::cout << "Q is" << Q << " and C2 is " << C2 << "\n";
            }
            Force_vector = Vector3{Grad_Q[3 * i], Grad_Q[3 * i + 1], Grad_Q[3 * i + 2]} * C2;
            Force[Vertices_triangle[i]] -= Force_vector;
        }
        Force_vector = Vector3{Grad_Q[9], Grad_Q[10], Grad_Q[11]} * C2;
        Bead_1->Total_force -= Force_vector; // Make sure the total force is 0 at the end of the interaction

        for (Vertex v : f.adjacentVertices())
        {

            // C2 = Energy_contribution[v]/distance_bead[v]/6.0;
            // if (distance_bead[v] > rc && rc > 0.0)
            // {
            //     // std::cout<<"Skipping vertex " << v.getIndex() << " because it is outside the cutoff distance\n";
            //     continue; // Skip this vertex if the bead is outside the cutoff distance or the normal is not pointing towards the bead
            // }
            Positions_r << geometry->inputVertexPositions[v].x, geometry->inputVertexPositions[v].y, geometry->inputVertexPositions[v].z,
                Bead_pos.x, Bead_pos.y, Bead_pos.z;

            Grad_r = geometry->gradient_r(Positions_r, distance_bead[v]);

            C1 = (d_Energy_contribution[v] / (distance_bead[v]) - Energy_contribution[v] / (distance_bead[v] * distance_bead[v])) * Q / 6.0;

            Force_vector = Vector3{Grad_r[0], Grad_r[1], Grad_r[2]} * C1;
            Force[v] -= Force_vector;
            Force_vector = Vector3{Grad_r[3], Grad_r[4], Grad_r[5]} * C1;
            Bead_1->Total_force -= Force_vector; // Make sure the total force is 0 at the end of the interaction
        }

        // std::cout<<" \n";
    }
    // std::cout<<"Th"
    return Force;
}

SparseMatrix<double> Normal_dot_Interaction::Hessian()
{
    // std::cout<<"Calculating Hessian for Normal dot interaction\n";
    // This matrix should be a (3N+3) by (3N+3) matrix but we really dont need all those number
    double rc = Energy_constants[2];
    int N_vert = mesh->nVertices();

    int outside = 1;
    if (Energy_constants.size() > 3)
        outside = static_cast<int>(Energy_constants[3]);

    // std::cout<<"Declaring the sparsematric\n";
    // std::cout<<"The number of total beads is " << Bead_1->Total_beads << "\n";
    SparseMatrix<double> Hessian(3 * N_vert + 3 * Bead_1->Total_beads, 3 * N_vert + 3 * Bead_1->Total_beads);

    // std::cout<<"The size of the Hessian is " << Hessian.rows() << " by " << Hessian.cols() << "\n";

    std::vector<T> tripletList = Hessian_bonds_triplet();

    Halfedge he;
    Eigen::Vector<double, 13> Positions_triangle;
    Eigen::Vector<double, 6> Positions_r;
    Eigen::Vector<double, 12> Grad_Q;
    Eigen::Vector<double, 6> Grad_r;
    std::array<int, 4> Vertices_triangle;
    std::array<int, 2> Vertices_r;

    Eigen::Matrix<double, 12, 12> Hessian_block_Q;
    Eigen::Matrix<double, 6, 6> Hessian_block_r;

    Eigen::Matrix<double, 6, 6> M_6_6;
    Eigen::Matrix<double, 6, 12> M_6_12;
    Eigen::Matrix<double, 12, 6> M_12_6;

    Vector3 Force_vector;
    Vector3 Positions_r_vec;
    VertexData<Vector3> vector_bead(*mesh);
    VertexData<double> distance_bead(*mesh, 0.0);
    VertexData<double> Energy_contribution(*mesh, 0.0);
    VertexData<double> d_Energy_contribution(*mesh, 0.0);
    VertexData<double> dd_Energy_contribution(*mesh, 0.0);
    Vector3 Bead_pos = Bead_1->Pos;

    for (Vertex v : mesh->vertices())
    {
        vector_bead[v] = Bead_pos - geometry->inputVertexPositions[v];
        distance_bead[v] = vector_bead[v].norm();
        Energy_contribution[v] = E_r(distance_bead[v], Energy_constants);
        d_Energy_contribution[v] = dE_r(distance_bead[v], Energy_constants);
        dd_Energy_contribution[v] = ddE_r(distance_bead[v], Energy_constants);
    }

    Vector3 Normal;
    double Q;
    double C;

    Vertices_triangle[3] = mesh->nVertices() + Bead_1->Bead_id; // This is the bead assigned index, we will add it at the end of the triangle
    Vertices_r[1] = mesh->nVertices() + Bead_1->Bead_id;        // This is the bead assigned index, we will add it at the end of the r vector
    // std::cout<<"THe bead indices has been used\n";
    for (Face f : mesh->faces())
    {
        // std::cout<<"\n";
        he = f.halfedge();
        Vertices_triangle[0] = he.vertex().getIndex();
        Vertices_triangle[1] = he.next().vertex().getIndex();
        Vertices_triangle[2] = he.next().next().vertex().getIndex();

        if (he.vertex().isBoundary() || he.next().vertex().isBoundary() || he.next().next().vertex().isBoundary())
            continue;

        Positions_triangle << geometry->inputVertexPositions[Vertices_triangle[0]].x, geometry->inputVertexPositions[Vertices_triangle[0]].y, geometry->inputVertexPositions[Vertices_triangle[0]].z,
            geometry->inputVertexPositions[Vertices_triangle[1]].x, geometry->inputVertexPositions[Vertices_triangle[1]].y, geometry->inputVertexPositions[Vertices_triangle[1]].z,
            geometry->inputVertexPositions[Vertices_triangle[2]].x, geometry->inputVertexPositions[Vertices_triangle[2]].y, geometry->inputVertexPositions[Vertices_triangle[2]].z,
            Bead_pos.x, Bead_pos.y, Bead_pos.z, 1;
        Normal = geometry->faceNormal_vec(f);

        Q = outside * dot(vector_bead[Vertices_triangle[0]], Normal);
        // double Q2 = dot( vector_bead[Vertices_triangle[1]], Normal );
        // double Q3 = dot( vector_bead[Vertices_triangle[2]], Normal );

        if (Q < 0)
        {
            continue;
        }
        Grad_Q = outside * geometry->gradient_triple_product(Positions_triangle);

        C = (Energy_contribution[Vertices_triangle[0]] / distance_bead[Vertices_triangle[0]] + Energy_contribution[Vertices_triangle[1]] / distance_bead[Vertices_triangle[1]] + Energy_contribution[Vertices_triangle[2]] / distance_bead[Vertices_triangle[2]]) / 6.0;

        Hessian_block_Q = outside * C * geometry->hessian_triple_product(Positions_triangle);
        // I have the first term

        for (size_t row = 0; row < 12; row++)
        {
            for (size_t col = 0; col < 12; col++)
            {
                tripletList.push_back(T(3 * Vertices_triangle[row / 3] + row % 3, 3 * Vertices_triangle[col / 3] + col % 3, Hessian_block_Q(row, col)));
            }
        }

        // That was the fastest term I still need to add a correction for the index of the Bead, this only works for one bead

        for (Vertex v : f.adjacentVertices())
        {
            // We need to do the whole thing now
            Vertices_r[0] = v.getIndex();

            // if (distance_bead[v] > rc && rc > 0.0)
            // {
            //     continue; // Skip this vertex if the bead is outside the cutoff distance or the normal is not pointing towards the bead
            // }
            Positions_r << geometry->inputVertexPositions[v].x, geometry->inputVertexPositions[v].y, geometry->inputVertexPositions[v].z,
                Bead_pos.x, Bead_pos.y, Bead_pos.z;
            Grad_r = geometry->gradient_r(Positions_r, distance_bead[v]);
            Hessian_block_r = geometry->hessian_r(Positions_r, distance_bead[v]);

            C = (dd_Energy_contribution[v] / (distance_bead[v]) - 2 * d_Energy_contribution[v] / (distance_bead[v] * distance_bead[v]) + 2 * Energy_contribution[v] / (distance_bead[v] * distance_bead[v] * distance_bead[v])) * Q / 6.0;

            M_6_6 = C * Grad_r * Grad_r.transpose();

            C = (d_Energy_contribution[v] / distance_bead[v] - Energy_contribution[v] / (distance_bead[v] * distance_bead[v])) * Q / 6.0;

            Hessian_block_r = C * Hessian_block_r;

            // M_6_6 = M_6_6 + Hessian_block_r;

            for (size_t row = 0; row < 6; row++)
            {
                for (size_t col = 0; col < 6; col++)
                {
                    tripletList.push_back(T(3 * Vertices_r[row / 3] + row % 3, 3 * Vertices_r[col / 3] + col % 3, Hessian_block_r(row, col) + M_6_6(row, col)));
                    // tripletList.push_back(T(3*Vertices_r[row/3]+row%3, 3*Vertices_r[col/3]+col%3, ));
                }
            }
            // I just added the 2 terms that are 6x6, which are the hessian of r and the outer product of the gradients

            C = (d_Energy_contribution[v] / (distance_bead[v]) - Energy_contribution[v] / (distance_bead[v] * distance_bead[v])) / 6.0;
            M_6_12 = C * Grad_r * Grad_Q.transpose();
            M_12_6 = C * Grad_Q * Grad_r.transpose();

            for (size_t row = 0; row < 6; row++)
            {
                for (size_t col = 0; col < 12; col++)
                {
                    tripletList.push_back(T(3 * Vertices_r[row / 3] + row % 3, 3 * Vertices_triangle[col / 3] + col % 3, M_6_12(row, col)));
                }
            }

            for (size_t row = 0; row < 12; row++)
            {
                for (size_t col = 0; col < 6; col++)
                {
                    tripletList.push_back(T(3 * Vertices_triangle[row / 3] + row % 3, 3 * Vertices_r[col / 3] + col % 3, M_12_6(row, col)));
                }
            }
            // I just added the 2 terms that are 6x12 and 12x6, which have the product of the gradients of r and Q
        }
    }

    // We have all the nonzeroterms, now i need to add the zeros.

    // std::cout<<"Setting from triplets\n";
    Hessian.setFromTriplets(tripletList.begin(), tripletList.end());
    // std::cout<<"Hessian calculated\n";
    return Hessian;
}

SparseMatrix<double> Normal_dot_Interaction::Hessian_IP()
{
    // std::cout << "Calculating Hessian for Normal dot interaction\n";
    // This matrix should be a (3N+3) by (3N+3) matrix but we really dont need all those number
    double rc = Energy_constants[2];
    int N_vert = mesh->nVertices();

    int outside = 1;
    // std::cout<<"Energy constants size is " << Energy_constants.size() << "\n";
    if (Energy_constants.size() > 3)
    {
        outside = static_cast<int>(Energy_constants[3]);
        // std::cout << "Outside is " << outside << "\n";
    }

    // std::cout<<"Declaring the sparsematric\n";
    // std::cout<<"The number of total beads is " << Bead_1->Total_beads << "\n";
    SparseMatrix<double> Hessian(3 * N_vert + 3 * Bead_1->Total_beads, 3 * N_vert + 3 * Bead_1->Total_beads);

    // std::cout<<"The size of the Hessian is " << Hessian.rows() << " by " << Hessian.cols() << "\n";

    std::vector<T> tripletList = Hessian_bonds_triplet();

    Halfedge he;
    Eigen::Vector<double, 13> Positions_triangle;
    Eigen::Vector<double, 6> Positions_r;
    Eigen::Vector<double, 12> Grad_Q;
    Eigen::Vector<double, 6> Grad_r;
    std::array<int, 4> Vertices_triangle;
    std::array<int, 2> Vertices_r;

    Eigen::Matrix<double, 12, 12> Hessian_block_Q;
    Eigen::Matrix<double, 6, 6> Hessian_block_r;

    Eigen::Matrix<double, 6, 6> M_6_6;
    Eigen::Matrix<double, 6, 12> M_6_12;
    Eigen::Matrix<double, 12, 6> M_12_6;

    Vector3 Force_vector;
    Vector3 Positions_r_vec;
    VertexData<Vector3> vector_bead(*mesh);
    VertexData<double> distance_bead(*mesh, 0.0);
    VertexData<double> Energy_contribution(*mesh, 0.0);
    VertexData<double> d_Energy_contribution(*mesh, 0.0);
    VertexData<double> dd_Energy_contribution(*mesh, 0.0);
    Vector3 Bead_pos = Bead_1->Pos;

    for (Vertex v : mesh->vertices())
    {
        vector_bead[v] = Bead_pos - geometry->inputVertexPositions[v];
        distance_bead[v] = vector_bead[v].norm();
        Energy_contribution[v] = E_r(distance_bead[v], Energy_constants);
        d_Energy_contribution[v] = dE_r(distance_bead[v], Energy_constants);
        dd_Energy_contribution[v] = ddE_r(distance_bead[v], Energy_constants);
    }

    Vector3 Normal;
    double Q;
    double C;

    Vertices_triangle[3] = mesh->nVertices() + Bead_1->Bead_id; // This is the bead assigned index, we will add it at the end of the triangle
    Vertices_r[1] = mesh->nVertices() + Bead_1->Bead_id;        // This is the bead assigned index, we will add it at the end of the r vector
    // std::cout<<"THe bead indices has been used\n";
    for (Face f : mesh->faces())
    {
        // std::cout<<"\n";
        he = f.halfedge();
        Vertices_triangle[0] = he.vertex().getIndex();
        Vertices_triangle[1] = he.next().vertex().getIndex();
        Vertices_triangle[2] = he.next().next().vertex().getIndex();

        if (he.vertex().isBoundary() || he.next().vertex().isBoundary() || he.next().next().vertex().isBoundary())
            continue;

        Positions_triangle << geometry->inputVertexPositions[Vertices_triangle[0]].x, geometry->inputVertexPositions[Vertices_triangle[0]].y, geometry->inputVertexPositions[Vertices_triangle[0]].z,
            geometry->inputVertexPositions[Vertices_triangle[1]].x, geometry->inputVertexPositions[Vertices_triangle[1]].y, geometry->inputVertexPositions[Vertices_triangle[1]].z,
            geometry->inputVertexPositions[Vertices_triangle[2]].x, geometry->inputVertexPositions[Vertices_triangle[2]].y, geometry->inputVertexPositions[Vertices_triangle[2]].z,
            Bead_pos.x, Bead_pos.y, Bead_pos.z, 1;
        Normal = geometry->faceNormal_vec(f);

        Q = outside * dot(vector_bead[Vertices_triangle[0]], Normal);
        // double Q2 = dot( vector_bead[Vertices_triangle[1]], Normal );
        // double Q3 = dot( vector_bead[Vertices_triangle[2]], Normal );

        if (Q < 0)
        {
            continue;
        }
        Grad_Q = outside * geometry->gradient_triple_product(Positions_triangle);

        C = (Energy_contribution[Vertices_triangle[0]] / distance_bead[Vertices_triangle[0]] + Energy_contribution[Vertices_triangle[1]] / distance_bead[Vertices_triangle[1]] + Energy_contribution[Vertices_triangle[2]] / distance_bead[Vertices_triangle[2]]) / 6.0;

        Hessian_block_Q = outside * C * geometry->hessian_triple_product(Positions_triangle);
        // I have the first term

        for (size_t row = 0; row < 12; row++)
        {
            for (size_t col = 0; col < 12; col++)
            {
                tripletList.push_back(T(3 * Vertices_triangle[row / 3] + row % 3, 3 * Vertices_triangle[col / 3] + col % 3, Hessian_block_Q(row, col)));
            }
        }

        // That was the fastest term I still need to add a correction for the index of the Bead, this only works for one bead

        for (Vertex v : f.adjacentVertices())
        {
            // We need to do the whole thing now
            Vertices_r[0] = v.getIndex();

            if (distance_bead[v] > rc && rc > 0.0)
            {
                continue; // Skip this vertex if the bead is outside the cutoff distance or the normal is not pointing towards the bead
            }
            Positions_r << geometry->inputVertexPositions[v].x, geometry->inputVertexPositions[v].y, geometry->inputVertexPositions[v].z,
                Bead_pos.x, Bead_pos.y, Bead_pos.z;
            Grad_r = geometry->gradient_r(Positions_r, distance_bead[v]);
            Hessian_block_r = geometry->hessian_r(Positions_r, distance_bead[v]);

            C = (dd_Energy_contribution[v] / (distance_bead[v]) - 2 * d_Energy_contribution[v] / (distance_bead[v] * distance_bead[v]) + 2 * Energy_contribution[v] / (distance_bead[v] * distance_bead[v] * distance_bead[v])) * Q / 6.0;

            M_6_6 = C * Grad_r * Grad_r.transpose();

            C = (d_Energy_contribution[v] / distance_bead[v] - Energy_contribution[v] / (distance_bead[v] * distance_bead[v])) * Q / 6.0;

            Hessian_block_r = C * Hessian_block_r;

            // M_6_6 = M_6_6 + Hessian_block_r;

            for (size_t row = 0; row < 6; row++)
            {
                for (size_t col = 0; col < 6; col++)
                {
                    tripletList.push_back(T(3 * Vertices_r[row / 3] + row % 3, 3 * Vertices_r[col / 3] + col % 3, Hessian_block_r(row, col) + M_6_6(row, col)));
                    // tripletList.push_back(T(3*Vertices_r[row/3]+row%3, 3*Vertices_r[col/3]+col%3, ));
                }
            }
            // I just added the 2 terms that are 6x6, which are the hessian of r and the outer product of the gradients

            C = (d_Energy_contribution[v] / (distance_bead[v]) - Energy_contribution[v] / (distance_bead[v] * distance_bead[v])) / 6.0;
            M_6_12 = C * Grad_r * Grad_Q.transpose();
            M_12_6 = C * Grad_Q * Grad_r.transpose();

            for (size_t row = 0; row < 6; row++)
            {
                for (size_t col = 0; col < 12; col++)
                {
                    tripletList.push_back(T(3 * Vertices_r[row / 3] + row % 3, 3 * Vertices_triangle[col / 3] + col % 3, M_6_12(row, col)));
                }
            }

            for (size_t row = 0; row < 12; row++)
            {
                for (size_t col = 0; col < 6; col++)
                {
                    tripletList.push_back(T(3 * Vertices_triangle[row / 3] + row % 3, 3 * Vertices_r[col / 3] + col % 3, M_12_6(row, col)));
                }
            }
            // I just added the 2 terms that are 6x12 and 12x6, which have the product of the gradients of r and Q
        }
    }

    // We have all the nonzeroterms, now i need to add the zeros.

    for (Vertex v : mesh->vertices())
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {

                tripletList.push_back(T(3 * v.getIndex() + j, 3 * (mesh->nVertices() + Bead_1->Bead_id) + i, 0.0));
                tripletList.push_back(T(3 * (mesh->nVertices() + Bead_1->Bead_id) + i, 3 * v.getIndex() + j, 0.0));
            }
        }
    }

    for (auto &t : tripletList)
    {
        assert(t.row() >= 0 && t.row() < 3 * (N_vert + Bead_1->Total_beads) && "Row index OOB Bead");
        assert(t.col() >= 0 && t.col() < 3 * (N_vert + Bead_1->Total_beads) && "Col index OOB Bead");
    }
    // std::cout<<"Setting from triplets\n";
    Hessian.setFromTriplets(tripletList.begin(), tripletList.end());
    // std::cout<<"Hessian calculated\n";
    return Hessian;
}

double Integrated_Interaction::Tot_Energy()
{
    double Total_E = 0.0;
    double rc = Energy_constants[2];
    Vector3 vec_r;
    double face_area;
    for (Face f : mesh->faces())
    {
        face_area = geometry->faceArea(f);
        for (Vertex v : f.adjacentVertices())
        {
            vec_r = Bead_1->Pos - geometry->inputVertexPositions[v];
            double r = vec_r.norm();
            if (r < rc || rc < 0.0)
            {
                // if(geometry->inputVertexPositions[v].x>5.0){
                // std::cout<<"THis vertex is far and contributes with "<< E_r(r, Energy_constants)*face_area/3 <<" ? \n";
                // std::cout<<"The position of the vertex is "<< geometry->inputVertexPositions[v] << " and the bead is at " << Bead_1->Pos << "\n";
                // std::cout<<"tHE distance is "<< r << "\n";
                // std::cout<<"The cutoff distance is" << rc << "\n";
                // }

                Total_E += E_r(r, Energy_constants) * face_area / 3;
            }
        }
    }

    // So I need to add the ENergy of the bonds
    Total_E += Bond_energy();

    return Total_E;
}
VertexData<double> Integrated_Interaction::V_Tot_Energy()
{
    VertexData<double> vertexE(*mesh, 0.0);
    return vertexE;
}
double Integrated_Interaction::V_Energy(Vertex v)
{
    return 0.0;
}

VertexData<Vector3> Integrated_Interaction::Gradient()
{
    VertexData<Vector3> Force(*mesh, {0.0, 0.0, 0.0});
    Bead_1->Prev_Total_force = Bead_1->Total_force;
    Bead_1->Total_force = {0.0, 0.0, 0.0};
    Bead_1->Total_force += Bond_force();
    Bead_1->Total_force += Bead_1->CoverageForce;
    // std::cout<<"Force initialied\n";
    // std::cout<<"THis is bead " << Bead_1->Bead_id << "\n";
    double rc = Energy_constants[2];
    Vector3 vec_r;
    double face_area;
    Halfedge he;
    Eigen::Vector<double, 6> Positions_r;
    Eigen::Vector<double, 9> Positions_f;
    Eigen::Vector<double, 6> Grad_r;
    Eigen::Vector<double, 9> Grad_f;

    std::array<int, 3> Vertices_triangle;

    VertexData<Vector3> vector_bead(*mesh);
    VertexData<double> distance_bead(*mesh, 0.0);
    VertexData<double> Energy_contribution(*mesh, 0.0);
    VertexData<double> d_Energy_contribution(*mesh, 0.0);
    Vector3 Bead_pos = Bead_1->Pos;
    Vector3 Force_vector;
    for (Vertex v : mesh->vertices())
    {
        vector_bead[v] = Bead_pos - geometry->inputVertexPositions[v];
        distance_bead[v] = vector_bead[v].norm();
        Energy_contribution[v] = E_r(distance_bead[v], Energy_constants);
        d_Energy_contribution[v] = dE_r(distance_bead[v], Energy_constants);
        // std::cout<<"The distance for vertex "<< v.getIndex() << " is " << distance_bead[v] << " and the energy contribution is " << Energy_contribution[v] << "\n";
    }

    for (Face f : mesh->faces())
    {
        face_area = geometry->faceArea(f);
        he = f.halfedge();
        Vertices_triangle[0] = he.vertex().getIndex();
        Vertices_triangle[1] = he.next().vertex().getIndex();
        Vertices_triangle[2] = he.next().next().vertex().getIndex();

        // Ineed to check here something

        Positions_f = stack_positions<3>(*geometry, Vertices_triangle);
        Grad_f = geometry->gradient_triangle_area(Positions_f);
        //
        for (int i = 0; i < 3; i++)
        {
            Force_vector = {Grad_f[3 * i], Grad_f[3 * i + 1], Grad_f[3 * i + 2]};
            Force_vector *= (Energy_contribution[Vertices_triangle[0]] + Energy_contribution[Vertices_triangle[1]] + Energy_contribution[Vertices_triangle[2]]) / 3.0;
            Force[Vertices_triangle[i]] -= Force_vector;
        }

        for (Vertex v : f.adjacentVertices())
        {

            double r = distance_bead[v];
            if (r < rc || rc < 0.0)
            {
                Positions_r << geometry->inputVertexPositions[v].x, geometry->inputVertexPositions[v].y, geometry->inputVertexPositions[v].z,
                    Bead_pos.x, Bead_pos.y, Bead_pos.z;
                Grad_r = geometry->gradient_r(Positions_r, distance_bead[v]);
                Force_vector = Vector3{Grad_r[0], Grad_r[1], Grad_r[2]};
                Force_vector *= (1.0 / 3.0) * dE_r(r, Energy_constants) * face_area;
                Force[v] -= Force_vector;
                Bead_1->Total_force += Force_vector;
            }
        }
    }

    // std::cout<<"The bead total force should be" << Bead_1->Total_force << "\n";
    return Force;
}

SparseMatrix<double> Integrated_Interaction::Hessian()
{

    std::cout << "Getting rc\n";
    double rc = Energy_constants[2];
    std::cout << "Building hessian\n";
    SparseMatrix<double> Hessian(3 * mesh->nVertices() + 3 * Bead_1->Total_beads, 3 * mesh->nVertices() + 3 * Bead_1->Total_beads);
    // typedef Eigen::Triplet<double> T;

    std::vector<T> tripletList = Hessian_bonds_triplet();
    Halfedge he;
    Eigen::Vector<double, 6> Positions_r;
    Eigen::Vector<double, 9> Positions_f;
    Eigen::Vector<double, 6> Grad_r;
    Eigen::Vector<double, 9> Grad_f;
    std::array<int, 4> Vertices_triangle;
    std::array<int, 2> Vertices_r;
    Eigen::Matrix<double, 9, 9> Hessian_block_f;
    Eigen::Matrix<double, 6, 6> Hessian_block_r;
    Eigen::Matrix<double, 6, 9> M_6_9;
    Eigen::Matrix<double, 9, 6> M_9_6;
    Eigen::Matrix<double, 6, 6> M_6_6;

    Vector3 Force_vector;
    VertexData<Vector3> vector_bead(*mesh);
    VertexData<double> distance_bead(*mesh, 0.0);
    VertexData<double> Energy_contribution(*mesh, 0.0);
    VertexData<double> d_Energy_contribution(*mesh, 0.0);
    VertexData<double> dd_Energy_contribution(*mesh, 0.0);
    Vector3 Bead_pos = Bead_1->Pos;

    double C;
    double face_area;

    for (Vertex v : mesh->vertices())
    {
        vector_bead[v] = Bead_pos - geometry->inputVertexPositions[v];
        distance_bead[v] = vector_bead[v].norm();
        Energy_contribution[v] = E_r(distance_bead[v], Energy_constants);
        d_Energy_contribution[v] = dE_r(distance_bead[v], Energy_constants);
        dd_Energy_contribution[v] = ddE_r(distance_bead[v], Energy_constants);
    }

    Vertices_triangle[3] = mesh->nVertices() + Bead_1->Bead_id;
    Vertices_r[1] = mesh->nVertices() + Bead_1->Bead_id; // This is the bead assigned index, we will add it at the end of the r vector

    for (Face f : mesh->faces())
    {

        face_area = geometry->faceArea(f);
        he = f.halfedge();
        Vertices_triangle[0] = he.vertex().getIndex();
        Vertices_triangle[1] = he.next().vertex().getIndex();
        Vertices_triangle[2] = he.next().next().vertex().getIndex();

        Positions_f = stack_positions<3>(*geometry, Vertices_triangle);
        Grad_f = geometry->gradient_triangle_area(Positions_f);

        Hessian_block_f = geometry->hessian_triangle_area(Positions_f);

        C = (Energy_contribution[Vertices_triangle[0]] + Energy_contribution[Vertices_triangle[1]] + Energy_contribution[Vertices_triangle[2]]) / 3.0;
        Hessian_block_f = C * Hessian_block_f;

        for (size_t row = 0; row < 9; row++)
        {
            for (size_t col = 0; col < 9; col++)
            {
                tripletList.push_back(T(3 * Vertices_triangle[row / 3] + row % 3, 3 * Vertices_triangle[col / 3] + col % 3, Hessian_block_f(row, col)));
            }
        }

        // Now we iterate over the triangles in the face

        for (int i = 0; i < 3; i++)
        {
            Vertices_r[0] = Vertices_triangle[i];
            if (distance_bead[Vertices_triangle[i]] > rc && rc > 0.0)
            {
                continue; // Skip this vertex if the bead is outside the cutoff distance
            }
            Positions_r << geometry->inputVertexPositions[Vertices_triangle[i]].x, geometry->inputVertexPositions[Vertices_triangle[i]].y, geometry->inputVertexPositions[Vertices_triangle[i]].z,
                Bead_pos.x, Bead_pos.y, Bead_pos.z;
            Grad_r = geometry->gradient_r(Positions_r, distance_bead[Vertices_triangle[i]]);
            Hessian_block_r = geometry->hessian_r(Positions_r, distance_bead[Vertices_triangle[i]]);

            // I have all the elements now its time to do the terms

            // We have  Ei' Grad_A Grad_r^T  + Ei' Grad_r Grad_A^T
            C = (d_Energy_contribution[Vertices_triangle[i]]) / 3.0;

            M_9_6 = C * Grad_f * Grad_r.transpose();
            M_6_9 = C * Grad_r * Grad_f.transpose();

            for (size_t row = 0; row < 9; row++)
            {
                for (size_t col = 0; col < 6; col++)
                {
                    tripletList.push_back(T(3 * Vertices_triangle[row / 3] + row % 3, 3 * Vertices_r[col / 3] + col % 3, M_9_6(row, col)));
                }
            }
            for (size_t row = 0; row < 6; row++)
            {
                for (size_t col = 0; col < 9; col++)
                {
                    tripletList.push_back(T(3 * Vertices_r[row / 3] + row % 3, 3 * Vertices_triangle[col / 3] + col % 3, M_6_9(row, col)));
                }
            }

            // Now we add A E'' Grad_r Grad_r^T

            C = face_area * dd_Energy_contribution[Vertices_triangle[i]] / 3.0;
            M_6_6 = C * Grad_r * Grad_r.transpose();
            for (size_t row = 0; row < 6; row++)
            {
                for (size_t col = 0; col < 6; col++)
                {
                    tripletList.push_back(T(3 * Vertices_r[row / 3] + row % 3, 3 * Vertices_r[col / 3] + col % 3, M_6_6(row, col)));
                }
            }

            // Now we add the term A Ei' Hess(r)

            C = face_area * d_Energy_contribution[Vertices_triangle[i]] / 3.0;
            Hessian_block_r = C * Hessian_block_r;
            for (size_t row = 0; row < 6; row++)
            {
                for (size_t col = 0; col < 6; col++)
                {
                    tripletList.push_back(T(3 * Vertices_r[row / 3] + row % 3, 3 * Vertices_r[col / 3] + col % 3, Hessian_block_r(row, col)));
                }
            }
        }
    }

    // Now we have all the terms, we just need to set the triplets
    Hessian.setFromTriplets(tripletList.begin(), tripletList.end());

    return Hessian;
}

SparseMatrix<double> Integrated_Interaction::Hessian_IP()
{

    std::cout << "Getting rc\n";
    double rc = Energy_constants[2];
    std::cout << "Building hessian\n";
    SparseMatrix<double> Hessian(3 * mesh->nVertices() + 3 * Bead_1->Total_beads, 3 * mesh->nVertices() + 3 * Bead_1->Total_beads);
    // typedef Eigen::Triplet<double> T;

    std::vector<T> tripletList = Hessian_bonds_triplet();
    Halfedge he;
    Eigen::Vector<double, 6> Positions_r;
    Eigen::Vector<double, 9> Positions_f;
    Eigen::Vector<double, 6> Grad_r;
    Eigen::Vector<double, 9> Grad_f;
    std::array<int, 4> Vertices_triangle;
    std::array<int, 2> Vertices_r;
    Eigen::Matrix<double, 9, 9> Hessian_block_f;
    Eigen::Matrix<double, 6, 6> Hessian_block_r;
    Eigen::Matrix<double, 6, 9> M_6_9;
    Eigen::Matrix<double, 9, 6> M_9_6;
    Eigen::Matrix<double, 6, 6> M_6_6;

    Vector3 Force_vector;
    VertexData<Vector3> vector_bead(*mesh);
    VertexData<double> distance_bead(*mesh, 0.0);
    VertexData<double> Energy_contribution(*mesh, 0.0);
    VertexData<double> d_Energy_contribution(*mesh, 0.0);
    VertexData<double> dd_Energy_contribution(*mesh, 0.0);
    Vector3 Bead_pos = Bead_1->Pos;

    double C;
    double face_area;

    for (Vertex v : mesh->vertices())
    {
        vector_bead[v] = Bead_pos - geometry->inputVertexPositions[v];
        distance_bead[v] = vector_bead[v].norm();
        Energy_contribution[v] = E_r(distance_bead[v], Energy_constants);
        d_Energy_contribution[v] = dE_r(distance_bead[v], Energy_constants);
        dd_Energy_contribution[v] = ddE_r(distance_bead[v], Energy_constants);
    }

    Vertices_triangle[3] = mesh->nVertices() + Bead_1->Bead_id;
    Vertices_r[1] = mesh->nVertices() + Bead_1->Bead_id; // This is the bead assigned index, we will add it at the end of the r vector

    for (Face f : mesh->faces())
    {

        face_area = geometry->faceArea(f);
        he = f.halfedge();
        Vertices_triangle[0] = he.vertex().getIndex();
        Vertices_triangle[1] = he.next().vertex().getIndex();
        Vertices_triangle[2] = he.next().next().vertex().getIndex();

        Positions_f = stack_positions<3>(*geometry, Vertices_triangle);
        Grad_f = geometry->gradient_triangle_area(Positions_f);

        Hessian_block_f = geometry->hessian_triangle_area(Positions_f);

        C = (Energy_contribution[Vertices_triangle[0]] + Energy_contribution[Vertices_triangle[1]] + Energy_contribution[Vertices_triangle[2]]) / 3.0;
        Hessian_block_f = C * Hessian_block_f;

        for (size_t row = 0; row < 9; row++)
        {
            for (size_t col = 0; col < 9; col++)
            {
                tripletList.push_back(T(3 * Vertices_triangle[row / 3] + row % 3, 3 * Vertices_triangle[col / 3] + col % 3, Hessian_block_f(row, col)));
            }
        }

        // Now we iterate over the triangles in the face

        for (int i = 0; i < 3; i++)
        {
            Vertices_r[0] = Vertices_triangle[i];
            if (distance_bead[Vertices_triangle[i]] > rc && rc > 0.0)
            {
                continue; // Skip this vertex if the bead is outside the cutoff distance
            }
            Positions_r << geometry->inputVertexPositions[Vertices_triangle[i]].x, geometry->inputVertexPositions[Vertices_triangle[i]].y, geometry->inputVertexPositions[Vertices_triangle[i]].z,
                Bead_pos.x, Bead_pos.y, Bead_pos.z;
            Grad_r = geometry->gradient_r(Positions_r, distance_bead[Vertices_triangle[i]]);
            Hessian_block_r = geometry->hessian_r(Positions_r, distance_bead[Vertices_triangle[i]]);

            // I have all the elements now its time to do the terms

            // We have  Ei' Grad_A Grad_r^T  + Ei' Grad_r Grad_A^T
            C = (d_Energy_contribution[Vertices_triangle[i]]) / 3.0;

            M_9_6 = C * Grad_f * Grad_r.transpose();
            M_6_9 = C * Grad_r * Grad_f.transpose();

            for (size_t row = 0; row < 9; row++)
            {
                for (size_t col = 0; col < 6; col++)
                {
                    tripletList.push_back(T(3 * Vertices_triangle[row / 3] + row % 3, 3 * Vertices_r[col / 3] + col % 3, M_9_6(row, col)));
                }
            }
            for (size_t row = 0; row < 6; row++)
            {
                for (size_t col = 0; col < 9; col++)
                {
                    tripletList.push_back(T(3 * Vertices_r[row / 3] + row % 3, 3 * Vertices_triangle[col / 3] + col % 3, M_6_9(row, col)));
                }
            }

            // Now we add A E'' Grad_r Grad_r^T

            C = face_area * dd_Energy_contribution[Vertices_triangle[i]] / 3.0;
            M_6_6 = C * Grad_r * Grad_r.transpose();
            for (size_t row = 0; row < 6; row++)
            {
                for (size_t col = 0; col < 6; col++)
                {
                    tripletList.push_back(T(3 * Vertices_r[row / 3] + row % 3, 3 * Vertices_r[col / 3] + col % 3, M_6_6(row, col)));
                }
            }

            // Now we add the term A Ei' Hess(r)

            C = face_area * d_Energy_contribution[Vertices_triangle[i]] / 3.0;
            Hessian_block_r = C * Hessian_block_r;
            for (size_t row = 0; row < 6; row++)
            {
                for (size_t col = 0; col < 6; col++)
                {
                    tripletList.push_back(T(3 * Vertices_r[row / 3] + row % 3, 3 * Vertices_r[col / 3] + col % 3, Hessian_block_r(row, col)));
                }
            }
        }
    }

    for (Vertex v : mesh->vertices())
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {

                tripletList.push_back(T(3 * v.getIndex() + j, 3 * (mesh->nVertices() + Bead_1->Bead_id) + i, 0.0));
                tripletList.push_back(T(3 * (mesh->nVertices() + Bead_1->Bead_id) + i, 3 * v.getIndex() + j, 0.0));
            }
        }
    }

    // Now we have all the terms, we just need to set the triplets
    Hessian.setFromTriplets(tripletList.begin(), tripletList.end());

    return Hessian;
}

VertexData<Vector3> No_mem_Inter::Gradient()
{
    VertexData<Vector3> Force(*mesh, {0.0, 0.0, 0.0});
    Bead_1->Prev_Total_force = Bead_1->Total_force;
    Bead_1->Total_force = Bond_force();
    Bead_1->Total_force += Bead_1->CoverageForce;
    return Force;
}

SparseMatrix<double> No_mem_Inter::Hessian()
{
    SparseMatrix<double> Hessian(3 * mesh->nVertices() + 3 * Bead_1->Total_beads, 3 * mesh->nVertices() + 3 * Bead_1->Total_beads);
    typedef Eigen::Triplet<double> T;

    std::vector<T> tripletList = Hessian_bonds_triplet();
    Hessian.setFromTriplets(tripletList.begin(), tripletList.end());
    return Hessian;
}

SparseMatrix<double> No_mem_Inter::Hessian_IP()
{
    return Hessian();
}

double Face_Integrated_Interaction::Tot_Energy()
{
    validateSetup();

    double energy = 0.0;

    for (Face f : mesh->faces())
    {
        energy += Face_Energy(f);
    }

    return energy + Bond_energy();
}

VertexData<double> Face_Integrated_Interaction::V_Tot_Energy()
{
    validateSetup();

    VertexData<double> vertexEnergy(*mesh, 0.0);

    for (Face f : mesh->faces())
    {
        const double faceEnergy = Face_Energy(f);

        Halfedge he = f.halfedge();

        Vertex v0 = he.vertex();
        Vertex v1 = he.next().vertex();
        Vertex v2 = he.next().next().vertex();

        /*
         * This is diagnostic bookkeeping only. The face energy has no
         * unique vertex decomposition, so distribute it equally.
         */
        vertexEnergy[v0] += faceEnergy / 3.0;
        vertexEnergy[v1] += faceEnergy / 3.0;
        vertexEnergy[v2] += faceEnergy / 3.0;
    }

    return vertexEnergy;
}

double Face_Integrated_Interaction::V_Energy(Vertex v)
{
    validateSetup();

    double energy = 0.0;

    for (Face f : v.adjacentFaces())
    {
        energy += Face_Energy(f) / 3.0;
    }

    return energy;
}

VertexData<Vector3> Face_Integrated_Interaction::Gradient()
{
    validateSetup();

    VertexData<Vector3> membraneForce(
        *mesh,
        Vector3({0.0, 0.0, 0.0}));

    // Same bead-force convention as the other interactions: this Gradient()
    // starts the bead's force for this evaluation, then adds bonds, the
    // Coverage energy's contribution and the face terms
    Bead_1->Prev_Total_force = Bead_1->Total_force;
    Bead_1->Total_force = Vector3{0.0, 0.0, 0.0};
    Bead_1->Total_force += Bond_force();
    Bead_1->Total_force += Bead_1->CoverageForce;

    Vector3 beadForce({0.0, 0.0, 0.0});
    for (Face f : mesh->faces())
    {
        Add_Face_Force(f, membraneForce, beadForce);
    }
    Bead_1->Total_force += beadForce;

    return membraneForce;
}

SparseMatrix<double> Face_Integrated_Interaction::Hessian()
{
    validateSetup();

    // No analytic Hessian for the face terms yet: only the bonds, in the
    // same layout as the other interactions (3 dofs per vertex, then per bead)
    const size_t nDof = 3 * mesh->nVertices() + 3 * Bead_1->Total_beads;
    SparseMatrix<double> hessian(nDof, nDof);
    std::vector<Eigen::Triplet<double>> tripletList = Hessian_bonds_triplet();
    hessian.setFromTriplets(tripletList.begin(), tripletList.end());
    return hessian;
}

SparseMatrix<double> Face_Integrated_Interaction::Hessian_IP()
{
    return Hessian();
}

double Adhesion::Face_Energy(Face f)
{
    const bead_geometry::FaceCoverageData face =
        bead_geometry::evaluateFaceCoverage(*geometry, f, Bead_1->Pos, sigma());
    if (!face.selected)
        return 0.0;
    return -strength() * sigma() * sigma() * face.weight * face.omegaUnsigned;
}

void Adhesion::Add_Face_Force(Face f, VertexData<Vector3> &membraneForce, Vector3 &beadForce)
{
    const bead_geometry::FaceCoverageData face =
        bead_geometry::evaluateFaceCoverage(*geometry, f, Bead_1->Pos, sigma());
    if (!face.selected)
        return;

    // a_f = sigma^2 w(r_f) |omega_f|, r_f measured to the centroid, so
    // dr_f/dp_i = radialDirection / 3 for each corner p_i
    const double sigma2 = sigma() * sigma();
    const Vector3 gradWeight = face.dWeightDr * face.radialDirection / 3.0;
    const Vector3 dA0 = sigma2 * (face.weight * face.gradOmega0 + face.omegaUnsigned * gradWeight);
    const Vector3 dA1 = sigma2 * (face.weight * face.gradOmega1 + face.omegaUnsigned * gradWeight);
    const Vector3 dA2 = sigma2 * (face.weight * face.gradOmega2 + face.omegaUnsigned * gradWeight);

    // E = -W a_f, force = -dE/dp = W da_f/dp
    const Vector3 force0 = strength() * dA0;
    const Vector3 force1 = strength() * dA1;
    const Vector3 force2 = strength() * dA2;

    Halfedge he = f.halfedge();
    membraneForce[he.vertex()] += force0;
    membraneForce[he.next().vertex()] += force1;
    membraneForce[he.next().next().vertex()] += force2;

    // a_f depends only on p_i - bead position, so da_f/dq = -sum_i da_f/dp_i
    beadForce -= force0 + force1 + force2;
}

Plane_Adhesion::Plane_Adhesion(ManifoldSurfaceMesh *inputMesh,
                               VertexPositionGeometry *inputGeometry,
                               std::vector<double> constants)
{
    mesh = inputMesh;
    geometry = inputGeometry;
    Energy_constants = constants;
    if (Energy_constants.size() < 6)
        throw std::invalid_argument("Plane_Adhesion: constants are {W, sigma, delta, nx, ny, nz}");
    if (!(width() > 0.0))
        throw std::invalid_argument("Plane_Adhesion: the contact width delta (rc) must be positive");
    Vector3 n{Energy_constants[3], Energy_constants[4], Energy_constants[5]};
    if (norm(n) == 0.0)
        throw std::invalid_argument("Plane_Adhesion: the plane normal (Z_Axis) is zero");
    normal_ = n / norm(n);
}

double Plane_Adhesion::Face_Energy(Face f)
{
    Halfedge he = f.halfedge();
    const Vector3 p0 = geometry->inputVertexPositions[he.vertex()];
    const Vector3 p1 = geometry->inputVertexPositions[he.next().vertex()];
    const Vector3 p2 = geometry->inputVertexPositions[he.next().next().vertex()];

    const double projected = -0.5 * dot(cross(p1 - p0, p2 - p0), normal_);
    if (projected <= 0.0)
        return 0.0;
    double dWdh;
    const double h = dot(normal_, (p0 + p1 + p2) / 3.0 - Bead_1->Pos);
    const double w = bead_geometry::planeContactWeight(h, width(), dWdh);
    return -strength() * w * projected;
}

void Plane_Adhesion::Add_Face_Force(Face f, VertexData<Vector3> &membraneForce, Vector3 &beadForce)
{
    Halfedge he = f.halfedge();
    const Vertex v0 = he.vertex(), v1 = he.next().vertex(), v2 = he.next().next().vertex();
    const Vector3 p0 = geometry->inputVertexPositions[v0];
    const Vector3 p1 = geometry->inputVertexPositions[v1];
    const Vector3 p2 = geometry->inputVertexPositions[v2];

    const double projected = -0.5 * dot(cross(p1 - p0, p2 - p0), normal_);
    if (projected <= 0.0)
        return;
    double dWdh;
    const double h = dot(normal_, (p0 + p1 + p2) / 3.0 - Bead_1->Pos);
    const double w = bead_geometry::planeContactWeight(h, width(), dWdh);
    if (w == 0.0)
        return;

    // d a_f / d p_i = 1/2 (p_(i+2) - p_(i+1)) x n,  d h_f / d p_i = n / 3
    // E = -W w a_f, force = -dE/dp = W (w da_f/dp + a_f dw/dh dh/dp)
    const Vector3 heightTerm = projected * dWdh * normal_ / 3.0;
    const Vector3 force0 = strength() * (w * 0.5 * cross(p2 - p1, normal_) + heightTerm);
    const Vector3 force1 = strength() * (w * 0.5 * cross(p0 - p2, normal_) + heightTerm);
    const Vector3 force2 = strength() * (w * 0.5 * cross(p1 - p0, normal_) + heightTerm);

    membraneForce[v0] += force0;
    membraneForce[v1] += force1;
    membraneForce[v2] += force2;

    // E depends only on p_i - P, so dE/dP = -sum_i dE/dp_i
    beadForce -= force0 + force1 + force2;
}

double Cilinder_Interaction::Tot_Energy()
{

    double Tot_E = 0;

    double rc = Energy_constants[2];
    double Zpos = Energy_constants[3];
    double Zc = Energy_constants[4];
    Vector3 Z_hat;
    Z_hat.x = Energy_constants[5];
    Z_hat.y = Energy_constants[6];
    Z_hat.z = Energy_constants[7];
    Eigen::Vector3d Z_axis;
    Eigen::Vector3d Z_orig;
    Z_orig << 0.0, 0.0, 1.0;
    Z_axis << Z_hat.x, Z_hat.y, Z_hat.z;
    // The next thing i need to do is to
    double s = (Z_orig.cross(Z_axis)).norm();
    double c = Z_axis[2];
    Eigen::Matrix3d Z_X = geometry->Cross_product_matrix(Z_orig.cross(Z_axis));
    Eigen::Matrix3d Rot = Eigen::Matrix3d::Identity() + Z_X + Z_X * Z_X * (1) / (1 + c);

    double faceArea;
    double rho;

    VertexData<Vector3> vector_bead(*mesh);
    VertexData<double> distance_bead(*mesh, 0.0);
    VertexData<double> r_Energy_contribution(*mesh, 0.0);
    VertexData<double> z_Energy_contribution(*mesh, 0.0);

    Vector3 Force_vector;
    Eigen::Vector3d VertPos;
    for (Vertex v : mesh->vertices())
    {

        vector_bead[v] = geometry->inputVertexPositions[v];
        VertPos << vector_bead[v].x, vector_bead[v].y, vector_bead[v].z;
        VertPos = Rot.transpose() * VertPos;

        distance_bead[v] = Zpos - VertPos[2];
        rho = sqrt(VertPos[0] * VertPos[0] + VertPos[1] * VertPos[1]);
        VertPos = Rot * Eigen::Vector3d(0.0, 0.0, distance_bead[v]);
        vector_bead[v] = Vector3({VertPos[0], VertPos[1], VertPos[2]});

        r_Energy_contribution[v] = E_r(rho, Energy_constants);

        z_Energy_contribution[v] = E_z(distance_bead[v], Energy_constants);

        Tot_E += geometry->vertexDualArea(v) * r_Energy_contribution[v] * z_Energy_contribution[v];

        // std::cout<<"The distance for vertex "<< v.getIndex() << " is " << distance_bead[v] << " and the energy contribution is " << Energy_contribution[v] << "\n";
    }
    return Tot_E;
}
VertexData<double> Cilinder_Interaction::V_Tot_Energy()
{
    VertexData<double> vertexE(*mesh, 0.0);
    return vertexE;
}

double Cilinder_Interaction::V_Energy(Vertex v)
{
    return 0.0;
}

VertexData<Vector3> Cilinder_Interaction::Gradient()
{
    VertexData<Vector3> Force(*mesh, {0.0, 0.0, 0.0});
    Bead_1->Prev_Total_force = Bead_1->Total_force;
    Bead_1->Total_force = {0.0, 0.0, 0.0};
    Bead_1->Total_force += Bead_1->CoverageForce;
    // Bead_1->Total_force +=Bond_force();
    // The first constant is the potential strength
    double rc = Energy_constants[2];
    double Zpos = Energy_constants[3];
    double Zc = Energy_constants[4];
    Vector3 Z_hat;
    Z_hat.x = Energy_constants[5];
    Z_hat.y = Energy_constants[6];
    Z_hat.z = Energy_constants[7];
    Eigen::Vector3d Z_axis;
    Eigen::Vector3d Z_orig;
    Z_orig << 0.0, 0.0, 1.0;
    Z_axis << Z_hat.x, Z_hat.y, Z_hat.z;
    // The next thing i need to do is to
    double s = (Z_orig.cross(Z_axis)).norm();
    double c = Z_axis[2];
    Eigen::Matrix3d Z_X = geometry->Cross_product_matrix(Z_orig.cross(Z_axis));
    Eigen::Matrix3d Rot = Eigen::Matrix3d::Identity() + Z_X + Z_X * Z_X * (1) / (1 + c);
    // So if i have the rotation i just need to

    // std::cout << "The rotation matrix is " << Rot << " \n";

    // Should i just start doing stuff? well i guess i can

    double faceArea;
    double rho;
    Halfedge he;
    Eigen::Vector<double, 6> Positions_r;
    Eigen::Vector<double, 9> Positions_f;
    Eigen::Vector<double, 6> Grad_r;
    Eigen::Vector<double, 9> Grad_f;
    std::array<int, 3> Vertices_triangle;
    VertexData<Vector3> vector_bead(*mesh);
    VertexData<double> distance_bead(*mesh, 0.0);
    VertexData<double> r_Energy_contribution(*mesh, 0.0);
    VertexData<double> dr_Energy_contribution(*mesh, 0.0);

    VertexData<double> z_Energy_contribution(*mesh, 0.0);
    VertexData<double> dz_Energy_contribution(*mesh, 0.0);

    Vector3 Force_vector;
    Eigen::Vector3d VertPos;
    // So the thing is that i
    // Rot.transposeInPlace()
    for (Vertex v : mesh->vertices())
    {

        vector_bead[v] = geometry->inputVertexPositions[v];
        VertPos << vector_bead[v].x, vector_bead[v].y, vector_bead[v].z;
        // VertPos = Rot.transpose() * VertPos;

        distance_bead[v] = Zpos - VertPos[2];
        rho = sqrt(VertPos[0] * VertPos[0] + VertPos[1] * VertPos[1]);
        VertPos = Rot * Eigen::Vector3d(0.0, 0.0, distance_bead[v]);
        vector_bead[v] = Vector3({VertPos[0], VertPos[1], VertPos[2]});

        r_Energy_contribution[v] = E_r(rho, Energy_constants);
        dr_Energy_contribution[v] = dE_r(rho, Energy_constants);

        z_Energy_contribution[v] = E_z(distance_bead[v], Energy_constants);
        dz_Energy_contribution[v] = dE_z(distance_bead[v], Energy_constants);
        // std::cout<<"The distance for vertex "<< v.getIndex() << " is " << distance_bead[v] << " and the energy contribution is " << Energy_contribution[v] << "\n";
    }

    for (Face f : mesh->faces())
    {
        faceArea = geometry->faceArea(f);
        he = f.halfedge();
        Vertices_triangle[0] = he.vertex().getIndex();
        Vertices_triangle[1] = he.next().vertex().getIndex();
        Vertices_triangle[2] = he.next().next().vertex().getIndex();
        Positions_f = stack_positions<3>(*geometry, Vertices_triangle);
        Grad_f = geometry->gradient_triangle_area(Positions_f);

        for (int i = 0; i < 3; i++)
        {
            // Now i need to
            Force_vector = {Grad_f[3 * i], Grad_f[3 * i + 1], Grad_f[3 * i + 2]};
            Force_vector *= (r_Energy_contribution[Vertices_triangle[0]] * z_Energy_contribution[Vertices_triangle[0]] + r_Energy_contribution[Vertices_triangle[1]] * z_Energy_contribution[Vertices_triangle[1]] + r_Energy_contribution[Vertices_triangle[2]] * z_Energy_contribution[Vertices_triangle[2]]) / 3.0;
            Force[Vertices_triangle[i]] -= Force_vector;
        }

        for (Vertex v : f.adjacentVertices())
        {
            if (v.isBoundary())
            {
                Force[v] = {0.0, 0.0, 0.0};
                continue;
            }
            double r = distance_bead[v];

            VertPos = {geometry->inputVertexPositions[v].x, geometry->inputVertexPositions[v].y, geometry->inputVertexPositions[v].z};

            VertPos = Rot.transpose() * VertPos;

            rho = sqrt(VertPos[0] * VertPos[0] + VertPos[1] * VertPos[1]);
            if (rho < rc || rc < 0.0)
            {

                // Positions_r << geometry->inputVertexPositions[v].x, geometry->inputVertexPositions[v].y, geometry->inputVertexPositions[v].z,
                // Bead_pos.x, Bead_pos.y, Bead_pos.z;
                Positions_r << 0.0, 0.0, 0.0,
                    vector_bead[v].x, vector_bead[v].y, vector_bead[v].z;
                // The positions r is a vector that has the coords of the vertex and the coords of the closes point
                //
                Grad_r = geometry->gradient_r(Positions_r, distance_bead[v]);
                Force_vector = Vector3{Grad_r[0], Grad_r[1], Grad_r[2]};
                Force_vector *= (1.0 / 3.0) * dz_Energy_contribution[v] * r_Energy_contribution[v] * faceArea;
                Force[v] -= Force_vector;
                Bead_1->Total_force += Force_vector;
            }
        }
    }

    return Force;
}