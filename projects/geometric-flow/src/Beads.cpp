

#include "Beads.h"
#include <fstream>
#include <omp.h>
#include <chrono>
#include <Eigen/Core>
#include "Interaction.h"

using namespace std;
/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
Bead::Bead(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, Vector3 Position, double input_sigma, double strg)
{

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
    Pos = Position;
    sigma = input_sigma;
    strength = strg;
    interaction = "Shifted_LJ_Normal_var";
    state = "default"; // The other two states are manual and frozee
    Velocity = Vector3({0, 0, 0});
    rc = sigma * 2.0;
    prev_force = 0.0;
    Total_force = {0, 0, 0};
    prev_E_stationary = 0;
}

Bead::Bead(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, Vector3 Position, double input_sigma, double strg, double input_rc)
{

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
    Pos = Position;
    sigma = input_sigma;
    strength = strg;
    interaction = "Shifted_LJ_Normal_var";
    state = "default"; // The other two states are manual and frozee
    Velocity = Vector3({0, 0, 0});
    rc = input_rc;
    prev_force = 0.0;
    Total_force = {0, 0, 0};
    prev_E_stationary = 0;
}
Bead::Bead(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, Vector3 Position, std::vector<double> Energy_constants, Interaction *Interact, int id, int Number_beads)
{
    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
    Pos = Position;
    sigma = Energy_constants[1];
    strength = Energy_constants[0];
    rc = Energy_constants[2];
    interaction = "";

    Bead_I = Interact;
    Bead_I->Bead_1 = this;
    Bead_id = id;
    Total_beads = Number_beads;

    state = "default"; // The other two states are manual and frozee
    Velocity = Vector3({0, 0, 0});

    prev_force = 0.0;
    Total_force = {0, 0, 0};
    prev_E_stationary = 0;
}

void Bead::Add_bead(Bead *bead, std::string Interaction, vector<double> Interaction_constants)
{

    Beads.push_back(bead);
    Bond_type.push_back(Interaction);
    Interaction_constants_vector.push_back(Interaction_constants);
}
void Bead::Reasign_mesh(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo)
{
    mesh = inputMesh;
    geometry = inputGeo;
    Bead_I->mesh = inputMesh;
    Bead_I->geometry = inputGeo;
    Bead_I->Bead_1 = this;
}


// The membrane gradient of this bead's interaction
VertexData<Vector3> Bead::Gradient()
{
    return Bead_I->Gradient();
}

void Bead::Reset_bead(Vector3 Actual_pos)
{
    this->Pos = Actual_pos;
    return;
}

void Bead::Move_bead(double dt, Vector3 center)
{
    if (state == "default")
    {
        if (Constraint == "None")
        {
            this->Pos = this->Pos + Total_force * dt - center;
        }
        if (Constraint == "Radial")
        {

            double theta = Constraint_constants[0];
            Vector3 unit_r = Vector3{sin(theta), 0.0, cos(theta)};
            // std::cout<<"THe angle is " << theta << " and the vector is " << unit_r<<" \n";
            // std::cout<<"THe displacement is in " << unit_r*(dot(Total_force,unit_r)) <<" \n";
            this->Pos = this->Pos + unit_r * (dot(Total_force, unit_r) * dt) - center;
            // std::cout<<"THe center is " << center << "\n";
        }
    }
    if (state == "froze")
        return;

    if (state == "manual")
    {
        if (dot(Total_force, Velocity) > 0)
        {
            // The force is in the direction of the velocity, lets move
            this->Pos += dot(Total_force, Velocity) * Velocity.unit() * dt;
        }
        else
        {
            this->Pos += Velocity * dt;
            // std::cout << "\t The distance to the final position is " << (FinalPos - Pos).norm() << "\n";
            // std::cout << "\t The dot product is " << dot(Total_force, Velocity) << " cute \n";
        }
        // If the force is not align, we dont move
        // I want do do manual in a better way something like  dot(Total_force,Velocity)*Velocity.unit()  cause it should be nicer to the integration
        // std::cout<<"Moving manually\n";
    }
    // std::cout<<"The bead is moving "<< norm(Total_force*dt -center)<<"\n";
    return;
}

void Bead::Move_bead(double dt, Vector3 center, Vector3 Force)
{
    if (state == "default")
    {
        if (Constraint == "None")
        {
            this->Pos = this->Pos + Force * dt - center;
        }
        if (Constraint == "Radial")
        {

            double theta = Constraint_constants[0];
            Vector3 unit_r = Vector3{sin(theta), 0.0, cos(theta)};
            // std::cout<<"THe angle is " << theta << " and the vector is " << unit_r<<" \n";
            // std::cout<<"THe displacement is in " << unit_r*(dot(Total_force,unit_r)) <<" \n";
            this->Pos = this->Pos + unit_r * (dot(Force, unit_r) * dt) - center;
            // std::cout<<"THe center is " << center << "\n";
        }
    }
    if (state == "froze")
        return;

    if (state == "manual")
    {

        if (dot(Total_force, Velocity) > 0)
        {
            // The force is in the direction of the velocity, lets move
            this->Pos += dot(Total_force, Velocity) * Velocity.unit() * dt;
        }
        else
        {
            this->Pos += Velocity * dt;
            // std::cout << "\t The distance to the final position is " << (FinalPos - Pos).norm() << "\n";
            // std::cout << "\t The dot product is " << dot(Total_force, Velocity) << " cute \n";
        }
    }

    // std::cout<<"The bead is moving "<< norm(Total_force*dt -center)<<"\n";
    return;
}


void Bead::update_state()
{
    // Now
    if (state == "manual")
    {
        if ((Pos - FinalPos).norm() < 0.01)
        {
            state = "froze";
            std::cout << "Bead has reached the final position, freezing it\n";
            // Now i want to also remove the bond energy
            Pos = FinalPos;
            if (Beads.size() > 0)
            {
                // aLSO BREAKING THE BONDS
                Beads.resize(0);
                // for (size_t i = 0; i < Beads.size(); i++)
                // {
                //     Interaction_constants_vector[i][0] = 0;
                // }
            }
        }
    }
}