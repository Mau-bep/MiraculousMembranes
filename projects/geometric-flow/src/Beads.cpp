

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


double Bead::Energy()
{

    double Total_E = 0.0;
    double r;

    // I need to add the energy of the interaction between beads
    for (size_t bead = 0; bead < Beads.size(); bead++)
    {
        vector<double> params = Interaction_constants_vector[bead];
        // Ok so i loaded de params of the interaction
        if (Bond_type[bead] == "Lineal" && state != "manual" && state != "froze")
        {
            Total_E += params[0] * (Pos - Beads[bead]->Pos).norm();
        }

        if (Bond_type[bead] == "Harmonic" && state != "manual" && state != "froze")
        {
            // The armonic interaction has two parameters (stiffness and rest length) assuming first is the sitffness and there is restlength -
            Total_E += params[0] * dot(Pos - Beads[bead]->Pos, Pos - Beads[bead]->Pos) / 2.0;
            // std::cout<<"ADDed harmonic energy"<< Total_E <<"\n";
        }
        if (Bond_type[bead] == "Shifted-LJ" && state != "manual" && state != "froze")
        {
            // The shifted LJ has two parameters (stiffness and rest length) assuming first is the sitffness and there is restlength -
            double dist = (Pos - Beads[bead]->Pos).norm();
            double threshold = sigma * pow(2.0, 1 / 6.0);
            if (dist < threshold)
            {
                Total_E += params[0] * 4 * (pow((sigma / dist), 12) - pow((sigma / dist), 6)) - params[0] * 4 * (pow((sigma / threshold), 12) - pow((sigma / threshold), 6));
            }
        }
    }
    // std::cout<<"THe interaction is " << interaction <<" \n";
    // std::cout<<"Total E is " << Total_E <<" \n";

    double dual_area;
    if (interaction == "Spring")
    {

        for (Vertex v : mesh->vertices())
        {

            r = (this->Pos - geometry->inputVertexPositions[v]).norm();

            // dual_area=geometry->barycentricDualArea(v);

            Total_E += 0.5 * strength * (r - sigma) * (r - sigma);
        }
        return Total_E;
    }

    if (interaction == "Linear_field")
    {
        Total_E = 0;
        for (Vertex v : mesh->vertices())
        {
            dual_area = geometry->barycentricDualArea(v);
            Total_E += dual_area * strength * ((Pos - geometry->inputVertexPositions[v]).norm());
        }

        return Total_E;
    }
    if (interaction == "Gravity")
    {
        Total_E = 0;
        for (Vertex v : mesh->vertices())
        {
            dual_area = geometry->barycentricDualArea(v);
            Total_E += dual_area * strength * (-1 * (geometry->inputVertexPositions[v].x));
        }

        return Total_E;
    }
    if (interaction == "One_over_r_x")
    {

        for (Vertex v : mesh->vertices())
        {
            dual_area = geometry->barycentricDualArea(v);
            Total_E += dual_area * strength * (1 / (geometry->inputVertexPositions[v].x - Pos.x));
        }

        return Total_E;
    }
    if (interaction == "Shifted-LJ")
    {
        // double dual_area;
        double alpha;
        double rc2 = rc * rc;
        double r2;
        for (Vertex v : mesh->vertices())
        {
            r2 = (this->Pos - geometry->inputVertexPositions[v]).norm2();
            if (r2 > rc2)
            {
                continue;
            }

            dual_area = geometry->barycentricDualArea(v);
            // dual_area=geometry->circumcentricDualArea(v);

            Total_E += dual_area * (4 * strength * (pow(sigma * sigma / r2, 6) - pow((sigma * sigma / r2), 3)) + strength);

            // Total_E += strength*alpha*( (sigma*sigma/r2)-1  )*pow( (rc2/r2)-1 ,2.0);
        }
        // std::cout<<"THe total energy is "<<Total_E <<" \n";
        return Total_E;
    }

    if (interaction == "LJ")
    {
        for (Vertex v : mesh->vertices())
        {
            r = (this->Pos - geometry->inputVertexPositions[v]).norm();
            // dual_area=geometry->circumcentricDualArea(v);
            dual_area = geometry->barycentricDualArea(v);

            // Total_E+= 4*strength*(pow(sigma/r,12)-pow(sigma/r,6));
            // Total_E+=2*dual_area;

            // if(v.getIndex()==3){
            Total_E += 4 * dual_area * strength * (pow(sigma / r, 12) - pow(sigma / r, 6));

            //     std::cout<<" The radial distance is "<<r <<"the dual area is "<< dual_area <<" and the total computed energy is "<<Total_E<<" \n";

            // }
        }
        return Total_E;
    }

    if (interaction == "test_angle")
    {

        Vertex v_1 = mesh->vertex(1);
        Corner c = v_1.halfedge().corner();

        Total_E += geometry->angle(c);
    }

    if (interaction == "test_normal")
    {

        // for(Vertex v : mesh->vertices()){
        for (Face f : mesh->faces())
        {
            Total_E += dot(geometry->faceNormal(f), sqrt(3) * Vector3({0, 1, 0}));
        }
        // }
    }
    if (interaction == "test_angle_normal")
    {

        // Total_E=dot(geometry->vertexNormalAngleWeighted(mesh->vertex(3)),sqrt(3)*Vector3({0,1,0}));
        // for(Vertex v : mesh->vertices()){
        for (Vertex v : mesh->vertices())
        {
            Total_E += dot(geometry->vertexNormalAngleWeighted(v), sqrt(3) * Vector3({0, 1, 0}));
        }
        return Total_E;
        // }
    }
    if (interaction == "test_angle_normal_r")
    {

        Vector3 unit_r;
        for (Vertex v : mesh->vertices())
        {
            unit_r = (Pos - geometry->inputVertexPositions[v]).unit();
            Total_E += dot(geometry->vertexNormalAngleWeighted(v), unit_r);
        }
        return Total_E;
        // }
    }

    if (interaction == "test_unit_r")
    {
        Vector3 unit_r;

        for (Vertex v : mesh->vertices())
        {

            unit_r = (Pos - geometry->inputVertexPositions[v]).unit();
            Total_E += dot(unit_r, sqrt(3) * Vector3({1, 1, 1}));
        }
        return Total_E;
    }

    if (interaction == "test_angle_normal_r_normalized")
    {
        double val;
        Vector3 unit_r;
        Vector3 Angle_normal;

        for (Vertex v : mesh->vertices())
        {
            // Vertex v = mesh->vertex(3);
            unit_r = (Pos - geometry->inputVertexPositions[v]).unit();

            Angle_normal = geometry->vertexNormalAngleWeighted(v);
            // std::cout<<"Dot product is "<< dot(unit_r,Angle_normal)<<"\n";
            if (dot(unit_r, Angle_normal) > 0)
            {

                Total_E += dot(unit_r, Angle_normal) / Angle_normal.norm();
            }
        }

        return Total_E;
    }

    if (interaction == "test_angle_normal_r_normalized_LJ")
    {
        double val;
        Vector3 unit_r;
        Vector3 Angle_normal;
        double r2;
        double rc2 = rc * rc;
        double alpha;
        for (Vertex v : mesh->vertices())
        {
            // Vertex v = mesh->vertex(3);

            unit_r = (Pos - geometry->inputVertexPositions[v]).unit();

            Angle_normal = geometry->vertexNormalAngleWeighted(v);
            // std::cout<<"Dot product is "<< dot(unit_r,Angle_normal)<<"\n";
            if (dot(unit_r, Angle_normal) > 0)
            {

                unit_r = Pos - geometry->inputVertexPositions[v];
                r2 = unit_r.norm2();
                unit_r = unit_r.unit();

                if (r2 < rc2)
                {

                    alpha = 2 * (rc2 / (sigma * sigma)) * pow(3 / (2 * ((rc2 / (sigma * sigma)) - 1)), 3.0);
                    dual_area = geometry->barycentricDualArea(v);

                    Total_E += dual_area * strength * alpha * ((sigma * sigma / r2) - 1) * pow((rc2 / r2) - 1, 2.0) * dot(unit_r, Angle_normal) / Angle_normal.norm();
                }
            }
        }

        return Total_E;
    }

    if (interaction == "Shifted_LJ_Normal")
    {
        double val;
        Vector3 unit_r;
        Vector3 Angle_normal;
        double r2;
        double r;
        double rc2 = rc * rc;
        double alpha;
        for (Vertex v : mesh->vertices())
        {

            Angle_normal = geometry->vertexNormalAngleWeighted(v);

            unit_r = Pos - geometry->inputVertexPositions[v];
            r2 = unit_r.norm2();
            r = unit_r.norm();
            unit_r = unit_r / r;

            if (r < rc && dot(unit_r, Angle_normal) > 0)
            {

                alpha = 2 * (rc2 / (sigma * sigma)) * pow(3 / (2 * ((rc2 / (sigma * sigma)) - 1)), 3);
                dual_area = geometry->barycentricDualArea(v);
                // std::cout<<"We are here at some point right?\n";
                Total_E += dual_area * strength * alpha * ((sigma * sigma / r2) - 1) * pow((rc2 / r2) - 1, 2) * dot(unit_r, Angle_normal) / Angle_normal.norm();
            }
        }
        return Total_E;
    }

    // std::cout<<"Energy is here still" << Total_E << " \n";
    if (interaction == "Shifted_LJ_Normal_nopush" || interaction == "Frenkel_Normal_nopush")
    {
        double val;
        Vector3 unit_r;
        Vector3 unit_r2;
        Vector3 unit_r3;

        Vector3 Angle_normal;
        Vector3 Face_normal;
        double face_area;
        double r2;
        double r;
        double rc2 = rc * rc;
        double alpha;
        for (Face f : mesh->faces())
        {
            Face_normal = geometry->faceNormal(f);
            face_area = geometry->faceArea(f);
            for (Vertex v : f.adjacentVertices())
            {
                unit_r = Pos - geometry->inputVertexPositions[v];
                r = unit_r.norm();
                unit_r = unit_r.unit();
                r2 = r * r;
                if (r < rc && dot(Face_normal, unit_r) > 0)
                {
                    // if(r<rc){

                    alpha = 2 * (rc2 / (sigma * sigma)) * pow(3 / (2 * ((rc2 / (sigma * sigma)) - 1)), 3);
                    // Total_E+=(face_area/3.0);
                    // Total_E+=(face_area/3.0)*strength*alpha*( (sigma*sigma/r2)-1  )*pow( (rc2/r2)-1 ,2);
                    Total_E += (face_area / 3.0) * strength * alpha * ((sigma * sigma / r2) - 1) * pow((rc2 / r2) - 1, 2) * dot(Face_normal, unit_r);
                    // Total_E+=dot(Face_normal,unit_r);
                }
            }
        }

        return Total_E;
    }

    if (interaction == "Shifted_LJ_Normal_nopush_inside" || interaction == "Frenkel_Normal_nopush_inside")
    {
        double val;
        Vector3 unit_r;
        Vector3 unit_r2;
        Vector3 unit_r3;

        Vector3 Angle_normal;
        Vector3 Face_normal;
        double face_area;
        double r2;
        double r;
        double rc2 = rc * rc;
        double alpha;
        for (Face f : mesh->faces())
        {
            Face_normal = geometry->faceNormal(f);
            face_area = geometry->faceArea(f);
            for (Vertex v : f.adjacentVertices())
            {
                unit_r = Pos - geometry->inputVertexPositions[v];
                r = unit_r.norm();
                unit_r = unit_r.unit();
                r2 = r * r;
                if (r < rc && dot(Face_normal, unit_r) < 0)
                {
                    // if(r<rc){

                    alpha = 2 * (rc2 / (sigma * sigma)) * pow(3 / (2 * ((rc2 / (sigma * sigma)) - 1)), 3);
                    // Total_E+=(face_area/3.0);
                    // Total_E+=(face_area/3.0)*strength*alpha*( (sigma*sigma/r2)-1  )*pow( (rc2/r2)-1 ,2);
                    Total_E += -1 * (face_area / 3.0) * strength * alpha * ((sigma * sigma / r2) - 1) * pow((rc2 / r2) - 1, 2) * dot(Face_normal, unit_r);
                    // Total_E+=dot(Face_normal,unit_r);
                }
            }
        }

        return Total_E;
    }
    // std::cout<<"Should never be here\n";

    if (interaction == "Shifted_LJ_Normal_var")
    {
        double val;
        Vector3 unit_r;
        Vector3 unit_r2;
        Vector3 unit_r3;

        Vector3 Angle_normal;
        Vector3 Face_normal;
        double face_area;
        double r2;
        double r;
        double rc2 = rc * rc;
        double alpha;

        for (Face f : mesh->faces())
        {
            Face_normal = geometry->faceNormal(f);
            face_area = geometry->faceArea(f);
            for (Vertex v : f.adjacentVertices())
            {
                unit_r = Pos - geometry->inputVertexPositions[v];
                r = unit_r.norm();
                unit_r = unit_r / r;
                r2 = r * r;
                // if(r<rc && dot(Face_normal,unit_r)>0){
                if (r < rc)
                {
                    // std::cout<<"THe POS of the vertex is "<< geometry->inputVertexPositions[v]<<" \n";

                    alpha = 2 * (rc2 / (sigma * sigma)) * pow(3 / (2 * ((rc2 / (sigma * sigma)) - 1)), 3);
                    // Total_E+=(face_area/3.0);
                    // Total_E+=(face_area/3.0)*strength*alpha*( (sigma*sigma/r2)-1  )*pow( (rc2/r2)-1 ,2);
                    Total_E += (face_area / 3.0) * strength * alpha * ((sigma * sigma / r2) - 1) * pow((rc2 / r2) - 1, 2) * dot(Face_normal, unit_r);
                    // Total_E+=dot(Face_normal,unit_r);
                }
            }
        }

        return Total_E;
    }

    if (interaction == "test_angle_normal_r_normalized_LJ_Full")
    {
        double val;
        Vector3 unit_r;
        Vector3 Angle_normal;
        double r2;
        double r;
        double rc2 = rc * rc;
        double alpha;
        for (Vertex v : mesh->vertices())
        {

            Angle_normal = geometry->vertexNormalAngleWeighted(v);
            // std::cout<<"Dot product is "<< dot(unit_r,Angle_normal)<<"\n";
            unit_r = Pos - geometry->inputVertexPositions[v];
            r2 = unit_r.norm2();
            r = unit_r.norm();
            unit_r = unit_r / r;

            if (r < rc)
            {

                // if(r2<rc2){

                alpha = 2 * (rc2 / (sigma * sigma)) * pow(3 / (2 * ((rc2 / (sigma * sigma)) - 1)), 3);
                dual_area = geometry->barycentricDualArea(v);

                Total_E += dual_area * strength * alpha * ((sigma * sigma / r2) - 1) * pow((rc2 / r2) - 1, 2) * dot(unit_r, Angle_normal) / Angle_normal.norm();

                // }
            }
        }

        return Total_E;
    }
    if (interaction == "pulling")
    {
        double r1 = 0.5;
        double r2 = 2.0;
        Vector3 unit_r;
        for (Vertex v : mesh->vertices())
        {
            unit_r = this->Pos - geometry->inputVertexPositions[v];
            r = unit_r.norm();
            dual_area = geometry->barycentricDualArea(v);
            if (r < r1)
            {
                Total_E += dual_area * r + strength * r1 - strength * (r1 - r2) / 2;
            }
            if (r > r1 && r < r2)
            {
                Total_E += dual_area * (strength / (r1 - r2)) * (r - r2) * (r - r2) / 2;
            }
        }
        return Total_E;
    }

    return Total_E;
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