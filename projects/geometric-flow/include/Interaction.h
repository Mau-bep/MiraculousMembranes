#pragma once

#include <Eigen/Core>
#include <omp.h>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "Beads.h"

class Bead;

using namespace geometrycentral;
using namespace geometrycentral::surface;

class Interaction
{

public:
    ManifoldSurfaceMesh *mesh;
    VertexPositionGeometry *geometry;

    // The interaction needs to know which bead it is
    Bead *Bead_1;
    // What are the parameters for the interaction

    std::vector<double> Energy_constants;

    std::vector<std::string> E_Features;
    std::vector<int> E_Features_val;

    virtual double Bond_energy();
    virtual Vector3 Bond_force();
    virtual double Tot_Energy() = 0;
    virtual VertexData<double> V_Tot_Energy() = 0;
    virtual double V_Energy(Vertex v) = 0;
    virtual VertexData<Vector3> Gradient() = 0;
    virtual SparseMatrix<double> Hessian() = 0;
    virtual SparseMatrix<double> Hessian_IP() = 0;
    virtual std::vector<Eigen::Triplet<double>> Hessian_bonds_triplet();

    virtual double E_r(double r, std::vector<double> Energy_constants) = 0;
    virtual double dE_r(double r, std::vector<double> Energy_constants) = 0;
    virtual double ddE_r(double r, std::vector<double> Energy_constants) = 0;

    virtual double E_z(double r, std::vector<double> Energy_constants) = 0;
    virtual double dE_z(double r, std::vector<double> Energy_constants) = 0;

    // virtual Vector3 F_r(double r,Vector3 r_vec, std::vector<double> Energy_constants) = 0;
};

class No_mem_Inter : public Interaction
{
public:
    No_mem_Inter() {};

    double Tot_Energy() override
    {
        return Bond_energy();
    }
    VertexData<double> V_Tot_Energy()
    {
        VertexData<double> V(*mesh, 0.0);
        return V;
    }
    double V_Energy(Vertex v)
    {
        return 0.0;
    }
    virtual VertexData<Vector3> Gradient();
    virtual SparseMatrix<double> Hessian();
    virtual SparseMatrix<double> Hessian_IP();

    virtual double E_r(double r, std::vector<double> Energy_constants) override
    {
        return 0.0;
    }
    virtual double dE_r(double r, std::vector<double> Energy_constants) override
    {
        return 0.0;
    }
    virtual double ddE_r(double r, std::vector<double> Energy_constants) override
    {
        return 0.0;
    }
    virtual double E_z(double r, std::vector<double> Energy_constants) override
    {
        return 0.0;
    }
    virtual double dE_z(double r, std::vector<double> Energy_constants) override
    {
        return 0.0;
    }
};

class Cilinder_Interaction : public Interaction
{
public:
    Cilinder_Interaction() {};
    virtual VertexData<double> V_Tot_Energy();
    virtual double V_Energy(Vertex v);
    virtual double Tot_Energy();
    virtual VertexData<Vector3> Gradient();

    // Mode selector for the radial energy function. Supported values:
    // "Uniform" (default) and "ParabolaSection".
    std::string cylinder_E_mode = "Uniform";

    void setCylinderEMode(const std::string &mode)
    {
        cylinder_E_mode = mode;
    }

    virtual double Uniform(double r, std::vector<double> Energy_constants)
    {
        double rc = Energy_constants[2];
        if (r < 0.0 || r > rc)
            return 0.0;
        return 1.0;
    }
    virtual double d_Uniform(double r, std::vector<double> Energy_constants)
    {
        double rc = Energy_constants[2];
        if (r < 0.0 || r > rc)
            return 0.0;
        return 0.0;
    }
    virtual double ParabolaSection(double r, std::vector<double> Energy_constants)
    {
        // double eps = Energy_constants[0];
        double rc = Energy_constants[2];

        if (r < 0.0)
            return 0.0;

        if (r > 0 && r < rc / 2.0)
            return r * r * 2 / (rc * rc) - 1;

        if (r > rc / 2.0 && r < rc)
            return -r * r * 2 / (rc * rc) + r * 4 / rc - 2;

        return 0.0;
    }

    virtual double d_ParabolaSection(double r, std::vector<double> Energy_constants)
    {
        // double eps = Energy_constants[0];
        double rc = Energy_constants[2];
        if (r > 0 && r < rc / 2.0)
            return r * 4 / (rc * rc);

        if (r > rc / 2.0 && r < rc)
            return -r * 4 / (rc * rc) + 4 / rc;

        return 0.0;
    }

    virtual double Cubic(double r, std::vector<double> Energy_constants)
    {
        double rc = Energy_constants[2];
        if (r < 0.0)
            return 0.0;
        if (r > 0 && r < rc)
            return r * r * r * 2 / (rc * rc * rc) - r * r * 3 / (rc * rc) + 1;
        return 0.0;
    }

    virtual double d_Cubic(double r, std::vector<double> Energy_constants)
    {
        double rc = Energy_constants[2];
        if (r < 0.0)
            return 0.0;
        if (r > 0 && r < rc)
            return 6 * r * r / (rc * rc * rc) - 6 * r / (rc * rc);
        return 0.0;
    }

    SparseMatrix<double> Hessian()
    {
        SparseMatrix<double> Hessian(3, 3);
        return Hessian;
    }
    SparseMatrix<double> Hessian_IP()
    {
        SparseMatrix<double> Hessian(3, 3);
        return Hessian;
    }
    // Default dispatcher: picks the radial energy according to `cylinder_E_mode`.
    double E_r(double r, std::vector<double> Energy_constants) = 0;
    double dE_r(double r, std::vector<double> Energy_constants) = 0;
    double E_z(double r, std::vector<double> Energy_constants) = 0;
    double dE_z(double r, std::vector<double> Energy_constants) = 0;
};

class Plane_Interaction : public Cilinder_Interaction
{
public:
    double E_r(double r, std::vector<double> Energy_constants) override
    {
        // return ParabolaSection(r, Energy_constants);
        return Uniform(r, Energy_constants);
    }
    double dE_r(double r, std::vector<double> Energy_constants) override
    {
        // return d_ParabolaSection(r, Energy_constants);
        return d_Uniform(r, Energy_constants);
    }
    double ddE_r(double r, std::vector<double> Energy_constants) override
    {
        return 0.0;
    }

    virtual double E_z(double r, std::vector<double> Energy_constants) = 0;
    virtual double dE_z(double r, std::vector<double> Energy_constants) = 0;
};

class Gravity_Plane : public Plane_Interaction
{

public:
    Gravity_Plane() {}

    Gravity_Plane(std::vector<double> params)
    {
        Energy_constants = params;
    }
    Gravity_Plane(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, std::vector<double> params)
    {
        Energy_constants = params;
        mesh = inputMesh;
        geometry = inputGeo;
    }
    double E_z(double r, std::vector<double> Energy_constants) override
    {
        // WHat is the energy
        double eps = Energy_constants[0];
        return -1.0 * eps * r;
    }
    double dE_z(double r, std::vector<double> Energy_constants) override
    {
        double eps = Energy_constants[0];
        return -1.0 * eps;
    }
};

class Pinch_Interaction : public Cilinder_Interaction
{
public:
    Pinch_Interaction();
    Pinch_Interaction(std::vector<double> params)
    {
        Energy_constants = params;
    }
    Pinch_Interaction(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, std::vector<double> params)
    {
        Energy_constants = params;
        mesh = inputMesh;
        geometry = inputGeo;
    }

    double E_r(double r, std::vector<double> Energy_constants) override
    {
        return Cubic(r, Energy_constants);
    }
    double dE_r(double r, std::vector<double> Energy_constants) override
    {
        return d_Cubic(r, Energy_constants);
    }
    double ddE_r(double r, std::vector<double> Energy_constants) override
    {
        return 0.0;
    }
    double E_z(double r, std::vector<double> Energy_constants) override
    {
        // WHat is the energy
        double eps = Energy_constants[0];
        return -1.0 * eps * r;
    }
    double dE_z(double r, std::vector<double> Energy_constants) override
    {
        double eps = Energy_constants[0];
        return -1.0 * eps;
    }
};

class Integrated_Interaction : public Interaction
{
public:
    Integrated_Interaction() {};

    virtual double Tot_Energy();
    virtual VertexData<double> V_Tot_Energy();
    virtual double V_Energy(Vertex v);
    virtual VertexData<Vector3> Gradient();
    virtual SparseMatrix<double> Hessian();
    virtual SparseMatrix<double> Hessian_IP();

    virtual double E_r(double r, std::vector<double> Energy_constants) = 0;
    virtual double dE_r(double r, std::vector<double> Energy_constants) = 0;
    virtual double E_z(double r, std::vector<double> Energy_constants) override
    {
        return 0.0;
    }
    virtual double dE_z(double r, std::vector<double> Energy_constants) override
    {
        return 0.0;
    }

    virtual double ddE_r(double r, std::vector<double> Energy_constants) = 0;
    // virtual Vector3 F_r(double r, Vector3 r_vec, std::vector<double> Energy_constants) = 0;
};

class Normal_dot_Interaction : public Interaction
{

public:
    Normal_dot_Interaction() {};

    virtual double Tot_Energy();
    virtual VertexData<double> V_Tot_Energy();
    virtual double V_Energy(Vertex v);
    virtual VertexData<Vector3> Gradient();
    virtual SparseMatrix<double> Hessian();
    virtual SparseMatrix<double> Hessian_IP();
    virtual double E_r(double r, std::vector<double> Energy_constants) = 0;
    virtual double dE_r(double r, std::vector<double> Energy_constants) = 0;
    virtual double ddE_r(double r, std::vector<double> Energy_constants) = 0;
    // virtual Vector3 F_r(double r, Vector3 r_vec, std::vector<double> Energy_constants) = 0;
    virtual double E_z(double r, std::vector<double> Energy_constants) override
    {
        return 0.0;
    }
    virtual double dE_z(double r, std::vector<double> Energy_constants) override
    {
        return 0.0;
    }
};

class Frenkel : public Integrated_Interaction
{

public:
    Frenkel() {}

    Frenkel(std::vector<double> params)
    {
        Energy_constants = params;
    }
    Frenkel(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, std::vector<double> params)
    {
        Energy_constants = params;
        mesh = inputMesh;
        geometry = inputGeo;
    }

    double E_r(double r, std::vector<double> Energy_constants) override
    {
        // r is the distance between the two beads
        // Energy_constants[0] is the strength of the interaction
        // Energy_constants[1] is the sigma of the interaction
        double sigma = Energy_constants[1];
        double epsilon = Energy_constants[0];
        double rc = Energy_constants[2]; // cutoff distance

        if (r >= rc)
        {
            return 0.0; // No interaction beyond cutoff
        }
        double rc2 = rc * rc;
        double r2 = r * r;
        double alpha = 2 * (rc2 / (sigma * sigma)) * pow(3 / (2 * ((rc2 / (sigma * sigma)) - 1)), 3.0);

        return epsilon * alpha * ((sigma * sigma / r2) - 1) * pow((rc2 / r2) - 1, 2.0);
    }
    double dE_r(double r, std::vector<double> Energy_constants) override
    {

        double epsilon = Energy_constants[0];
        double sigma = Energy_constants[1];
        double rc = Energy_constants[2]; // cutoff distance
        if (r >= rc)
        {
            return 0.0; // No interaction beyond cutoff
        }
        double rc2 = rc * rc;
        double r2 = r * r;
        double alpha = 2 * (rc2 / (sigma * sigma)) * pow(3 / (2 * ((rc2 / (sigma * sigma)) - 1)), 3.0);

        // Ok so now i need to calculate the derivative
        double Q1 = -2 * alpha / (r2 * r);
        double Q2 = rc2 / r2 - 1;
        double Q3 = sigma * sigma * (rc2 / r2 - 1) + 2 * rc2 * (sigma * sigma / r2 - 1);

        return epsilon * Q1 * Q2 * Q3;
    }
    double ddE_r(double r, std::vector<double> Energy_constants) override
    {

        double sigma = Energy_constants[1];
        double epsilon = Energy_constants[0];
        double rc = Energy_constants[2]; // cutoff distance
        if (r >= rc)
        {
            return 0.0; // No interaction beyond cutoff
        }
        double rc2 = rc * rc;
        double r2 = r * r;
        double alpha = 2 * (rc2 / (sigma * sigma)) * pow(3 / (2 * ((rc2 / (sigma * sigma)) - 1)), 3.0);

        // Ok so now i need to calculate the second derivative
        double Q1 = -2 * alpha / (r2 * r);
        double Q2 = rc2 / r2 - 1;
        double Q3 = sigma * sigma * (rc2 / r2 - 1) + 2 * rc2 * (sigma * sigma / r2 - 1);

        double dQ1 = 6 * alpha / (r2 * r2);
        double dQ2 = -2 * rc2 / (r2 * r);
        double dQ3 = sigma * sigma * (-2) * rc2 / (r2 * r) - 4 * rc2 * sigma * sigma / (r2 * r);

        return dQ1 * Q2 * Q3 + Q1 * dQ2 * Q3 + Q1 * Q2 * dQ3;
    }
};

class Frenkel_Normal : public Normal_dot_Interaction
{

public:
    Frenkel_Normal() {}

    Frenkel_Normal(std::vector<double> params)
    {
        Energy_constants = params;
    }
    Frenkel_Normal(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, std::vector<double> params)
    {
        Energy_constants = params;
        mesh = inputMesh;
        geometry = inputGeo;
    }

    double E_r(double r, std::vector<double> Energy_constants) override
    {
        // r is the distance between the two beads
        // Energy_constants[0] is the strength of the interaction
        // Energy_constants[1] is the sigma of the interaction
        double sigma = Energy_constants[1];
        double epsilon = Energy_constants[0];
        double rc = Energy_constants[2]; // cutoff distance

        if (r >= rc)
        {
            return 0.0; // No interaction beyond cutoff
        }
        double rc2 = rc * rc;
        double r2 = r * r;
        double alpha = 2 * (rc2 / (sigma * sigma)) * pow(3 / (2 * ((rc2 / (sigma * sigma)) - 1)), 3.0);

        return epsilon * alpha * ((sigma * sigma / r2) - 1) * pow((rc2 / r2) - 1, 2.0);
    }
    double dE_r(double r, std::vector<double> Energy_constants) override
    {

        double epsilon = Energy_constants[0];
        double sigma = Energy_constants[1];
        double rc = Energy_constants[2]; // cutoff distance
        if (r >= rc)
        {
            return 0.0; // No interaction beyond cutoff
        }
        double rc2 = rc * rc;
        double r2 = r * r;
        double alpha = 2 * (rc2 / (sigma * sigma)) * pow(3 / (2 * ((rc2 / (sigma * sigma)) - 1)), 3.0);

        // Ok so now i need to calculate the derivative
        double Q1 = -2 * alpha / (r2 * r);
        double Q2 = rc2 / r2 - 1;
        double Q3 = sigma * sigma * (rc2 / r2 - 1) + 2 * rc2 * (sigma * sigma / r2 - 1);

        return epsilon * Q1 * Q2 * Q3;
    }
    double ddE_r(double r, std::vector<double> Energy_constants) override
    {

        double sigma = Energy_constants[1];
        double epsilon = Energy_constants[0];
        double rc = Energy_constants[2]; // cutoff distance
        if (r >= rc)
        {
            return 0.0; // No interaction beyond cutoff
        }
        double rc2 = rc * rc;
        double r2 = r * r;
        double alpha = 2 * (rc2 / (sigma * sigma)) * pow(3 / (2 * ((rc2 / (sigma * sigma)) - 1)), 3.0);

        // Ok so now i need to calculate the second derivative
        double Q1 = -2 * alpha / (r2 * r);
        double Q2 = rc2 / r2 - 1;
        double Q3 = sigma * sigma * (rc2 / r2 - 1) + 2 * rc2 * (sigma * sigma / r2 - 1);

        double dQ1 = 6 * alpha / (r2 * r2);
        double dQ2 = -2 * rc2 / (r2 * r);
        double dQ3 = sigma * sigma * (-2) * rc2 / (r2 * r) - 4 * rc2 * sigma * sigma / (r2 * r);

        return dQ1 * Q2 * Q3 + Q1 * dQ2 * Q3 + Q1 * Q2 * dQ3;
    }
};

class Constant_Normal : public Normal_dot_Interaction
{
public:
    Constant_Normal(std::vector<double> params)
    {
        Energy_constants = params;
    }
    double E_r(double r, std::vector<double> Energy_constants) override
    {

        return 1.0;
    }
    double dE_r(double r, std::vector<double> Energy_constants) override
    {

        return 0.0; // Constant interaction, no gradient
    }
    double ddE_r(double r, std::vector<double> Energy_constants) override
    {

        return 0.0; // Constant interaction, no second derivative
    }
};
class One_over_r_Normal : public Normal_dot_Interaction
{
public:
    One_over_r_Normal() {}

    One_over_r_Normal(std::vector<double> params)
    {
        Energy_constants = params;
    }

    One_over_r_Normal(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, std::vector<double> params)
    {
        Energy_constants = params;
        mesh = inputMesh;
        geometry = inputGeo;
    }

    double E_r(double r, std::vector<double> Energy_constants) override
    {
        //  1/r =
        // Which are the energy constants
        double epsilon = Energy_constants[0];

        return epsilon / r + fabs(epsilon);
    }

    double dE_r(double r, std::vector<double> Energy_constants) override
    {
        double epsilon = Energy_constants[0];

        return -epsilon / (r * r);
    }

    double ddE_r(double r, std::vector<double> Energy_constants) override
    {
        double epsilon = Energy_constants[0];

        return 2 * epsilon / (r * r * r);
    }
};

class Linear_Normal : public Normal_dot_Interaction
{
public:
    Linear_Normal() {}

    Linear_Normal(std::vector<double> params)
    {
        Energy_constants = params;
    }

    Linear_Normal(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, std::vector<double> params)
    {
        Energy_constants = params;
        mesh = inputMesh;
        geometry = inputGeo;
    }

    double E_r(double r, std::vector<double> Energy_constants) override
    {
        //  1/r =
        // Which are the energy constants
        double epsilon = Energy_constants[0];

        return epsilon * r + fabs(epsilon);
    }

    double dE_r(double r, std::vector<double> Energy_constants) override
    {
        double epsilon = Energy_constants[0];

        return epsilon;
    }

    double ddE_r(double r, std::vector<double> Energy_constants) override
    {
        double epsilon = Energy_constants[0];

        return 0;
    }
};

class One_over_r : public Integrated_Interaction
{
public:
    One_over_r() {}

    One_over_r(std::vector<double> params)
    {
        Energy_constants = params;
    }

    One_over_r(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, std::vector<double> params)
    {
        Energy_constants = params;
        mesh = inputMesh;
        geometry = inputGeo;
    }

    double E_r(double r, std::vector<double> Energy_constants) override
    {
        //  1/r =
        // Which are the energy constants
        double epsilon = Energy_constants[0];

        return epsilon / r;
    }

    double dE_r(double r, std::vector<double> Energy_constants) override
    {
        double epsilon = Energy_constants[0];

        return -epsilon / (r * r);
    }

    double ddE_r(double r, std::vector<double> Energy_constants) override
    {
        double epsilon = Energy_constants[0];

        return 2 * epsilon / (r * r * r);
    }
};

class Linear : public Integrated_Interaction
{
public:
    Linear() {}

    Linear(std::vector<double> params)
    {
        Energy_constants = params;
    }

    Linear(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, std::vector<double> params)
    {
        Energy_constants = params;
        mesh = inputMesh;
        geometry = inputGeo;
    }

    double E_r(double r, std::vector<double> Energy_constants) override
    {
        //  1/r =
        // Which are the energy constants
        double epsilon = Energy_constants[0];

        return epsilon * r;
    }

    double dE_r(double r, std::vector<double> Energy_constants) override
    {
        double epsilon = Energy_constants[0];

        return epsilon;
    }

    double ddE_r(double r, std::vector<double> Energy_constants) override
    {
        double epsilon = Energy_constants[0];

        return 0;
    }
};

class LJ_Normal : public Normal_dot_Interaction
{

public:
    LJ_Normal() {}

    LJ_Normal(std::vector<double> params)
    {
        Energy_constants = params;
    }
    LJ_Normal(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, std::vector<double> params)
    {
        Energy_constants = params;
        mesh = inputMesh;
        geometry = inputGeo;
    }

    double E_r(double r, std::vector<double> Energy_constants) override
    {
        // r is the distance between the two beads
        // Energy_constants[0] is the strength of the interaction
        // Energy_constants[1] is the sigma of the interaction
        double sigma = Energy_constants[1];
        double epsilon = Energy_constants[0];
        double rc = Energy_constants[2];    // cutoff distance
        double shift = Energy_constants[3]; // shift value
        if (r >= rc)
        {
            return 0.0; // No interaction beyond cutoff
        }
        return 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6)) + shift;
    }
    double dE_r(double r, std::vector<double> Energy_constants) override
    {

        // r is the distance between the two beads
        // Energy_constants[0] is the strength of the interaction
        // Energy_constants[1] is the sigma of the interaction

        double sigma = Energy_constants[1];
        double epsilon = Energy_constants[0];
        double rc = Energy_constants[2]; // cutoff distance
        if (r >= rc)
        {
            return 0.0; // No interaction beyond cutoff
        }
        return -24 * epsilon * (2 * pow(sigma / r, 12) - pow(sigma / r, 6)) / r;
    }
    double ddE_r(double r, std::vector<double> Energy_constants) override
    {
        // r is the distance between the two beads
        // Energy_constants[0] is the strength of the interaction
        // Energy_constants[1] is the sigma of the interaction
        double sigma = Energy_constants[1];
        double epsilon = Energy_constants[0];
        double rc = Energy_constants[2]; // cutoff distance
        if (r >= rc)
        {
            return 0.0; // No interaction beyond cutoff
        }
        return 24 * epsilon * (2 * pow(sigma / r, 12) - pow(sigma / r, 6)) / (r * r) + 24 * epsilon * (24 * pow(sigma / r, 12) - 6 * pow(sigma / r, 6)) / (r * r);
    }
};

class LJ : public Integrated_Interaction
{
public:
    LJ() {}

    LJ(std::vector<double> params)
    {
        Energy_constants = params;
    }
    LJ(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, std::vector<double> params)
    {
        Energy_constants = params;
        mesh = inputMesh;
        geometry = inputGeo;
    }

    double E_r(double r, std::vector<double> Energy_constants) override
    {
        // r is the distance between the two beads
        double epsilon = Energy_constants[0];
        double sigma = Energy_constants[1];
        double rc = Energy_constants[2];    // cutoff distance
        double shift = Energy_constants[3]; // shift value
        if (r >= rc && rc > 0.0)
        {
            return 0.0; // No interaction beyond cutoff
        }
        return 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6)) + shift;
    }
    double dE_r(double r, std::vector<double> Energy_constants) override
    {

        // r is the distance between the two beads
        // Energy_constants[0] is the strength of the interaction
        // Energy_constants[1] is the sigma of the interaction
        double sigma = Energy_constants[1];
        double epsilon = Energy_constants[0];
        double rc = Energy_constants[2]; // cutoff distance
        if (r >= rc && rc > 0)
        {
            return 0.0; // No interaction beyond cutoff
        }
        return -24 * epsilon * (2 * pow(sigma / r, 12) - pow(sigma / r, 6)) / r;
    }
    double ddE_r(double r, std::vector<double> Energy_constants) override
    {
        // r is the distance between the two beads
        // Energy_constants[0] is the strength of the interaction
        // Energy_constants[1] is the sigma of the interaction
        double sigma = Energy_constants[1];
        double epsilon = Energy_constants[0];
        double rc = Energy_constants[2]; // cutoff distance
        if (r >= rc && rc > 0)
        {
            return 0.0; // No interaction beyond cutoff
        }
        return 24 * epsilon * (2 * pow(sigma / r, 12) - pow(sigma / r, 6)) / (r * r) + 24 * epsilon * (24 * pow(sigma / r, 12) - 6 * pow(sigma / r, 6)) / (r * r);
    }
};
