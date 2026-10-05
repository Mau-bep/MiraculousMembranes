#pragma once

#include <Eigen/Core>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "Beads.h"
// #include "Mem-3dg.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class E_Handler
{
public:
    ManifoldSurfaceMesh *mesh = nullptr;
    VertexPositionGeometry *geometry = nullptr;

    std::vector<std::string> Energies;
    std::vector<std::vector<double>> Energy_constants;
    std::vector<std::string> Constraints;
    Eigen::VectorXd Lagrange_mult;
    int N_constraints = 0;
    double Trgt_vol = 0.0;
    double Trgt_area = 0.0;
    bool boundary = false;
    std::vector<Bead *> Beads;

    std::vector<double> Energy_values;

    VertexData<Vector3> Previous_grad;
    VertexData<Vector3> Current_grad;
    VertexData<Vector3> Vertex_normals;
    FaceData<double> Face_reference;

    Eigen::MatrixXd Jacobian_constraints;

    std::vector<double> Gradient_norms;

    // Cost accounting for Timing.txt: number of calls, and wall time in total
    // and per entry of Energies
    long energy_calls = 0;
    long gradient_calls = 0;
    double energy_ms = 0.0;
    double gradient_ms = 0.0;
    std::vector<double> term_energy_ms;
    std::vector<double> term_gradient_ms;

    // Memo of Calculate_energies: when nothing the energies read has changed
    // since the last call (vertex positions, connectivity, beads, energy and
    // interaction constants), the stored result is returned. The line searches
    // ask for the energy of the same state up to three times per step.
    // Mode: 0 off, 1 on, 2 verify (recompute on every hit, abort on a mismatch).
    int energy_memo_mode = 1;
    long energy_memo_hits = 0;

    // Constructors
    E_Handler() {};
    E_Handler(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo);
    E_Handler(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, std::vector<std::string> inputEnergyNames, std::vector<std::vector<double>> inputEnergyConstants);

    void Add_Bead(Bead *bead);



    virtual void update_face_reference();
    virtual void update_vertex_normals();

    virtual double E_Volume_constraint(const std::vector<double> &Constants) const;
    virtual double E_Area_constraint(const std::vector<double> &Constants) const;
    virtual double E_SurfaceTension(const std::vector<double> &Constants) const;

    virtual double E_MembraneTension(const std::vector<double> &Constants) const;
    virtual double E_ExcessTension(const std::vector<double> &Constants) const;
    virtual double E_Bending(const std::vector<double> &Constants) const;
    virtual double E_Bending_tan(const std::vector<double> &Constants) const;
    // virtual double E_Bending_2(const std::vector<double> &Constants) const;
    virtual double E_Laplace(const std::vector<double> &Constants) const;
    virtual double E_Edge_reg(const std::vector<double> &Constants) const;
    virtual double E_Face_reg(const std::vector<double> &Constants) const;

    virtual VertexData<double> Ev_Bending(const std::vector<double> &Constants) const;
    virtual double V_Bending(const std::vector<double> &Constants, Vertex v) const;

    virtual VertexData<double> Ev_Bending_tan(const std::vector<double> &Constants) const;
    virtual double V_Bending_tan(const std::vector<double> &Constants, Vertex v) const;
    virtual double E_Coverage(const std::vector<double> &Constants) const;

    // We do the precomputed ones now
    // virtual double E_Area_constraint_precomp(const std::vector<double> &Constants) const;
    // virtual double E_SurfaceTension_precomp(const std::vector<double> &Constants) const;
    // virtual double E_Bending_precomp(const std::vector<double> &Constants) const;
    // virtual double E_Edge_reg_precomp(const std::vector<double> &Constants) const;
    virtual void Debug_Coverage(const std::vector<double> &Constants) const;

    virtual VertexData<Vector3> F_Coverage(const std::vector<double> &Constants);
    virtual VertexData<Vector3> F_Volume_constraint(const std::vector<double> &Constants) const;
    virtual VertexData<Vector3> F_SurfaceTension(const std::vector<double> &Constants) const;
    virtual VertexData<Vector3> F_MembraneTension(const std::vector<double> &Constants) const;
    virtual VertexData<Vector3> F_ExcessTension(const std::vector<double> &Constants) const;
    virtual VertexData<Vector3> F_Area_constraint(const std::vector<double> &Constants) const;
    virtual VertexData<Vector3> F_Bending(const std::vector<double> &Constants) const;
    virtual VertexData<Vector3> F_Bending_tan(const std::vector<double> &Constants) const;
    // virtual VertexData<Vector3> F_Bending_2(const std::vector<double> &Constants) const;
    virtual VertexData<Vector3> F_Volume_constraint_2(const std::vector<double> &Constants) const;
    virtual VertexData<Vector3> F_Volume(const std::vector<double> &Constants) const;
    virtual VertexData<Vector3> F_Laplace(const std::vector<double> &Constants) const;
    virtual VertexData<Vector3> F_Edge_reg(const std::vector<double> &Constants) const;
    virtual VertexData<Vector3> F_Edge_reg_2(const std::vector<double> &Constants) const;
    virtual VertexData<Vector3> F_Face_reg(const std::vector<double> &Constants) const;

    // virtual VertexData<Vector3> F_Volume_constraint_precomp(const std::vector<double> &Constants) const;
    // virtual VertexData<Vector3> F_Bending_precomp(const std::vector<double> &Constants) const;
    // virtual VertexData<Vector3> F_Edge_reg_precomp(const std::vector<double> &Constants) const;

    // virtual VertexData<Vector3> F_SurfaceTension_precomp(const std::vector<double> &Constants) const;

    virtual SparseMatrix<double> H_SurfaceTension(const std::vector<double> &Constants);
    virtual SparseMatrix<double> H_Bending(const std::vector<double> &Constants);
    SparseMatrix<double> H_Bending_impl(const std::vector<double> &Constants, bool tan_version);
    virtual SparseMatrix<double> H_Bending_tan(const std::vector<double> &Constants);
    virtual SparseMatrix<double> H_Volume(const std::vector<double> &Constants);

    virtual SparseMatrix<double> H_Laplace(const std::vector<double> &Constants);
    virtual SparseMatrix<double> H_Edge_reg(const std::vector<double> &Constants);
    virtual SparseMatrix<double> H_Face_reg(const std::vector<double> &Constants);

    virtual void Calculate_energies(double *E);
    virtual void Calculate_Lag_norm(double *Norm);
    virtual void Calculate_Lag_norm_Normal(double *Norm);

    // THis function saves the value of the gradient to current gradient but before saves current gradient to previous gradient
    virtual void Calculate_gradient();
    virtual void Calculate_Jacobian();
    virtual void Calculate_Jacobian_Normal();

    // virtual void Calculate_energies_precomp(double *E);
    // virtual void Calculate_gradient_precomp();
    // virtual void Calculate_Jacobian_precomp();

    virtual SparseMatrix<double> Calculate_Hessian();
    virtual SparseMatrix<double> Calculate_Hessian_Normal();
    virtual SparseMatrix<double> Calculate_Hessian_Normal_clipped();
    virtual SparseMatrix<double> Calculate_Hessian_E();
    virtual SparseMatrix<double> Calculate_Hessian_Constraints();
    virtual SparseMatrix<double> Calculate_Hessian_E_Normal();

private:
    void compute_energies(double *E);
    // Everything the energy terms read, as bytes (see energy_memo_mode)
    void energy_state_key(std::vector<char> &key) const;
    std::vector<char> memo_key, memo_scratch;
    bool memo_valid = false;
    double memo_E = 0.0;
    std::vector<double> memo_values;
};