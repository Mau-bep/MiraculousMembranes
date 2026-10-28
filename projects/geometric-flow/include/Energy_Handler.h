#pragma once

#include <cstdint>
#include <Eigen/Core>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "Beads.h"
// #include "Mem-3dg.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// Face normals and areas at the current positions, and the total area summed
// as VertexPositionGeometry::totalArea() does (same order, nan faces skipped)
struct FaceGeometry
{
    std::vector<Vector3> normal;
    std::vector<double> area;
    double total_area = 0.0;
    void compute(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geometry);
};

// The per-element quantities of the bending terms (see compute() in
// Energy_Handler.cpp): per vertex the dual area and the scalar mean curvature
// (or its Tan version), per edge the dihedral angle, 2 tan(dihedral / 2) and the
// length. The vectors are indexed by element index and reused between calls.
struct BendingGeometry
{
    std::vector<double> dual_area;
    std::vector<double> mean_curv;
    std::vector<double> dihedral; // 0 on the boundary
    std::vector<double> tan_term; // only filled by the Tan version
    std::vector<double> length;
    void compute(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geometry, const FaceGeometry &faces,
                 bool tan_version);
};

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

    // The force of energy term i (Energies[i]) at the current positions, as Calculate_gradient adds it to
    // Current_grad. Returns false, leaving `force` alone, for a term with no force: a constraint switched off
    // by a zero constant, or a name that is not known. bead_count is the index of the next "Bead" term and
    // goes up by one for each of them.
    bool Term_force(size_t i, size_t &bead_count, VertexData<Vector3> &force);
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
    // Inside Calculate_energies and Calculate_gradient the terms share one
    // FaceGeometry of the current positions (computed on first use); anywhere
    // else faces() computes it afresh
    friend struct FaceCacheScope;
    mutable FaceGeometry face_cache;
    mutable bool face_cache_fresh = false;
    bool face_cache_active = false;
    const FaceGeometry &faces() const;
    // The same for the bending quantities, [0] plain and [1] Tan version
    mutable BendingGeometry bending_cache[2];
    mutable bool bending_fresh[2] = {false, false};
    const BendingGeometry &bending(bool tan_version) const;
    // Both caches stay valid from one evaluation to the next while the mesh and
    // the positions are bit for bit the same: the gradient at an accepted line
    // search point reuses those of the energy evaluated there. Called when an
    // evaluation starts; drops the caches if anything changed.
    void track_geometry();
    const void *geo_mesh = nullptr, *geo_geometry = nullptr;
    uint64_t geo_counts[5] = {0, 0, 0, 0, 0};
    std::vector<Vector3> geo_positions;
    double total_area() const { return faces().total_area; }

    void compute_energies(double *E);
    // Everything the energy terms read, as bytes (see energy_memo_mode)
    void energy_state_key(std::vector<char> &key) const;
    std::vector<char> memo_key, memo_scratch;
    bool memo_valid = false;
    double memo_E = 0.0;
    std::vector<double> memo_values;
};