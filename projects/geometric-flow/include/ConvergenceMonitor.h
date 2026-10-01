#pragma once

// Convergence tests of the L-BFGS integrators ("stopping" in the input file).
//
// Steps are grouped in consecutive windows of `window` steps. At the end of
// a window two quantities are checked:
//
//   energy : |sum of dE over the window| < tol_E * max(|E|, E_floor)
//            with dE the change made by the integrator step alone
//            (Mem3DG::E_step_end - E_step_start), so the jumps from remeshing
//            and from moving manual beads are not counted.
//   force  : sqrt( 1/N sum_i ((F_i . n_i) / A_i)^2 ) < tol_g * KB / R^3
//            the RMS normal force density, with n_i the normals the
//            BFGS-Normal step moves along, A_i a third of the adjacent face
//            areas, R = sqrt(A / 4pi) and KB the bending modulus (1 when no
//            bending energy is on). Dividing by A_i makes it a pressure that
//            does not depend on the mesh resolution; KB / R^3 is the pressure
//            scale of a bent vesicle of that size.
//
// A test holds once `patience` windows in a row pass it.
//
// The functions below only read positions and forces, so evaluating them does
// not change the simulation.

#include <string>
#include <vector>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

struct StoppingParams
{
    bool bfgs_switch = true; // BFGS -> BFGS-Normal once the energy test holds
    int window = 100;
    int patience = 3;
    double E_floor = 1.0;       // energy scale used when |E| is smaller
    double bfgs_tol_E = 1e-5;   // BFGS: relative energy change per window
    double normal_tol_E = 1e-7; // BFGS-Normal: relative energy change per window
    double normal_tol_g = 0.1;  // BFGS-Normal: RMS normal force density in units of KB / R^3
};

class ConvergenceMonitor
{
public:
    explicit ConvergenceMonitor(const StoppingParams &params) : p(params) {}

    // Start over, e.g. after a change of integrator
    void reset()
    {
        n_steps = 0;
        sum_dE = 0.0;
        passed = 0;
    }

    // Account one integrator step; true when it closes a window
    bool add_step(double dE)
    {
        sum_dE += dE;
        return ++n_steps >= p.window;
    }

    // |sum dE| / max(|E|, E_floor) of the window being filled
    double window_dE_rel(double E) const;

    // Record whether the window that just closed passed and start the next
    // one. True once `patience` windows in a row have passed.
    bool close_window(bool pass)
    {
        n_steps = 0;
        sum_dE = 0.0;
        passed = pass ? passed + 1 : 0;
        return passed >= p.patience;
    }

    int passed_windows() const { return passed; }
    const StoppingParams &params() const { return p; }

private:
    StoppingParams p;
    int n_steps = 0;
    double sum_dE = 0.0;
    int passed = 0;
};

// sqrt( 1/N sum_i ((force_i . normal_i) / A_i)^2 ), A_i = a third of the
// areas of the faces around vertex i. Without normals |force_i| is used.
double force_density_rms(ManifoldSurfaceMesh &mesh, const VertexPositionGeometry &geometry,
                         const VertexData<Vector3> &force, const VertexData<Vector3> *normals);

// KB / R^3 with R = sqrt(A / 4pi); KB = 1 when bending_modulus is 0
double force_density_scale(ManifoldSurfaceMesh &mesh, const VertexPositionGeometry &geometry, double bending_modulus);

// First constant of the first bending energy that is on (the energy handler
// skips them below 1e-5), 0 if there is none
double bending_modulus(const std::vector<std::string> &energies, const std::vector<std::vector<double>> &constants);
