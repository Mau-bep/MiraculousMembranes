#pragma once

// Convergence tests of the L-BFGS integrators ("stopping" in the input file).
//
// Steps are grouped in consecutive windows of `window` steps. At the end of
// a window two quantities are checked:
//
//   energy : |E_end - E_start| < tol_E * max(|E|, E_floor), the net change
//            over the window: from the energy before its first integrator
//            step to the energy after its last one, remeshing and moved beads
//            included. With remeshing on, the integrator keeps removing what
//            each remesh adds back; only the net change shows when BFGS has
//            nothing left to gain (cylinder -> sphere: the integrator-only
//            change stayed at 2.5e-3 per 100 steps while the net change fell
//            to ~5e-5). BFGS-Normal does not remesh, so there the two agree.
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
    bool normal_stop = true; // end the run once both BFGS-Normal tests hold
    int window = 100;
    int patience = 3;
    double E_floor = 1.0;       // energy scale used when |E| is smaller
    double bfgs_tol_E = 1e-3;   // BFGS: relative energy change per window
    double normal_tol_E = 1e-6; // BFGS-Normal: relative energy change per window
    double normal_tol_g = 0.5;  // BFGS-Normal: RMS normal force density in units of KB / R^3
    // Line search steps below 1e-10 this many times in a row reset the L-BFGS
    // history; twice as many switch BFGS to BFGS-Normal or end a BFGS-Normal run.
    // 0 (default) turns this off: runs with beads (the bead displacement cap)
    // take long stretches of such steps and still make progress.
    int stall_steps = 0;
};

class ConvergenceMonitor
{
public:
    explicit ConvergenceMonitor(const StoppingParams &params) : p(params) {}

    // Start over, e.g. after a change of integrator
    void reset()
    {
        n_steps = 0;
        passed = 0;
    }

    // Account one integrator step with the energy before and after it; true
    // when it closes a window
    bool add_step(double E_before, double E_after)
    {
        if (n_steps == 0)
            E_start = E_before;
        E_end = E_after;
        return ++n_steps >= p.window;
    }

    // |E_end - E_start| / max(|E_end|, E_floor) of the window being filled
    double window_dE_rel() const;

    // Record whether the window that just closed passed and start the next
    // one. True once `patience` windows in a row have passed.
    bool close_window(bool pass)
    {
        n_steps = 0;
        passed = pass ? passed + 1 : 0;
        return passed >= p.patience;
    }

    int passed_windows() const { return passed; }
    const StoppingParams &params() const { return p; }

private:
    StoppingParams p;
    int n_steps = 0;
    double E_start = 0.0; // before the first step of the window
    double E_end = 0.0;   // after the last step so far
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
