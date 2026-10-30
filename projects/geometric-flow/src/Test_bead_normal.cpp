// Checks for the "Normal" constraint of a bead, the bead of a planar membrane that only moves along the normal of the plane:
//   1. Moves_along_axis is true for a "default" bead with the constraint while Axis_active (BFGS-Normal) only;
//      outside it the bead moves freely
//   2. Free_force keeps the component of the force along the axis; unconstrained, "Radial", "froze" and "rigid" beads are untouched
//   3. Move_bead (both overloads) moves the bead by dt (F.n) n and nothing else, a frozen bead stays
//
// Usage: Test_bead_normal

#include <cmath>
#include <cstdio>

#include "Beads.h"

namespace
{
    bool check(bool ok, const char *what)
    {
        std::printf("%s: %s\n", ok ? "ok  " : "FAIL", what);
        return ok;
    }

    bool same(Vector3 a, Vector3 b)
    {
        return (a - b).norm() < 1e-14;
    }
} // namespace

int main()
{
    bool ok = true;
    const Vector3 axis{0.6, 0.0, 0.8};
    const Vector3 F{1.0, 2.0, 3.0};
    const Vector3 kept = 3.0 * axis; // F.axis = 0.6 + 2.4

    Bead bead;
    bead.state = "default";
    bead.Constraint = "Normal";
    bead.Constraint_axis = axis;

    // Outside BFGS-Normal the constraint is off: the bead moves like an unconstrained one
    ok &= check(!bead.Moves_along_axis() && same(bead.Free_force(F), F), "inactive: Free_force keeps the whole force");
    bead.Pos = Vector3{1.0, 1.0, 1.0};
    bead.Total_force = F;
    bead.Move_bead(0.5, Vector3{0.0, 0.0, 0.0}, F);
    const bool free_3d = same(bead.Pos, Vector3{1.5, 2.0, 2.5});
    bead.Pos = Vector3{1.0, 1.0, 1.0};
    bead.Move_bead(0.5, Vector3{0.0, 0.0, 0.0});
    ok &= check(free_3d && same(bead.Pos, Vector3{1.5, 2.0, 2.5}), "inactive: both Move_bead overloads move along the whole force");

    bead.Axis_active = true;
    ok &= check(bead.Moves_along_axis(), "active: a default bead with the Normal constraint moves along its axis");
    ok &= check(same(bead.Free_force(F), kept), "Free_force keeps F.n n");
    ok &= check(same(bead.Free_force(axis), axis) && same(bead.Free_force(Vector3{0.0, 1.0, 0.0}), Vector3{0.0, 0.0, 0.0}),
                "Free_force: the axis passes, the perpendicular direction is removed");

    bead.Pos = Vector3{1.0, 1.0, 1.0};
    bead.Move_bead(0.5, Vector3{0.0, 0.0, 0.0}, F);
    ok &= check(same(bead.Pos, Vector3{1.0, 1.0, 1.0} + 0.5 * kept), "Move_bead(dt, center, F) moves by dt F.n n");

    bead.Pos = Vector3{1.0, 1.0, 1.0};
    bead.Total_force = F;
    bead.Move_bead(0.25, Vector3{0.0, 0.0, 0.0});
    ok &= check(same(bead.Pos, Vector3{1.0, 1.0, 1.0} + 0.25 * kept), "Move_bead(dt, center) moves by dt F.n n with Total_force");

    bead.Pos = Vector3{1.0, 1.0, 1.0};
    bead.Move_bead(0.5, Vector3{0.1, 0.2, 0.3}, F);
    ok &= check(same(bead.Pos, Vector3{1.0, 1.0, 1.0} + 0.5 * kept - Vector3{0.1, 0.2, 0.3}), "the recentering shift is subtracted as for the other beads");

    bead.Pos = Vector3{1.0, 1.0, 1.0};
    bead.Move_bead(1.0, Vector3{0.0, 0.0, 0.0}, Vector3{0.0, 5.0, 0.0});
    ok &= check(same(bead.Pos, Vector3{1.0, 1.0, 1.0}), "a force perpendicular to the axis does not move the bead");

    Bead frozen = bead;
    frozen.state = "froze";
    frozen.Pos = Vector3{1.0, 1.0, 1.0};
    frozen.Move_bead(1.0, Vector3{0.0, 0.0, 0.0}, F);
    ok &= check(!frozen.Moves_along_axis() && same(frozen.Free_force(F), F) && same(frozen.Pos, Vector3{1.0, 1.0, 1.0}),
                "a frozen bead is not a variable and does not move");

    Bead free_bead;
    free_bead.state = "default";
    free_bead.Pos = Vector3{1.0, 1.0, 1.0};
    ok &= check(!free_bead.Moves_along_axis() && same(free_bead.Free_force(F), F), "an unconstrained bead keeps its whole force");
    free_bead.Move_bead(0.5, Vector3{0.0, 0.0, 0.0}, F);
    ok &= check(same(free_bead.Pos, Vector3{1.5, 2.0, 2.5}), "an unconstrained bead moves along the whole force");

    Bead radial = free_bead;
    radial.Constraint = "Radial";
    radial.Constraint_constants = {0.0};
    ok &= check(!radial.Moves_along_axis() && same(radial.Free_force(F), F), "a Radial bead is not affected by the Normal constraint");

    std::printf("%s\n", ok ? "All checks passed" : "SOME CHECKS FAILED");
    return ok ? 0 : 1;
}
