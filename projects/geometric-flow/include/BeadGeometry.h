#pragma once

#include <cmath>
#include <string>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace bead_geometry
{

    using namespace geometrycentral;
    using namespace geometrycentral::surface;

    constexpr double EPS = 1e-12;
    constexpr double PI_VALUE = 3.141592653589793238462643383279502884;

    /*
     * The signed solid angle and its derivatives with respect to the
     * three triangle vertices.
     *
     * omegaSigned is:
     *
     *   2 atan2(N, D)
     *
     * where:
     *
     *   N = a . (b x c)
     *   D = |a||b||c| + (a.b)|c| + (b.c)|a| + (c.a)|b|
     *
     * with:
     *
     *   a = p0 - beadPosition
     *   b = p1 - beadPosition
     *   c = p2 - beadPosition
     */
    struct SolidAngleResult
    {
        double omegaSigned = 0.0;

        Vector3 grad0 = Vector3({0.0, 0.0, 0.0});
        Vector3 grad1 = Vector3({0.0, 0.0, 0.0});
        Vector3 grad2 = Vector3({0.0, 0.0, 0.0});

        bool valid = false;
    };

    /*
     * The complete contribution of one face to bead coverage.
     *
     * `selected == true` only if:
     *
     * 1. The face centroid lies in the shell (cosine, or the support of the profile):
     *
     *      sigma - 0.25 sigma < distance < sigma + 0.25 sigma
     *
     * 2. The outward face normal points toward the bead:
     *
     *      dot(bead -> centroid, outwardNormal) < 0
     *
     * The facing selection is intentionally not differentiated.
     */
    struct FaceCoverageData
    {
        bool selected = false;

        double centroidDistance = 0.0;

        // Unit vector from bead center toward the face centroid.
        Vector3 radialDirection = Vector3({0.0, 0.0, 0.0});

        // Cosine shell weight and radial derivative.
        double weight = 0.0;
        double dWeightDr = 0.0;

        double omegaSigned = 0.0;
        double omegaUnsigned = 0.0;

        /*
         * Derivatives of abs(omega) with respect to p0, p1, and p2.
         *
         * That is:
         *
         * gradOmega0 = d |omega| / d p0
         */
        Vector3 gradOmega0 = Vector3({0.0, 0.0, 0.0});
        Vector3 gradOmega1 = Vector3({0.0, 0.0, 0.0});
        Vector3 gradOmega2 = Vector3({0.0, 0.0, 0.0});
    };

    /*
     * Cosine shell centered at the physical bead radius sigma.
     *
     * delta = shellWidth * sigma  (shellWidth = 0.25 unless given)
     *
     * weight is nonzero only for:
     *
     *   sigma - delta < r < sigma + delta
     *
     * The function returns weight(r), and writes dweight/dr in dWeightDr.
     * shellWidth must be in (0, 1).
     *
     * shellPower p >= 1 (default 1, the plain shell) raises the shell to the
     * power p: weight = [ 1/2 (1 + cos(pi x)) ]^p = cos^(2p)(pi x / 2),
     * x = (r - sigma) / delta. The maximum (1 at r = sigma) and the support
     * are unchanged, the shell gets narrower around its maximum (p = 2, 4:
     * cos^4, cos^8 of pi x / 2). p = 1 takes exactly the old code path.
     * p < 1 is rejected: the slope would diverge at the edges of the shell
     * (weight ~ (1 - |x|)^(2p)), so the force would no longer be continuous.
     */
    double coverageShellWeight(double r, double sigma, double &dWeightDr, double shellWidth = 0.25,
                               double shellPower = 1.0);

    /*
     * Radial profile of the adhesion weight w(r). The default is the cosine shell above
     * (shellPower included). The others are radial pair potentials V(r) turned into a
     * weight, w = -V_shifted / depth, so that w = 1 at r = sigma (the bead radius, the
     * minimum of the well) for every profile and W = inter_str stays the adhesion energy
     * per unit of contact area, the same as the "w" of Deserno's phase diagram:
     *
     *   Frenkel:  V = alpha (s^2/r^2 - 1) (rc^2/r^2 - 1)^2 (zero value and slope at rc, depth 1),
     *             s solved so that the minimum is at r = sigma: s^2 = 2 rc^2 sigma^2 / (3 rc^2 - sigma^2)
     *   LJ:       V = (sigma/r)^12 - 2 (sigma/r)^6, minimum -1 at r = sigma
     *   Morse:    V = (1 - exp(-alpha (r - sigma) / sigma))^2 - 1, minimum -1 at r = sigma
     *
     * All three keep their repulsive core (w < 0 for r below about sigma) and are cut at
     * rc = sigma (1 + width). LJ and Morse are shifted by V(rc) so that the energy is zero
     * at the cutoff, w = (V(rc) - V(r)) / (1 + V(rc)); the force then jumps at rc (Frenkel
     * does not). A narrow width makes the shifted depth 1 + V(rc) small and the wall
     * relatively stiff, since the depth is renormalised to 1.
     *
     * `width` is the shell width of the profile (the cutoff), in (0, MAX_PROFILE_WIDTH];
     * `power` is only used by Cosine and `alpha` only by Morse.
     */
    enum class ShellKind
    {
        Cosine,
        Frenkel,
        LJ,
        Morse
    };

    constexpr double MAX_PROFILE_WIDTH = 4.0;

    struct ShellProfile
    {
        ShellKind kind = ShellKind::Cosine;
        double width = 0.25;
        double power = 1.0;
        double alpha = 6.0;
    };

    // Name of the profile as written in the input file ("cosine", "frenkel", "lj", "morse")
    const char *shellKindName(ShellKind kind);
    // The inverse; throws std::invalid_argument for an unknown name
    ShellKind shellKindFromName(const std::string &name);

    /*
     * Weight of the profile at distance r from the bead centre (radius sigma) and dw/dr.
     * Cosine takes exactly the coverageShellWeight path. Throws std::invalid_argument for
     * a non-positive sigma, a width outside (0, MAX_PROFILE_WIDTH] (Cosine: (0, 1)),
     * a Morse alpha <= 0 or a Cosine power < 1.
     */
    double shellProfileWeight(double r, double sigma, double &dWeightDr, const ShellProfile &profile);

    /*
     * Angular patch of the bead (a pinned contact line). Only the part of the
     * bead within the polar angle theta_Wr of an axis can adhere, as in the
     * adhesive patch of Saric & Cacciuto (theta_Wr of van der Wel et al., Sci.
     * Rep. 6, 32825): the membrane wraps the patch and its edge pins the
     * contact line, whatever the shape of the free membrane beyond it.
     *
     * With u the unit vector from the bead centre to the face centroid and
     *
     *   t = dot(u, axis),   cosEdge = cos(theta_Wr),
     *   dCos = cos(theta_Wr - width) - cosEdge > 0,
     *
     * the patch weight is 0 for t <= cosEdge, 1 for t >= cosEdge + dCos and
     * 1/2 (1 - cos(pi (t - cosEdge) / dCos)) in between: a C1 ramp of angular
     * width `width` just inside theta_Wr (the energy and its slope are
     * continuous, so the line searches see a smooth function). The axis points
     * from the bead centre towards the pole of the patch (where the membrane
     * touches first), so the inactive cap is opposite to it. The axis is fixed
     * in space (use it with frozen beads or beads that do not turn).
     *
     * Returns the weight and writes d weight / dt in dWeightDt.
     */
    double patchWeight(double t, double cosEdge, double dCos, double &dWeightDt);

    // The same with the second derivative too (a jump at the two ends of the ramp), for the Newton iteration of the free axis
    double patchWeight2(double t, double cosEdge, double dCos, double &dWeightDt, double &d2WeightDt2);

    /*
     * Compute the signed solid angle of triangle p0,p1,p2 as seen from
     * beadPosition, together with derivatives with respect to p0,p1,p2.
     */
    SolidAngleResult triangleSolidAngleGradient(
        const Vector3 &p0,
        const Vector3 &p1,
        const Vector3 &p2,
        const Vector3 &beadPosition);

    /*
     * Evaluate whether face f contributes to external-bead coverage and,
     * if so, calculate its shell weight, unsigned solid angle, and all
     * derivatives needed by Coverage and Adhesion.
     *
     * Assumes mesh normals are outward-oriented.
     *
     * side = +1 (default): the bead is outside, a face counts if its outward
     * normal points to the bead (dot(bead -> centroid, normal) < 0).
     * side = -1: the bead is inside the closed membrane, a face counts if its
     * outward normal points away from the bead (dot(bead -> centroid, normal) > 0).
     *
     * profile = nullptr or Cosine: the cosine shell of shellWidth and shellPower, as always.
     * Another kind replaces the weight by that profile (shellWidth and shellPower are then
     * ignored); its weight can be negative (the repulsive core), a face is selected whenever
     * the weight is not exactly zero.
     */
    FaceCoverageData evaluateFaceCoverage(
        const VertexPositionGeometry &geometry,
        Face f,
        const Vector3 &beadPosition,
        double sigma,
        double shellWidth = 0.25,
        double shellPower = 1.0,
        double side = 1.0,
        const ShellProfile *profile = nullptr);

    // Coverage constants are [K, cov_0, ..., cov_(N-1)], one target per bead.
    // Throws std::invalid_argument if the count does not match.
    void validateCoverageConstants(const std::vector<double> &Constants, size_t nBeads);

    // Valid coverage targets: -1 (bead disabled) or a value in [0, 1].
    bool validCoverageTarget(double cov);

    /*
     * Contact weight of a face at signed height h above a plane, centred on
     * the plane:
     *
     *   w(h) = 0.5 * [1 + cos(pi h / width)]   for |h| < width, else 0
     *
     * 1 at h = 0, 0 with zero slope at h = +-width. Writes dw/dh in dWdh.
     */
    double planeContactWeight(double h, double width, double &dWdh);

    /*
     * Contact weight of a face at height h above a plane for a pair-potential profile
     * (Frenkel, LJ or Morse, see ShellProfile): the planar version of the sphere bead, with
     * the plane in place of the bead surface. sigma is the "radius", the height of the minimum
     * of the well (w = 1 there), profile.width sets the cutoff sigma (1 + width), and the
     * repulsive core (w < 0 for h below about sigma) keeps the membrane on the +n side.
     * The potentials diverge at h = 0, so below h = PLANE_PROFILE_FLOOR sigma the weight is
     * continued linearly (value and slope matched), still a very stiff wall but finite.
     * Writes dw/dh in dWdh.
     */
    constexpr double PLANE_PROFILE_FLOOR = 0.25;

    double planeProfileWeight(double h, double sigma, double &dWdh, const ShellProfile &profile);

    /*
     * Repulsive wall for the faces that look away from the plane (the back of a membrane
     * sheet): the repulsive branch of the same profile, shifted to zero at its minimum,
     *
     *   g(h) = 1 - w(h)   for h < sigma,   0 beyond,
     *
     * so g and its slope vanish at h = sigma (w has its maximum 1 there) and g grows without
     * a bump up to the core, monotonically decreasing in h: a back face is only ever pushed
     * away from the plane. Same linear continuation as planeProfileWeight below the floor.
     * Writes dg/dh in dGdh.
     */
    double planeProfileRepulsion(double h, double sigma, double &dGdh, const ShellProfile &profile);

} // namespace bead_geometry