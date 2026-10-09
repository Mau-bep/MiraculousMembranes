#pragma once

#include <cmath>

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
     * 1. The face centroid lies in the cosine shell:
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
     */
    FaceCoverageData evaluateFaceCoverage(
        const VertexPositionGeometry &geometry,
        Face f,
        const Vector3 &beadPosition,
        double sigma,
        double shellWidth = 0.25,
        double shellPower = 1.0);

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

} // namespace bead_geometry