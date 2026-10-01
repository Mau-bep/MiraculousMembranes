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
     * delta = 0.25 * sigma
     *
     * weight is nonzero only for:
     *
     *   sigma - delta < r < sigma + delta
     *
     * The function returns weight(r), and writes dweight/dr in dWeightDr.
     */
    double coverageShellWeight(double r, double sigma, double &dWeightDr);

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
        double sigma);

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