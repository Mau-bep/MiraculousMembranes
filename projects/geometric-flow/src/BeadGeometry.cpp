#include "BeadGeometry.h"

#include <cmath>
#include <stdexcept>

namespace bead_geometry
{

    double coverageShellWeight(double r, double sigma, double &dWeightDr)
    {
        if (sigma <= 0.0)
        {
            throw std::invalid_argument(
                "coverageShellWeight(): sigma must be strictly positive.");
        }

        const double delta = 0.25 * sigma;
        const double x = (r - sigma) / delta;

        dWeightDr = 0.0;

        /*
         * Outside:
         *
         * sigma - delta <= r <= sigma + delta
         */
        if (std::abs(x) >= 1.0)
        {
            return 0.0;
        }

        /*
         * w(r) = 0.5 * [1 + cos(pi * x)]
         *
         * w(sigma) = 1
         * w(sigma - delta) = w(sigma + delta) = 0
         */
        const double weight =
            0.5 * (1.0 + std::cos(PI_VALUE * x));

        /*
         * dw/dr = -pi / (2 delta) sin(pi x)
         */
        dWeightDr =
            -0.5 * PI_VALUE * std::sin(PI_VALUE * x) / delta;

        return weight;
    }

    SolidAngleResult triangleSolidAngleGradient(
        const Vector3 &p0,
        const Vector3 &p1,
        const Vector3 &p2,
        const Vector3 &beadPosition)
    {
        SolidAngleResult result;

        const Vector3 a = p0 - beadPosition;
        const Vector3 b = p1 - beadPosition;
        const Vector3 c = p2 - beadPosition;

        const double ra = norm(a);
        const double rb = norm(b);
        const double rc = norm(c);

        if (ra < EPS || rb < EPS || rc < EPS)
        {
            return result;
        }

        const double ab = dot(a, b);
        const double bc = dot(b, c);
        const double ca = dot(c, a);

        /*
         * N = det(a,b,c) = a . (b x c)
         */
        const double N = dot(a, cross(b, c));

        /*
         * D = |a||b||c| + (a.b)|c| + (b.c)|a| + (c.a)|b|
         */
        const double D =
            ra * rb * rc + ab * rc + bc * ra + ca * rb;

        const double denominator = D * D + N * N;

        if (denominator < EPS)
        {
            return result;
        }

        result.omegaSigned = 2.0 * std::atan2(N, D);

        /*
         * Derivatives of N.
         */
        const Vector3 gradN_a = cross(b, c);
        const Vector3 gradN_b = cross(c, a);
        const Vector3 gradN_c = cross(a, b);

        /*
         * Derivatives of D.
         */
        const Vector3 gradD_a =
            (rb * rc + bc) * a / ra + rc * b + rb * c;

        const Vector3 gradD_b =
            (ra * rc + ca) * b / rb + rc * a + ra * c;

        const Vector3 gradD_c =
            (ra * rb + ab) * c / rc + ra * b + rb * a;

        /*
         * omega = 2 atan2(N,D)
         *
         * d omega = 2 [D dN - N dD] / [D^2 + N^2]
         */
        const double prefactor = 2.0 / denominator;

        result.grad0 =
            prefactor * (D * gradN_a - N * gradD_a);

        result.grad1 =
            prefactor * (D * gradN_b - N * gradD_b);

        result.grad2 =
            prefactor * (D * gradN_c - N * gradD_c);

        result.valid = true;

        return result;
    }

    FaceCoverageData evaluateFaceCoverage(
        const VertexPositionGeometry &geometry,
        Face f,
        const Vector3 &beadPosition,
        double sigma)
    {
        FaceCoverageData result;

        Halfedge he = f.halfedge();

        const Vertex v0 = he.vertex();
        const Vertex v1 = he.next().vertex();
        const Vertex v2 = he.next().next().vertex();

        const Vector3 p0 = geometry.inputVertexPositions[v0];
        const Vector3 p1 = geometry.inputVertexPositions[v1];
        const Vector3 p2 = geometry.inputVertexPositions[v2];

        const Vector3 centroid = (p0 + p1 + p2) / 3.0;

        const Vector3 beadToCentroid =
            centroid - beadPosition;

        const double centroidDistance = norm(beadToCentroid);

        if (centroidDistance < EPS)
        {
            return result;
        }

        result.centroidDistance = centroidDistance;
        result.radialDirection =
            beadToCentroid / centroidDistance;

        /*
         * Radial cosine shell selection.
         */
        result.weight = coverageShellWeight(
            centroidDistance,
            sigma,
            result.dWeightDr);

        if (result.weight <= 0.0)
        {
            return result;
        }

        /*
         * Outward mesh normal.
         */
        Vector3 faceNormal = geometry.faceNormal(f);
        const double normalLength = norm(faceNormal);

        if (normalLength < EPS)
        {
            return result;
        }

        faceNormal /= normalLength;

        /*
         * External-bead convention:
         *
         * bead -> centroid must be antiparallel to outward normal.
         *
         * Keep only:
         *
         * dot(radialDirection, faceNormal) < 0
         *
         * This criterion is deliberately selection-only. Its derivative
         * is not included.
         */
        if (dot(result.radialDirection, faceNormal) >= 0.0)
        {
            return result;
        }

        const SolidAngleResult solidAngle =
            triangleSolidAngleGradient(
                p0, p1, p2, beadPosition);

        if (!solidAngle.valid)
        {
            return result;
        }

        result.omegaSigned = solidAngle.omegaSigned;
        result.omegaUnsigned = std::abs(solidAngle.omegaSigned);

        /*
         * d|omega|/dx = sign(omega) d omega/dx
         *
         * The omega == 0 case is non-differentiable, but it should not
         * normally occur for a selected nondegenerate face.
         */
        const double omegaSign =
            (solidAngle.omegaSigned >= 0.0) ? 1.0 : -1.0;

        result.gradOmega0 = omegaSign * solidAngle.grad0;
        result.gradOmega1 = omegaSign * solidAngle.grad1;
        result.gradOmega2 = omegaSign * solidAngle.grad2;

        result.selected = true;

        return result;
    }

    void validateCoverageConstants(const std::vector<double> &Constants, size_t nBeads)
    {
        if (Constants.size() != nBeads + 1)
        {
            throw std::invalid_argument(
                "Coverage: Constants must have exactly Beads.size() + 1 values: "
                "[K, cov_0, cov_1, ..., cov_(N-1)].");
        }
    }

    bool validCoverageTarget(double cov)
    {
        return cov == -1.0 || (cov >= 0.0 && cov <= 1.0);
    }

} // namespace bead_geometry