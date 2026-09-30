#pragma once

#include <Eigen/Core>
#include "geometrycentral/surface/vertex_position_geometry.h"

// Stack the positions of the N vertices of a stencil into
// (x0, y0, z0, x1, y1, z1, ...), the layout the Eigen gradient/Hessian kit
// in core/src/geometry.cpp works with. verts can hold Vertex handles or
// vertex indices. Fixed size, so no allocation in the energy loops.
template <int N, class Container>
Eigen::Matrix<double, 3 * N, 1> stack_positions(const geometrycentral::surface::VertexPositionGeometry &geometry,
                                                const Container &verts)
{
    Eigen::Matrix<double, 3 * N, 1> out;
    for (int i = 0; i < N; i++)
    {
        const geometrycentral::Vector3 &p = geometry.inputVertexPositions[verts[i]];
        out(3 * i) = p.x;
        out(3 * i + 1) = p.y;
        out(3 * i + 2) = p.z;
    }
    return out;
}
