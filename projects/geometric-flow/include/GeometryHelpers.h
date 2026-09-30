#pragma once

#include <Eigen/Core>
#include "geometrycentral/surface/vertex_position_geometry.h"

// Stack the positions of the first n vertices of a stencil into
// (x0, y0, z0, x1, y1, z1, ...), the layout the Eigen gradient/Hessian kit
// in core/src/geometry.cpp works with. verts can hold Vertex handles or
// vertex indices.
template <class Container>
Eigen::VectorXd stack_positions(const geometrycentral::surface::VertexPositionGeometry &geometry,
                                const Container &verts, size_t n)
{
    Eigen::VectorXd out(3 * n);
    for (size_t i = 0; i < n; i++)
    {
        const geometrycentral::Vector3 &p = geometry.inputVertexPositions[verts[i]];
        out(3 * i) = p.x;
        out(3 * i + 1) = p.y;
        out(3 * i + 2) = p.z;
    }
    return out;
}
