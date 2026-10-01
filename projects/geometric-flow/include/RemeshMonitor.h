#pragma once

// Read-only mesh quality check for the "quality" adaptive remeshing mode.
//
// Counts the edges the remesher (deps/geometry-central remeshing.cpp) would
// act on if it ran now, using its own tests:
//   long  : split test of splitWorstEdges, L * sqrt((s_a + s_b)/2) > 1 and
//           L > 2 * min_absolute_length
//   short : collapse test of improveFaces on interior edges shorter than
//           size_max: an endpoint of valence 3 or 4, or the metric part of
//           shouldCollapse (every edge around the collapsed vertex stays below
//           0.9 in the sizing metric). The foldover and aspect vetoes are left
//           out; edges they block show up in the baseline after each remesh.
//   flip  : isDelaunay_improv fails (cos a + cos b < -0.05)
// with the sizing field recomputed from the current positions the way
// remesh() does (curvature sizing / refine_angle^2, clamped to
// [1/size_max^2, 1/size_min^2], area averaged to the vertices).
//
// Everything is computed from inputVertexPositions with local buffers: no
// geometry quantity is required or refreshed, so calling it does not change
// the simulation.

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/remeshing.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

struct MeshQuality
{
    size_t n_edges = 0;
    size_t n_long = 0;  // would be split
    size_t n_short = 0; // would be collapsed
    size_t n_flip = 0;  // not Delaunay, would be flipped
    size_t n_bad = 0;   // edges in at least one of the three sets

    double bad_fraction() const { return n_edges ? double(n_bad) / n_edges : 0.0; }
    double fraction(size_t n) const { return n_edges ? double(n) / n_edges : 0.0; }
};

// The vertex sizing remesh() would use for the current positions
VertexData<double> remesh_vertex_sizing(ManifoldSurfaceMesh &mesh, const VertexPositionGeometry &geometry,
                                        const RemeshOptions &options);

MeshQuality measure_mesh_quality(ManifoldSurfaceMesh &mesh, const VertexPositionGeometry &geometry,
                                 const RemeshOptions &options);
