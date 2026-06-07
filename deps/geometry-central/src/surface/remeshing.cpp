
#include "geometrycentral/surface/remeshing.h"

namespace geometrycentral {
namespace surface {

// The default trace options
const RemeshOptions defaultRemeshOptions;


int remesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, RemeshOptions options) {
  MutationManager mm(mesh, geom);
  return remesh(mesh, geom, mm, options);
}

int remesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm, RemeshOptions options) {


  geom.requireFaceSizing();
  for (Face f : mesh.faces()) {
    geom.faceSizing[f] = clamp(geom.faceSizing[f] / (options.refine_angle * options.refine_angle),
                               1.0 / (options.max_absolute_length * options.max_absolute_length),
                               1.0 / (options.min_absolute_length * options.min_absolute_length));
  }

  geom.requireVertexSizing();

  bool doConnectivityChanges = true;

  options.maxIterations = 1;
  for (size_t iIt = 0; iIt < options.maxIterations; iIt++) {
    size_t nFlips = 10;
    if (doConnectivityChanges) {

      splitWorstEdges(mesh, geom, mm, options);
      nFlips = fixDelaunay(mesh, geom, mm);
      options.numberOp += nFlips;
      improveFaces(mesh, geom, mm, options);
    }

    nFlips = fixDelaunay(mesh, geom, mm);
    // std::cout<<"The number of flips is "<< nFlips<<"\n";
    options.numberOp += nFlips;


    geom.inputVertexPositions = geom.vertexPositions;
    mesh.compress();
    geom.refreshQuantities();
    double smoothing = smoothByLaplacian(mesh, geom, mm);
    while (smoothing > 1e-4) smoothing = smoothByLaplacian(mesh, geom, mm);
  }
  geom.unrequireFaceSizing();
  geom.unrequireVertexSizing();
  geom.inputVertexPositions = geom.vertexPositions;
  geom.purgeQuantities();

  return options.numberOp;
}


void remesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, RemeshOptions options,
            std::vector<Face> activeFaces) {

  MutationManager mm(mesh, geom);
  std::vector<Face> activeFacesSplit = splitSubset(activeFaces, mesh, geom, mm, options);
  collapseSubset(activeFacesSplit, mesh, geom, mm, options);
  mesh.compress();
}

void dynamic_remesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, RemeshOptions options) {
  MutationManager mm(mesh, geom);
  return dynamic_remesh(mesh, geom, mm, options);
}

void dynamic_remesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm,
                    RemeshOptions options) {
  if (options.targetEdgeLength < 0) {
    double meanLength = 0;
    geom.requireEdgeLengths();
    for (Edge e : mesh.edges()) {
      meanLength += geom.edgeLengths[e];
    }
    geom.unrequireEdgeLengths();
    meanLength /= mesh.nEdges();

    options.targetEdgeLength = abs(options.targetEdgeLength * meanLength);
  }
  // std::cout << "THis is being called \n";
  bool doConnectivityChanges = true;

  for (size_t iIt = 0; iIt < options.maxIterations; iIt++) {
    if (doConnectivityChanges) {
      doConnectivityChanges = adjustEdgeLengths(mesh, geom, mm, options);
    }

    size_t nFlips = 10;

    fixDelaunay(mesh, geom, mm);
    // std::cout<<"The number of flips is "<< nFlips<<"\n";
    double flowDist = 1;
    // flowDist=1;
    // switch (options.smoothStyle) {
    // case RemeshSmoothStyle::Circumcentric:
    //   flowDist = smoothByCircumcenter(mesh, geom, mm, 1, options.boundaryCondition);
    //   break;
    // case RemeshSmoothStyle::Laplacian:
    //   flowDist = smoothByLaplacian(mesh, geom, mm, 1, options.boundaryCondition);
    //   break;
    // }

    // std::cout << iIt << " : " << changedConnectivity << " " << nFlips << " " << flowDist << std::endl;
    if ((nFlips == 0) && (flowDist < 0.01)) break;
  }
  geom.refreshQuantities();
}


void remesh_smoothing(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, RemeshOptions options) {
  MutationManager mm(mesh, geom);

  // std::cout<<"The number of flips is "<< nFlips<<"\n";
  double flowDist = 1;
  // flowDist=1;
  switch (options.smoothStyle) {
  case RemeshSmoothStyle::Circumcentric:
    flowDist = smoothByCircumcenter(mesh, geom, mm, 1, options.boundaryCondition);
    break;
  case RemeshSmoothStyle::Laplacian:
    flowDist = smoothByLaplacian(mesh, geom, mm, 1, options.boundaryCondition);
    break;
  }

  std::cout << "The smoothing distance is " << flowDist << " ... huh \n ";
  // std::cout << iIt << " : " << changedConnectivity << " " << nFlips << " " << flowDist << std::endl;


  geom.refreshQuantities();
}

Vector3 vertexNormal(VertexPositionGeometry& geom, Vertex v, MutationManager& mm) {
  Vector3 totalNormal = Vector3::zero();
  for (Corner c : v.adjacentCorners()) {
    Vector3 cornerNormal = geom.cornerAngle(c) * geom.faceNormal(c.face());
    totalNormal += cornerNormal;
  }
  return normalize(totalNormal);
}

Vector3 boundaryVertexTangent(VertexPositionGeometry& geom, Vertex v, MutationManager& mm) {
  if (v.isBoundary()) {
    auto edgeVec = [&](Edge e) -> Vector3 {
      return (geom.vertexPositions[e.halfedge().tipVertex()] - geom.vertexPositions[e.halfedge().tailVertex()])
          .normalize();
    };

    Vector3 totalTangent = Vector3::zero();
    for (Edge e : v.adjacentEdges()) {
      if (e.isBoundary()) {
        totalTangent += edgeVec(e);
      }
    }
    return totalTangent.normalize();
  } else {
    return Vector3::zero();
  }
}

inline Vector3 projectToPlane(Vector3 v, Vector3 norm) { return v - norm * dot(norm, v); }
inline Vector3 projectToLine(Vector3 v, Vector3 tangent) { return tangent * dot(tangent, v); }

inline Vector3 edgeMidpoint(SurfaceMesh& mesh, VertexPositionGeometry& geom, Edge e) {
  Vector3 endPos1 = geom.vertexPositions[e.halfedge().tailVertex()];
  Vector3 endPos2 = geom.vertexPositions[e.halfedge().tipVertex()];
  return (endPos1 + endPos2) / 2;
}
inline Vector3 edgeButterfly(SurfaceMesh& mesh, VertexPositionGeometry& geom, Edge e) {
  Halfedge he = e.halfedge();
  Vector3 P1 = geom.vertexPositions[he.vertex()];
  Vector3 P2 = geom.vertexPositions[he.next().vertex()];
  Vector3 Q1 = geom.vertexPositions[he.next().next().vertex()];
  Vector3 Q2 = geom.vertexPositions[he.twin().next().next().vertex()];
  Vector3 R1 = geom.vertexPositions[he.next().next().twin().next().next().vertex()];
  Vector3 R2 = geom.vertexPositions[he.next().twin().next().next().vertex()];
  Vector3 R3 = geom.vertexPositions[he.twin().next().twin().next().next().vertex()];
  Vector3 R4 = geom.vertexPositions[he.twin().next().next().twin().next().next().vertex()];

  return (8 * (P1 + P2) + 2 * (Q1 + Q2) - (R1 + R2 + R3 + R4)) / 16.0;
}

Vector3 findCircumcenter(Vector3 p1, Vector3 p2, Vector3 p3) {
  // barycentric coordinates of circumcenter
  double a = (p3 - p2).norm();
  double b = (p3 - p1).norm();
  double c = (p2 - p1).norm();
  double a2 = a * a;
  double b2 = b * b;
  double c2 = c * c;
  Vector3 circumcenterLoc{a2 * (b2 + c2 - a2), b2 * (c2 + a2 - b2), c2 * (a2 + b2 - c2)};
  // normalize to sum of 1
  circumcenterLoc = normalizeBarycentric(circumcenterLoc);

  // change back to space
  return circumcenterLoc[0] * p1 + circumcenterLoc[1] * p2 + circumcenterLoc[2] * p3;
}

Vector3 findCircumcenter(VertexPositionGeometry& geom, Face f) {
  // retrieve the face's vertices
  int index = 0;
  Vector3 p[3];
  for (Vertex v0 : f.adjacentVertices()) {
    p[index] = geom.vertexPositions[v0];
    index++;
  }
  return findCircumcenter(p[0], p[1], p[2]);
}

// Returns the barycenter for faces incident on a nonflippable edge (e.g. a boundary edge), and the circumcenter for all
// other faces
Vector3 findODTCenter(VertexPositionGeometry& geom, Face f, MutationManager& mm) {
  Vector3 p0 = geom.vertexPositions[f.halfedge().tailVertex()];
  Vector3 p1 = geom.vertexPositions[f.halfedge().tipVertex()];
  Vector3 p2 = geom.vertexPositions[f.halfedge().next().tipVertex()];

  for (Edge e : f.adjacentEdges()) {
    if (e.isBoundary() || !mm.mayFlipEdge(e)) {
      // e is not flippable. return barycenter
      return (p0 + p1 + p2) / 3.;
    }
  }
  return findCircumcenter(p0, p1, p2);
}

bool isDelaunay(VertexPositionGeometry& geom, Edge e) {
  float angle1 = geom.cornerAngle(e.halfedge().next().next().corner());
  float angle2 = geom.cornerAngle(e.halfedge().twin().next().next().corner());
  return angle1 + angle2 <= PI;
}

bool isDelaunay_improv(VertexPositionGeometry& geom, Edge e) {
  Halfedge he = e.halfedge();
  Vector3 p0 = geom.vertexPositions[he.vertex()];
  Vector3 p1 = geom.vertexPositions[he.next().vertex()];
  Vector3 p2 = geom.vertexPositions[he.next().next().vertex()];
  Vector3 p3 = geom.vertexPositions[he.twin().next().next().vertex()];

  double la = (p0 - p1).norm();
  double lb = (p1 - p2).norm();
  double lc = (p2 - p0).norm();
  double ld = (p3 - p1).norm();
  double le = (p3 - p0).norm();

  return (lb * lb + lc * lc - la * la) / (lb * lc) + (ld * ld + le * le - la * la) / (ld * le) >= -0.1;
}


inline double diamondAngle(Vector3 a, Vector3 b, Vector3 c, Vector3 d) // dihedral angle at edge a-b
{
  Vector3 n1 = cross(b - a, c - a);
  Vector3 n2 = cross(b - d, a - d);
  return PI - angle(n1, n2);
}

inline bool checkFoldover(Vector3 a, Vector3 b, Vector3 c, Vector3 x, double angle) {
  return diamondAngle(a, b, c, x) < angle;
}

bool shouldCollapse(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, Edge e) {
  std::vector<Halfedge> edgesToCheck;
  Vertex v1 = e.halfedge().vertex();
  Vertex v2 = e.halfedge().twin().vertex();

  // find (halfedge) link around the edge, starting with those surrounding v1
  for (Halfedge he : v1.outgoingHalfedges()) {
    if (he.next().tailVertex() != v2 && he.next().tipVertex() != v2) {
      edgesToCheck.push_back(he.next());
    }
  }

  // link around v2
  for (Halfedge he : v2.outgoingHalfedges()) {
    if (he.next().tailVertex() != v1 && he.next().tipVertex() != v1) {
      edgesToCheck.push_back(he.next());
    }
  }

  // see if the point that would form after a collapse would cause a major foldover with surrounding edges
  Vector3 midpoint = edgeMidpoint(mesh, geom, e);
  // Vector3 butterfly = edgeButterfly(mesh,geom,e);
  for (Halfedge he0 : edgesToCheck) {
    Halfedge heT = he0.twin();
    Vertex v1 = heT.tailVertex();
    Vertex v2 = heT.tipVertex();
    Vertex v3 = heT.next().tipVertex();
    Vector3 a = geom.vertexPositions[v1];
    Vector3 b = geom.vertexPositions[v2];
    Vector3 c = geom.vertexPositions[v3];
    if (checkFoldover(a, b, c, midpoint, 2)) {
      return false;
    }
  }


  return true;
}

bool shouldCollapse(VertexData<double> SubsetSizingV, ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, Edge e,
                    RemeshOptions options) {
  std::vector<Halfedge> edgesToCheck;
  Vertex v1 = e.halfedge().vertex();
  Vertex v2 = e.halfedge().twin().vertex();

  // find (halfedge) link around the edge, starting with those surrounding v1
  for (Halfedge he : v1.outgoingHalfedges()) {
    if (he.next().tailVertex() != v2 && he.next().tipVertex() != v2) {
      edgesToCheck.push_back(he.next());
    }
  }

  // link around v2
  for (Halfedge he : v2.outgoingHalfedges()) {
    if (he.next().tailVertex() != v1 && he.next().tipVertex() != v1) {
      edgesToCheck.push_back(he.next());
    }
  }

  // see if the point that would form after a collapse would cause a major foldover with surrounding edges
  double area;
  double perimeter;
  double aspect;
  double a0;
  Vector3 midpoint = edgeMidpoint(mesh, geom, e);
  double E_metric;
  double NewSizing = std::max(SubsetSizingV[v1], SubsetSizingV[v2]);
  // Vector3 butterfly = edgeButterfly(mesh,geom,e);
  for (Halfedge he0 : edgesToCheck) {
    Halfedge heT = he0.twin();
    Vertex v1 = heT.tailVertex();
    Vertex v2 = heT.tipVertex();
    Vertex v3 = heT.next().tipVertex();
    Vector3 a = geom.vertexPositions[v1];
    Vector3 b = geom.vertexPositions[v2];
    Vector3 c = geom.vertexPositions[v3];
    if (checkFoldover(a, b, c, midpoint, 2)) {
      return false;
    }

    v1 = he0.next().tipVertex();
    v2 = he0.tailVertex();
    v3 = he0.tipVertex();

    a = midpoint;
    b = geom.vertexPositions[v2];
    c = geom.vertexPositions[v3];

    a0 = geom.faceArea(he0.face());
    area = 0.5 * norm(cross(b - a, c - a));
    perimeter = norm(b - a) + norm(c - a) + norm(c - b);
    aspect = 12 * sqrt(3) * area / (perimeter * perimeter);
    if ((area < a0 && area < 0.1 * options.min_absolute_length * options.min_absolute_length) ||
        aspect < options.aspect_min)
      return false;

    // We check the metric is not too big
    heT = he0;
    if (geom.edgeLength(heT.edge()) < 1e-10 || geom.edgeLength(heT.next().edge()) < 1e-10 ||
        geom.edgeLength(heT.next().next().edge()) < 1e-10)
      return false;

    // Ok we need to change this

    heT = heT.next();

    E_metric = norm(c - a) * sqrt((NewSizing + SubsetSizingV[heT.tailVertex()]) / 2.0);
    if (E_metric > 0.9) return false;

    heT = heT.next();
    E_metric = norm(b - a) * sqrt((NewSizing + SubsetSizingV[heT.tipVertex()]) / 2.0);
    if (E_metric > 0.9) return false;

    heT = heT.next();
    E_metric = norm(b - c) * sqrt((SubsetSizingV[heT.tailVertex()] + SubsetSizingV[heT.tipVertex()]) / 2.0);
    if (E_metric > 0.9) return false;
  }


  return true;
}

bool shouldCollapse(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, Edge e, RemeshOptions options) {
  std::vector<Halfedge> edgesToCheck;
  Vertex v1 = e.halfedge().vertex();
  Vertex v2 = e.halfedge().twin().vertex();

  // find (halfedge) link around the edge, starting with those surrounding v1
  for (Halfedge he : v1.outgoingHalfedges()) {
    if (he.next().tailVertex() != v2 && he.next().tipVertex() != v2) {
      edgesToCheck.push_back(he.next());
    }
  }

  // link around v2
  for (Halfedge he : v2.outgoingHalfedges()) {
    if (he.next().tailVertex() != v1 && he.next().tipVertex() != v1) {
      edgesToCheck.push_back(he.next());
    }
  }

  // see if the point that would form after a collapse would cause a major foldover with surrounding edges
  double area;
  double perimeter;
  double aspect;
  double a0;
  Vector3 midpoint = edgeMidpoint(mesh, geom, e);
  double E_metric;
  double NewSizing = std::max(geom.vertexSizing[v1], geom.vertexSizing[v2]);
  // Vector3 butterfly = edgeButterfly(mesh,geom,e);
  for (Halfedge he0 : edgesToCheck) {
    Halfedge heT = he0.twin();
    Vertex v1 = heT.tailVertex();
    Vertex v2 = heT.tipVertex();
    Vertex v3 = heT.next().tipVertex();
    Vector3 a = geom.vertexPositions[v1];
    Vector3 b = geom.vertexPositions[v2];
    Vector3 c = geom.vertexPositions[v3];
    if (checkFoldover(a, b, c, midpoint, 2)) {
      return false;
    }

    v1 = he0.next().tipVertex();
    v2 = he0.tailVertex();
    v3 = he0.tipVertex();

    a = midpoint;
    b = geom.vertexPositions[v2];
    c = geom.vertexPositions[v3];

    a0 = geom.faceAreas[he0.face()];
    area = 0.5 * norm(cross(b - a, c - a));
    perimeter = norm(b - a) + norm(c - a) + norm(c - b);
    aspect = 12 * sqrt(3) * area / (perimeter * perimeter);
    if ((area < a0 && area < 0.1 * options.min_absolute_length * options.min_absolute_length) ||
        aspect < options.aspect_min)
      return false;

    // We check the metric is not too big
    heT = he0;
    if (geom.edgeLengths[heT.edge()] < 1e-10 || geom.edgeLengths[heT.next().edge()] < 1e-10 ||
        geom.edgeLengths[heT.next().next().edge()] < 1e-10)
      return false;

    // Ok we need to change this

    heT = heT.next();

    E_metric = norm(c - a) * sqrt((NewSizing + geom.vertexSizing[heT.tailVertex()]) / 2.0);
    if (E_metric > 0.9) return false;

    heT = heT.next();
    E_metric = norm(b - a) * sqrt((NewSizing + geom.vertexSizing[heT.tipVertex()]) / 2.0);
    if (E_metric > 0.9) return false;

    heT = heT.next();
    E_metric = norm(b - c) * sqrt((geom.vertexSizing[heT.tailVertex()] + geom.vertexSizing[heT.tipVertex()]) / 2.0);
    if (E_metric > 0.9) return false;
  }


  return true;
}

bool shouldCollapseSimple(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, Edge e, RemeshOptions options) {
  std::vector<Halfedge> edgesToCheck;
  Vertex v1 = e.halfedge().vertex();
  Vertex v2 = e.halfedge().twin().vertex();

  // find (halfedge) link around the edge, starting with those surrounding v1
  for (Halfedge he : v1.outgoingHalfedges()) {
    if (he.next().tailVertex() != v2 && he.next().tipVertex() != v2) {
      edgesToCheck.push_back(he.next());
    }
  }
  // link around v2
  for (Halfedge he : v2.outgoingHalfedges()) {
    if (he.next().tailVertex() != v1 && he.next().tipVertex() != v1) {
      edgesToCheck.push_back(he.next());
    }
  }

  // see if the point that would form after a collapse would cause a major foldover with surrounding edges
  Vector3 midpoint = edgeMidpoint(mesh, geom, e);

  // Vector3 butterfly = edgeButterfly(mesh,geom,e);
  for (Halfedge he0 : edgesToCheck) {
    Halfedge heT = he0.twin();
    Vertex v1 = heT.tailVertex();
    Vertex v2 = heT.tipVertex();
    Vertex v3 = heT.next().tipVertex();
    Vector3 a = geom.vertexPositions[v1];
    Vector3 b = geom.vertexPositions[v2];
    Vector3 c = geom.vertexPositions[v3];
    if (checkFoldover(a, b, c, midpoint, 2)) {
      return false;
    }
  }

  return true;
}

// Warning: requires that geom.vertexDualAreas and geom.vertexMeanCurvatures are filled in with accurate data
double getSmoothMeanCurvature(VertexPositionGeometry& geom, Vertex v) {
  double A = geom.vertexDualAreas[v];
  double S = geom.vertexMeanCurvatures[v];
  double K = S / A;
  return K;
}

// flatLength: specifies how long the target edge length should be in flat regions
// curvatureAdaptation: controls how much variation in target length occurs due to curvature
double findMeanTargetL(SurfaceMesh& mesh, VertexPositionGeometry& geom, Edge e, double flatLength,
                       double curvatureAdaptation) {
  double averageH = 0;
  for (Vertex v : e.adjacentVertices()) {
    averageH += getSmoothMeanCurvature(geom, v);
  }
  averageH /= 2;
  double L = flatLength / (fabs(averageH) * curvatureAdaptation + 1);
  return L;
}


size_t fixDelaunay(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom) {
  MutationManager mm(mesh, geom);
  return fixDelaunay(mesh, geom, mm);
}

size_t fixDelaunay(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm) {
  // Logic duplicated from surface/intrinsic_triangulation.cpp

  std::deque<Edge> edgesToCheck;      // queue of edges to check if Delaunay
  EdgeData<bool> inQueue(mesh, true); // true if edge is currently in edgesToCheck

  // start with all edges
  for (Edge e : mesh.edges()) {
    edgesToCheck.push_back(e);
  }

  // counter and limit for number of flips
  size_t flipMax = 100 * mesh.nVertices();
  size_t nFlips = 0;
  while (!edgesToCheck.empty() && nFlips < flipMax) {
    Edge e = edgesToCheck.front();
    edgesToCheck.pop_front();
    inQueue[e] = false;

    if (e.isBoundary() || isDelaunay_improv(geom, e)) continue;

    // if not Delaunay, try to flip edge
    bool wasFlipped = mm.flipEdge(e);

    if (!wasFlipped) continue;

    nFlips++;

    // Add neighbors to queue, as they may need flipping now
    Halfedge he = e.halfedge();
    std::array<Edge, 4> neighboringEdges{he.next().edge(), he.next().next().edge(), he.twin().next().edge(),
                                         he.twin().next().next().edge()};
    for (Edge nE : neighboringEdges) {
      if (!inQueue[nE]) {
        edgesToCheck.push_back(nE);
        inQueue[nE] = true;
      }
    }
  }
  return nFlips;
}

size_t fixDelaunay(std::vector<Edge> activeEdges, ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
                   MutationManager& mm) {
  // Logic duplicated from surface/intrinsic_triangulation.cpp

  std::deque<Edge> edgesToCheck;      // queue of edges to check if Delaunay
  EdgeData<bool> inQueue(mesh, true); // true if edge is currently in edgesToCheck

  // start with all edges
  for (Edge e : activeEdges) {
    edgesToCheck.push_back(e);
  }

  // counter and limit for number of flips
  size_t flipMax = 100 * mesh.nVertices();
  size_t nFlips = 0;
  while (!edgesToCheck.empty() && nFlips < flipMax) {
    Edge e = edgesToCheck.front();
    edgesToCheck.pop_front();
    inQueue[e] = false;

    if (e.isBoundary() || isDelaunay_improv(geom, e)) continue;

    // if not Delaunay, try to flip edge
    bool wasFlipped = mm.flipEdge(e);

    if (!wasFlipped) continue;

    nFlips++;

    // Add neighbors to queue, as they may need flipping now
    Halfedge he = e.halfedge();
    std::array<Edge, 4> neighboringEdges{he.next().edge(), he.next().next().edge(), he.twin().next().edge(),
                                         he.twin().next().next().edge()};
    for (Edge nE : neighboringEdges) {
      if (!inQueue[nE]) {
        edgesToCheck.push_back(nE);
        inQueue[nE] = true;
      }
    }
  }
  return nFlips;
}

double smoothByLaplacian(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, double stepSize,
                         RemeshBoundaryCondition bc) {
  MutationManager mm(mesh, geom);
  return smoothByLaplacian(mesh, geom, mm, stepSize, bc);
}

// doub;e smoothByLaplacian(ManifoldSurfaceMesh& mesh, VertexPosi)

double smoothByLaplacian(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm, double stepSize,
                         RemeshBoundaryCondition bc) {
  VertexData<Vector3> vertexOffsets(mesh);
  geom.requireVertexNormals();
  geom.requireVertexPositions();

  for (Vertex v : mesh.vertices()) {
    // calculate average of surrounding vertices
    Vector3 avgNeighbor = Vector3::zero();
    for (Vertex j : v.adjacentVertices()) {
      avgNeighbor += geom.vertexPositions[j];
    }
    avgNeighbor /= v.degree();

    Vector3 updateDirection = avgNeighbor - geom.vertexPositions[v];

    // project updateDirection onto space of allowed movements
    Vector3 stepDir;
    if (v.isBoundary()) {
      switch (bc) {
      case RemeshBoundaryCondition::Fixed:
        stepDir = Vector3::zero();
        break;
      case RemeshBoundaryCondition::Tangential:
        // for free boundary vertices, project the average to the boundary tangent line
        stepDir = projectToLine(updateDirection, boundaryVertexTangent(geom, v, mm));
        break;
      case RemeshBoundaryCondition::Free:
        // for free boundary vertices, project the average to the surface tangent plane
        stepDir = projectToPlane(updateDirection, vertexNormal(geom, v, mm));
        break;
      }
    } else {
      // for interior vertices, project the average to the tangent plane
      stepDir = projectToPlane(updateDirection, vertexNormal(geom, v, mm));
    }
    vertexOffsets[v] = stepSize * stepDir;
  }

  // update final vertices
  double totalMovement = 0;
  for (Vertex v : mesh.vertices()) {
    bool didMove = mm.repositionVertex(v, vertexOffsets[v]);
    if (didMove) {
      totalMovement += vertexOffsets[v].norm();
    }
  }
  // geom.inputVertexPositions = geom.vertexPositions;
  geom.unrequireVertexNormals();
  geom.unrequireVertexPositions();
  // What i would like to check is if there is any nAN VERTEX AFTERTHIS
  return totalMovement / mesh.nVertices();
}

double smoothByCircumcenter(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, double stepSize,
                            RemeshBoundaryCondition bc) {
  MutationManager mm(mesh, geom);
  return smoothByCircumcenter(mesh, geom, mm, stepSize, bc);
}

double smoothByCircumcenter(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm,
                            double stepSize, RemeshBoundaryCondition bc) {
  geom.requireFaceAreas();
  VertexData<Vector3> vertexOffsets(mesh);
  for (Vertex v : mesh.vertices()) {
    Vector3 updateDirection = Vector3::zero();
    for (Face f : v.adjacentFaces()) {
      // add the circumcenter weighted by face area to the update direction
      Vector3 circum = findODTCenter(geom, f, mm);
      updateDirection += geom.faceArea(f) * (circum - geom.vertexPositions[v]);
    }
    updateDirection /= (6 * geom.vertexDualArea(v));

    // project updateDirection onto space of allowed movements
    Vector3 stepDir;
    if (v.isBoundary()) {
      switch (bc) {
      case RemeshBoundaryCondition::Fixed:
        stepDir = Vector3::zero();
        break;
      case RemeshBoundaryCondition::Tangential:
        // for free boundary vertices, project the average to the boundary tangent line
        stepDir = projectToLine(updateDirection, boundaryVertexTangent(geom, v, mm));
        break;
      case RemeshBoundaryCondition::Free:
        // for free boundary vertices, project the average to the surface tangent plane
        stepDir = projectToPlane(updateDirection, vertexNormal(geom, v, mm));
        break;
      }
    } else {
      // for interior vertices, project the average to the tangent plane
      stepDir = projectToPlane(updateDirection, vertexNormal(geom, v, mm));
    }
    vertexOffsets[v] = stepSize * stepDir;
  }

  // update final vertices
  double totalMovement = 0;
  for (Vertex v : mesh.vertices()) {
    bool didMove = mm.repositionVertex(v, vertexOffsets[v]);
    if (didMove) {
      totalMovement += vertexOffsets[v].norm();
    }
  }
  return totalMovement / mesh.nVertices();
} // namespace surface


bool adjustEdgeLengths(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, RemeshOptions options) {
  MutationManager mm(mesh, geom);
  return adjustEdgeLengths(mesh, geom, mm, options);
}

bool adjustEdgeLengths(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm,
                       RemeshOptions options) {
  geom.requireVertexDualAreas();
  geom.requireVertexSizing();
  geom.requireEdgeLengths();
  geom.requireVertexPositions();

  bool didSplitOrCollapse = false;
  // queues of edges to CHECK to change
  std::vector<Edge> toSplit;
  std::vector<Edge> toCollapse;

  if (options.no_remesh_list) {
    for (Edge e : mesh.edges()) {
      if (options.No_remesh_list[e] == 0) {
        toSplit.push_back(e);
      }
    }
  } else {
    for (Edge e : mesh.edges()) {
      toSplit.push_back(e);
    }
  }

  size_t counter_vertex = 0;
  // actually splitting
  double E_sizing;
  double newSizing;
  while (!toSplit.empty()) {
    Edge e = toSplit.back();
    toSplit.pop_back();
    double length_e = geom.edgeLength(e);

    newSizing = (geom.vertexSizing[e.halfedge().vertex()] + geom.vertexSizing[e.halfedge().twin().vertex()]) / (2.0);
    E_sizing = length_e * sqrt(newSizing);

    if (E_sizing > 1) {

      Vector3 newPos = edgeMidpoint(mesh, geom, e);
      Halfedge he = mm.splitEdge(e, newPos);
      if (he != Halfedge()) {

        geom.vertexSizing[he.vertex()] = newSizing;

        Halfedge heround = he;
        int counter = 0;
        do {
          counter += 1;
          if (geom.edgeLengths[heround.edge()] < 1e-10) {
            geom.edgeLengths[heround.edge()] =
                norm(geom.vertexPositions[heround.vertex()] - geom.vertexPositions[heround.tipVertex()]);
          }
          heround = heround.twin().next();
        } while (heround != he);

        geom.vertexSizing[he.vertex()] = newSizing;
        geom.inputVertexPositions[he.vertex()] = newPos;
      }

      didSplitOrCollapse = true;
    }

    else {
      if (options.no_remesh_list) {

        if (options.No_remesh_list_v[e.firstVertex()] == 0 && options.No_remesh_list_v[e.secondVertex()] == 0) {
          toCollapse.push_back(e);
        } else {
          counter_vertex += 1;
        }

      } else {
        toCollapse.push_back(e);
      }
    }
  }

  while (!toCollapse.empty()) {

    Edge e = toCollapse.back();
    toCollapse.pop_back();
    if (e == Edge() || e.isDead()) continue; // make sure it exists

    // Now we do
    if (geom.edgeLengths[e] < options.max_absolute_length) {
      if (geom.edgeLengths[e] < 1e-4) std::cout << "THe edgelength check says " << geom.edgeLengths[e] << " \n";
      Vector3 newPos = edgeMidpoint(mesh, geom, e);
      newSizing = std::max(geom.vertexSizing[e.halfedge().tipVertex()], geom.vertexSizing[e.halfedge().tailVertex()]);
      if (shouldCollapse(mesh, geom, e, options)) {
        Vertex v = mm.collapseEdge(e, newPos);
        if (v != Vertex()) {
          geom.vertexSizing[v] = newSizing;
          didSplitOrCollapse = true;
        }
      }
    }
  }
  geom.unrequireEdgeLengths();
  geom.unrequireVertexDualAreas();
  geom.unrequireVertexSizing();
  geom.unrequireVertexPositions();
  mesh.compress();
  return didSplitOrCollapse;
}

double Face_sizing(Face f, VertexPositionGeometry& geom) {
  Eigen::Matrix2d Sizing;
  Sizing << 0.0, 0.0, 0.0, 0.0;
  geom.requireVertexPositions();
  Eigen::Matrix3d Rotation;
  Vector3 Normal1 = geom.faceNormal(f);
  Vector3 Axis1({0.0, 1.0, 0.0});
  Eigen::Vector3d Normal{Normal1.x, Normal1.y, Normal1.z};
  Eigen::Vector3d Axis{Axis1.x, Axis1.y, Axis1.z};

  Eigen::Vector3d uvw = Normal.cross(Axis);

  double rcos = Normal.dot(Axis);
  double rsin = uvw.norm();

  if (rsin > 1e-7) {
    uvw /= rsin;
  }

  Eigen::Matrix3d V_x;

  V_x << 0, -1 * uvw(2), uvw(1), uvw(2), 0, -1 * uvw(0), -1 * uvw(1), uvw(0), 0;

  Rotation = rcos * Eigen::Matrix3d::Identity() + rsin * V_x + uvw * uvw.transpose() * (1 - rcos);
  Vector3 Pos;
  Eigen::Vector3d spare_vec;
  std::vector<Eigen::Vector2d> Positions(0);
  std::vector<int> index(0);

  for (Vertex v : f.adjacentVertices()) {

    Pos = geom.vertexPositions[v];
    spare_vec = Eigen::Vector3d{Pos.x, Pos.y, Pos.z};
    spare_vec = Rotation * spare_vec;
    index.push_back(v.getIndex());
    Positions.push_back({spare_vec(0), spare_vec(2)});
  }
  int counter = 0;
  for (Edge e : f.adjacentEdges()) {
    Eigen::Vector2d e_mat = Positions[(counter + 1) % 3] - Positions[counter % 3];
    Eigen::Vector2d t_mat{-1 * e_mat(1), e_mat(0)};
    t_mat.normalize();
    double theta = geom.dihedralAngle(e.halfedge()); // TODO adapt
    if (e.halfedge().face() != f) {
      theta = geom.dihedralAngle(e.halfedge().twin()); // TOOD adapt
    }
    Sizing -= 1 / 2. * theta * e_mat.norm() * t_mat * t_mat.transpose();
    counter += 1;
  }

  Sizing /= geom.faceArea(f);
  Sizing = Sizing.transpose() * Sizing; // normalize by refine angle, so that the resulting sizing is in
  Eigen::EigenSolver<Eigen::Matrix2d> Solve;
  Eigen::Vector2d eigenvalues;
  Solve.compute(Sizing);
  eigenvalues = Solve.eigenvalues().real();

  return std::max(fabs(eigenvalues[0]), fabs(eigenvalues[1]));
}

bool independent(const Edge edge, const std::vector<Edge> edges) {
  for (size_t i = 0; i < edges.size(); i++) {
    Edge edge1 = edges[i];
    if (edge.firstVertex() == edge1.firstVertex() || edge.firstVertex() == edge1.secondVertex() ||
        edge.secondVertex() == edge1.firstVertex() || edge.secondVertex() == edge1.secondVertex())
      return false;
  }
  return true;
}
std::vector<Edge> independent_edges(const std::vector<Edge> edges) {
  std::vector<Edge> iedges;
  for (size_t e = 0; e < edges.size(); e++)
    if (independent(edges[e], iedges)) iedges.push_back(edges[e]);
  return iedges;
}

bool collapseSubset(std::vector<Face> activeFaces, ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
                    MutationManager& mm, RemeshOptions options) {


  FaceData<double> SubsetSizingF(mesh, 0.0);
  VertexData<double> SubsetSizingV(mesh, 0.0);
  std::vector<Vertex> activeVertices;
  std::vector<Face> neededFaces;
  std::vector<Edge> toCollapse;
  std::map<Face, int> fMap;
  std::map<Vertex, int> vMap;
  std::map<Edge, int> eMap;

  int idx = 0;
  int idx2 = 0;
  int idx3 = 0;

  bool didCollapse = false;
  // queues of edges to CHECK to change
  for (Face f : activeFaces) {
    for (Vertex v : f.adjacentVertices()) {
      vMap[v] = idx++;
      // I NEED TO DO AN EMAP SO WE DONT OVERRUN
      for (Face f1 : v.adjacentFaces()) {
        fMap[f1] = idx++;
      }
    }
  }
  for (Face f : activeFaces) {
    for (Edge e : f.adjacentEdges()) {
      eMap[e] = idx3++;
    }
  }
  for (auto& search : fMap) {
    neededFaces.push_back(search.first);
  }
  for (auto& search : vMap) {
    activeVertices.push_back(search.first);
  }
  for (auto& search : eMap) {
    toCollapse.push_back(search.first);
  }

  // std::cout << "Passed all the maps\n";

  for (Face f : neededFaces)
    SubsetSizingF[f] = clamp(Face_sizing(f, geom) / (options.refine_angle * options.refine_angle),
                             1.0 / (options.max_absolute_length * options.max_absolute_length),
                             1.0 / (options.min_absolute_length * options.min_absolute_length));


  double sizing = 0;
  for (Vertex v : activeVertices) {
    sizing = 0;
    for (Face f : v.adjacentFaces()) {
      sizing += geom.faceArea(f) * SubsetSizingF[f] / 3.0;
      if (SubsetSizingF[f] < 1e-1) std::cout << "THis sizing is " << SubsetSizingF[f] << " \n";
    }
    SubsetSizingV[v] = sizing / geom.vertexDualArea(v);
  }

  std::cout << "Passed the sizings\n";

  size_t counter_vertex = 0;

  double E_sizing;
  double newSizing;


  while (!toCollapse.empty()) {

    Edge e = toCollapse.back();
    toCollapse.pop_back();
    if (e == Edge() || e.isDead()) continue; // make sure it exists

    // Now we do
    if (geom.edgeLength(e) < options.max_absolute_length) {
      // std::cout << "Lengthcheck\n";
      Vector3 newPos = edgeMidpoint(mesh, geom, e);
      newSizing = (SubsetSizingV[e.halfedge().vertex()] + SubsetSizingV[e.halfedge().twin().vertex()]) / (2.0);
      // std::cout << "Calling should collapse\n";
      if (shouldCollapse(SubsetSizingV, mesh, geom, e, options)) {
        Vertex v = mm.collapseEdge(e, newPos);
        if (v != Vertex()) {
          options.numberOp += 1;
          // geom.vertexSizing[v] = newSizing;
          std::vector<Face> active_faces;
          // std::cout << "Flipping subset\n";
          for (Face f : v.adjacentFaces()) active_faces.push_back(f);
          flipSubset(active_faces, mesh, geom, mm, options);
          didCollapse = true;
        }
      }
    }
  }

  return didCollapse;
}

std::vector<Face> splitSubset(std::vector<Face> activeFaces, ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
                              MutationManager& mm, RemeshOptions options) {

  FaceData<double> SubsetSizingF(mesh, 0.0);
  VertexData<double> SubsetSizingV(mesh, 0.0);
  std::vector<Vertex> activeVertices;
  std::vector<Face> neededFaces;
  std::vector<Edge> toSplit;
  std::map<Face, int> fMap;
  std::map<Vertex, int> vMap;
  std::map<Edge, int> eMap;
  std::vector<Face> facesForward;
  int idx = 0;
  int idx2 = 0;
  int idx3 = 0;
  geom.requireVertexPositions();
  for (Face f : activeFaces) {
    for (Vertex v : f.adjacentVertices()) {
      vMap[v] = idx++;
      // I NEED TO DO AN EMAP SO WE DONT OVERRUN
      for (Face f1 : v.adjacentFaces()) {
        fMap[f1] = idx++;
      }
    }
  }
  for (Face f : activeFaces) {
    for (Edge e : f.adjacentEdges()) {
      eMap[e] = idx3++;
    }
  }

  for (auto& search : fMap) {
    neededFaces.push_back(search.first);
  }
  for (auto& search : vMap) {
    activeVertices.push_back(search.first);
  }
  for (auto& search : eMap) {
    toSplit.push_back(search.first);
  }
  // Ok so now i need to calculate the sizings

  for (Face f : neededFaces) {
    SubsetSizingF[f] = clamp(Face_sizing(f, geom) / (options.refine_angle * options.refine_angle),
                             1.0 / (options.max_absolute_length * options.max_absolute_length),
                             1.0 / (options.min_absolute_length * options.min_absolute_length));
    // std::cout << "THe sizing is" << SubsetSizingF[f] << " \n";
  }
  // NOw that we have all the needed faces, we need to do the vertexsizing
  double sizing = 0;
  for (Vertex v : activeVertices) {
    sizing = 0;
    for (Face f : v.adjacentFaces()) {
      sizing += geom.faceArea(f) * SubsetSizingF[f] / 3.0;
      if (SubsetSizingF[f] < 1e-1) std::cout << "THis sizing is " << SubsetSizingF[f] << " \n";
    }
    SubsetSizingV[v] = sizing / geom.vertexDualArea(v);
    // std::cout << "THe vertex sizing is" << SubsetSizingV[v] << " \n";
  }
  // We have all the things we need

  // OK so we should have unique faces now

  bool didSplit = false;

  double newSizing;
  while (!toSplit.empty()) {
    Edge e = toSplit.back();
    toSplit.pop_back();

    double length_e = geom.edgeLength(e);
    // Ok i have the edgelength
    // std::cout << "got edge length\n";

    newSizing = (SubsetSizingV[e.halfedge().vertex()] + SubsetSizingV[e.halfedge().twin().vertex()]) / (2.0);

    double E_sizing = length_e * sqrt(newSizing);
    if (E_sizing > 1.0) {
      Vector3 newPos = edgeMidpoint(mesh, geom, e);
      Halfedge he = mm.splitEdge(e, newPos);
      if (he != Halfedge()) {
        // SubsetertexSizing[he.vertex()] = newSizing;
        SubsetSizingV[he.vertex()] = newSizing;
        Halfedge heround = he;
        int counter = 0;
        do {
          facesForward.push_back(heround.face());
          counter += 1;
          // if (geom.edgeLength(heround.edge()) < 1e-10) {
          //   geom.edgeLength(heround.edge()] =
          //       norm(geom.vertexPositions[heround.vertex()] -
          //            geom.vertexPositions[heround.tipVertex()]); // do we really need this?
          // }
          heround = heround.twin().next();
        } while (heround != he);

        didSplit = true;
        flipSubset(facesForward, mesh, geom, mm, options);
      }
    }
    {
      // If the sizing is smaller than 1, we can still collapse
      facesForward.push_back(e.halfedge().face());
      facesForward.push_back(e.halfedge().twin().face());
    }
  }
  // geom.unrequireEdgeLengths();
  geom.unrequireVertexPositions();

  return facesForward;
}
int flipSubset(std::vector<Face> active, ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm,
               RemeshOptions options) {
  std::vector<Edge> edges;
  std::map<Edge, int> eMap;
  int idx = 0;
  Halfedge he;
  for (size_t f = 0; f < active.size(); f++) {
    he = active[f].halfedge();
    for (int e = 0; e < 3; e++) {
      eMap[he.edge()] = idx++;
      he = he.next();
    }
  }
  for (auto& search : eMap) {
    edges.push_back(search.first);
  }
  std::vector<Edge> indepEdges = independent_edges(edges);

  int flips = fixDelaunay(indepEdges, mesh, geom, mm);
  options.numberOp += flips;
  return flips;
}


struct Deterministic_sort {
  inline bool operator()(const std::pair<double, Edge>& left, const std::pair<double, Edge>& right) {
    return left.first < right.first;
  }
} deterministic_sort;


std::vector<Edge> findBadEdges(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm,
                               RemeshOptions options) {

  std::vector<std::pair<double, Edge>> edgems;
  for (Edge e : mesh.edges()) {
    // TODO add the no remesh option here (not necessary yet) This is the spirit
    // if (options.no_remesh_list) {
    //   for (Edge e : mesh.edges()) {

    //     if (options.No_remesh_list[e] == 0) {

    //       toSplit.push_back(e);
    //     }
    //   }
    // } else {
    //   for (Edge e : mesh.edges()) {
    //     toSplit.push_back(e);
    //   }
    // }
    double E_sizing = geom.edgeLengths[e] *
                      sqrt(geom.vertexSizing[e.halfedge().vertex()] + geom.vertexSizing[e.halfedge().twin().vertex()]) /
                      (2.0);
    if (E_sizing > 1) edgems.push_back(std::make_pair(E_sizing, e));
  }

  std::sort(edgems.begin(), edgems.end(), deterministic_sort);
  std::vector<Edge> edges(edgems.size());
  for (size_t e = 0; e < edgems.size(); e++) {
    edges[e] = edgems[edgems.size() - e - 1].second;
  }
  return edges;
}

bool splitWorstEdges(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm,
                     RemeshOptions options) {
  geom.requireVertexDualAreas();
  geom.requireVertexSizing();
  geom.requireEdgeLengths();

  bool didSplit = false;
  std::vector<Edge> toSplit = findBadEdges(mesh, geom, mm, options);
  std::vector<Face> activeFaces;

  double newSizing;
  while (!toSplit.empty()) {
    Edge e = toSplit.back();
    toSplit.pop_back();
    double length_e = geom.edgeLength(e);

    newSizing = (geom.vertexSizing[e.halfedge().vertex()] + geom.vertexSizing[e.halfedge().twin().vertex()]) / (2.0);

    Vector3 newPos = edgeMidpoint(mesh, geom, e);
    Halfedge he = mm.splitEdge(e, newPos);
    if (he != Halfedge()) {
      geom.vertexSizing[he.vertex()] = newSizing;
      Halfedge heround = he;
      int counter = 0;
      do {
        activeFaces.push_back(heround.face());
        counter += 1;
        if (geom.edgeLengths[heround.edge()] < 1e-10) {
          geom.edgeLengths[heround.edge()] =
              norm(geom.vertexPositions[heround.vertex()] - geom.vertexPositions[heround.tipVertex()]);
        }
        heround = heround.twin().next();
      } while (heround != he);

      didSplit = true;
      options.numberOp += 1;
      // flipSubset(activeFaces, mesh, geom, mm, options);
    }
  }

  // We splitted everything
  // return;
  // std::cout<<"THe number of edges that i shall not collapse are"<<counter_vertex<<" this number shouuld be
  // constant?\n"; actually collapsing
  geom.unrequireEdgeLengths();
  geom.unrequireVertexDualAreas();
  geom.unrequireVertexSizing();
  mesh.compress();
  return didSplit;
}

bool improveFaces(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm, RemeshOptions options) {
  geom.requireVertexDualAreas();
  geom.requireVertexSizing();
  geom.requireEdgeLengths();

  bool didCollapse = false;
  // queues of edges to CHECK to change
  std::vector<Edge> toCollapse;

  if (options.no_remesh_list) {
    for (Edge e : mesh.edges()) {

      if (options.No_remesh_list[e] == 0) {
        // std::cout<<"Not remeshing this edge "<<e.getIndex() <<"\n";
        toCollapse.push_back(e);
      }
    }
  } else {
    for (Edge e : mesh.edges()) {
      toCollapse.push_back(e);
    }
  }

  size_t counter_vertex = 0;
  // actually splitting
  double E_sizing;
  double newSizing;


  while (!toCollapse.empty()) {

    Edge e = toCollapse.back();
    toCollapse.pop_back();
    if (e == Edge() || e.isDead()) continue; // make sure it exists

    // Now we do
    if (geom.edgeLengths[e] < options.max_absolute_length) {
      if (geom.edgeLengths[e] < 1e-4) std::cout << "THe edgelength check says " << geom.edgeLengths[e] << " \n";
      Vector3 newPos = edgeMidpoint(mesh, geom, e);
      newSizing = std::max(geom.vertexSizing[e.halfedge().tipVertex()], geom.vertexSizing[e.halfedge().tailVertex()]);
      if (shouldCollapse(mesh, geom, e, options)) {
        Vertex v = mm.collapseEdge(e, newPos);
        if (v != Vertex()) {
          options.numberOp += 1;
          geom.vertexSizing[v] = newSizing;
          // std::vector<Face> active_faces;
          // for (Face f : v.adjacentFaces()) active_faces.push_back(f);
          // flipSubset(active_faces, mesh, geom, mm, options);
          didCollapse = true;
        }
      }
    }
  }
  geom.unrequireEdgeLengths();
  geom.unrequireVertexDualAreas();
  geom.unrequireVertexSizing();
  mesh.compress();
  return didCollapse;
}
void remeshSmallAngles(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, RemeshOptions options) {

  MutationManager mm = MutationManager(mesh, geom);
  bool didCollapse = false;
  geom.requireCornerAngles();
  geom.requireEdgeLengths();
  std::vector<Edge> toCollapse;
  Halfedge he;
  Edge e;
  for (Corner c : mesh.corners()) {

    if (geom.cornerAngles[c] < options.angleThresh) {
      he = c.halfedge();
      e = he.edge();
      he = he.next();
      if (geom.edgeLengths[e] > geom.edgeLengths[he.edge()]) e = he.edge();
      he = he.next();
      if (geom.edgeLengths[e] > geom.edgeLengths[he.edge()]) e = he.edge();
      // Now lets do the following
      toCollapse.push_back(e);
    }

    // I iterate over the corners
  }
  std::vector<Edge> toCollapseIndep = independent_edges(toCollapse);

  while (!toCollapseIndep.empty()) {
    e = toCollapseIndep.back();
    toCollapseIndep.pop_back();
    Vector3 newPos = edgeMidpoint(mesh, geom, e);
    if (e == Edge() || e.isDead()) continue;
    if (shouldCollapseSimple(mesh, geom, e, options)) {
      Vertex v = mm.collapseEdge(e, newPos);
      if (v != Vertex()) {
        didCollapse = true;
      }
    }
  }
  geom.unrequireEdgeLengths();
  geom.unrequireCornerAngles();
  if (didCollapse) mesh.compress();
}

} // namespace surface
} // namespace geometrycentral