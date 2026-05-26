#include "geometrycentral/surface/extrinsic_geometry_interface.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <limits>

namespace geometrycentral {
namespace surface {

// clang-format off
ExtrinsicGeometryInterface::ExtrinsicGeometryInterface(SurfaceMesh& mesh_) : 
  IntrinsicGeometryInterface(mesh_),

  edgeDihedralAnglesQ                  (&edgeDihedralAngles,                   std::bind(&ExtrinsicGeometryInterface::computeEdgeDihedralAngles, this),                 quantities),
  vertexMeanCurvaturesQ                (&vertexMeanCurvatures,                 std::bind(&ExtrinsicGeometryInterface::computeVertexMeanCurvatures, this),               quantities),
  vertexMinPrincipalCurvaturesQ        (&vertexMinPrincipalCurvatures,         std::bind(&ExtrinsicGeometryInterface::computeVertexMinPrincipalCurvatures, this),       quantities),
  vertexMaxPrincipalCurvaturesQ        (&vertexMaxPrincipalCurvatures,         std::bind(&ExtrinsicGeometryInterface::computeVertexMaxPrincipalCurvatures, this),       quantities),
  vertexPrincipalCurvatureDirectionsQ  (&vertexPrincipalCurvatureDirections,   std::bind(&ExtrinsicGeometryInterface::computeVertexPrincipalCurvatureDirections, this), quantities),
  facePrincipalCurvatureDirectionsQ    (&facePrincipalCurvatureDirections,     std::bind(&ExtrinsicGeometryInterface::computeFacePrincipalCurvatureDirections, this),   quantities),
  faceSizingQ                          (&faceSizing,                           std::bind(&ExtrinsicGeometryInterface::computeFaceSizing, this),                         quantities),
  vertexSizingQ                        (&vertexSizing,                         std::bind(&ExtrinsicGeometryInterface::computeVertexSizing, this),                     quantities)
  
  {
  }
// clang-format on

// Edge dihedral angle
void ExtrinsicGeometryInterface::requireEdgeDihedralAngles() { edgeDihedralAnglesQ.require(); }
void ExtrinsicGeometryInterface::unrequireEdgeDihedralAngles() { edgeDihedralAnglesQ.unrequire(); }


void ExtrinsicGeometryInterface::computeVertexMeanCurvatures() {
  edgeLengthsQ.ensureHave();
  edgeDihedralAnglesQ.ensureHave();

  vertexMeanCurvatures = VertexData<double>(mesh);

  // Here we use the Steiner approximation of mean curvature: if we
  // imagine each edge is an arc of a small cylinder connecting the
  // two adjacent faces, and let the radius of this cylinder go to
  // zero, then the mean curvature is one half the dihedral angle
  // times the edge length.


  for (Vertex v : mesh.vertices()) {
    double meanCurvature = 0.;
    for (Halfedge he : v.outgoingHalfedges()) {
      double len = edgeLengths[he.edge()];
      double alpha = edgeDihedralAngles[he.edge()];
      meanCurvature += alpha * len / 2.;
    }

    // Since each edge is shared by two vertices, we give half of
    // the curvature to each vertex (equivalently, we integrate
    // curvature over a dual cell associated with the vertex).
    // The resulting vertex curvatures are equal to 1 for a unit sphere.
    vertexMeanCurvatures[v] = meanCurvature / 2.;
  }
}
void ExtrinsicGeometryInterface::requireVertexMeanCurvatures() { vertexMeanCurvaturesQ.require(); }
void ExtrinsicGeometryInterface::unrequireVertexMeanCurvatures() { vertexMeanCurvaturesQ.unrequire(); }


void ExtrinsicGeometryInterface::computeVertexMinPrincipalCurvatures() {
  computePrincipalCurvatures(1, vertexMinPrincipalCurvatures);
}
void ExtrinsicGeometryInterface::requireVertexMinPrincipalCurvatures() { vertexMinPrincipalCurvaturesQ.require(); }
void ExtrinsicGeometryInterface::unrequireVertexMinPrincipalCurvatures() { vertexMinPrincipalCurvaturesQ.unrequire(); }

void ExtrinsicGeometryInterface::computeVertexMaxPrincipalCurvatures() {
  computePrincipalCurvatures(2, vertexMaxPrincipalCurvatures);
}
void ExtrinsicGeometryInterface::requireVertexMaxPrincipalCurvatures() { vertexMaxPrincipalCurvaturesQ.require(); }
void ExtrinsicGeometryInterface::unrequireVertexMaxPrincipalCurvatures() { vertexMaxPrincipalCurvaturesQ.unrequire(); }

void ExtrinsicGeometryInterface::computePrincipalCurvatures(int whichCurvature, VertexData<double>& kappa) {
  vertexGaussianCurvaturesQ.ensureHave();
  vertexMeanCurvaturesQ.ensureHave();
  vertexDualAreasQ.ensureHave();

  kappa = VertexData<double>(mesh);

  for (Vertex v : mesh.vertices()) {
    // Vertex mean and Gaussian curvatures are integrated
    // values; need to divide them by dual areas to get
    // pointwise quantities.
    double A = vertexDualAreas[v];
    double H = vertexMeanCurvatures[v] / A;
    double K = vertexGaussianCurvatures[v] / A;

    // The two principal curvatures are given by
    //    H +/- sqrt( H^2 - K )
    double c = std::sqrt(std::max(0., H * H - K));
    double k1 = H - c;
    double k2 = H + c;

    if (whichCurvature == 1)
      kappa[v] = std::min(k1, k2);
    else
      kappa[v] = std::max(k1, k2);
  }
}


void ExtrinsicGeometryInterface::computeVertexPrincipalCurvatureDirections() {
  edgeLengthsQ.ensureHave();
  halfedgeVectorsInVertexQ.ensureHave();
  edgeDihedralAnglesQ.ensureHave();

  vertexPrincipalCurvatureDirections = VertexData<Vector2>(mesh);

  for (Vertex v : mesh.vertices()) {
    Vector2 principalDir{0.0, 0.0};
    for (Halfedge he : v.outgoingHalfedges()) {
      double len = edgeLengths[he.edge()];
      double alpha = edgeDihedralAngles[he.edge()];
      Vector2 vec = halfedgeVectorsInVertex[he];
      principalDir += -vec * vec / len * alpha;
    }

    vertexPrincipalCurvatureDirections[v] = principalDir / 4;
  }
}
void ExtrinsicGeometryInterface::requireVertexPrincipalCurvatureDirections() {
  vertexPrincipalCurvatureDirectionsQ.require();
}
void ExtrinsicGeometryInterface::unrequireVertexPrincipalCurvatureDirections() {
  vertexPrincipalCurvatureDirectionsQ.unrequire();
}


void ExtrinsicGeometryInterface::computeFacePrincipalCurvatureDirections() {
  edgeLengthsQ.ensureHave();
  halfedgeVectorsInFaceQ.ensureHave();
  edgeDihedralAnglesQ.ensureHave();

  facePrincipalCurvatureDirections = FaceData<Vector2>(mesh);

  for (Face f : mesh.faces()) {
    Vector2 principalDir{0.0, 0.0};
    for (Halfedge he : f.adjacentHalfedges()) {
      double len = edgeLengths[he.edge()];
      double alpha = edgeDihedralAngles[he.edge()];
      Vector2 vec = halfedgeVectorsInFace[he];
      principalDir += -vec * vec / len * alpha;
    }

    facePrincipalCurvatureDirections[f] = principalDir / 4;
  }
}

void ExtrinsicGeometryInterface::requireFacePrincipalCurvatureDirections() {
  facePrincipalCurvatureDirectionsQ.require();
}
void ExtrinsicGeometryInterface::unrequireFacePrincipalCurvatureDirections() {
  facePrincipalCurvatureDirectionsQ.unrequire();
}

void ExtrinsicGeometryInterface::computeFaceSizing() {
  edgeLengthsQ.ensureHave();
  edgeDihedralAnglesQ.ensureHave();
  faceAreasQ.ensureHave();
  halfedgeVectorsInFaceQ.ensureHave();
  Eigen::Matrix2d Sizing_mat;

  faceSizing = FaceData<double>(mesh, 0.0);
  for (Face f : mesh.faces()) {
    Sizing_mat = Eigen::Matrix2d::Zero();
    for (Halfedge he : f.adjacentHalfedges()) {
      double len = edgeLengths[he.edge()];
      double alpha = edgeDihedralAngles[he.edge()];
      Vector2 vec = halfedgeVectorsInFace[he];
      Eigen::Vector2d vec_eig{-vec.y, vec.x};
      vec_eig.normalize();
      // std::cout << "THe vector vec eig is " << vec_eig.transpose() << " \n";
      Sizing_mat += -0.5 * alpha * len * vec_eig * vec_eig.transpose();
    }
    Sizing_mat /= faceAreas[f];

    Sizing_mat = Sizing_mat.transpose() * Sizing_mat; // normalize by refine angle, so that the resulting sizing is in
                                                      // units of edge length rather than curvature
    // We take the bigges eigenvalue of the sizing matrix as the face sizing; this corresponds to the direction of
    // maximum curvature, and the value of the eigenvalue is related to the curvature in that direction.
    Eigen::EigenSolver<Eigen::Matrix2d> Solve;
    Eigen::Vector2d eigenvalues;
    Solve.compute(Sizing_mat);
    eigenvalues = Solve.eigenvalues().real();
    // std::cout << "THe eigenvalues are " << eigenvalues << "\n";
    faceSizing[f] = std::max(fabs(eigenvalues[0]), fabs(eigenvalues[1]));
  }
}
void ExtrinsicGeometryInterface::requireFaceSizing() { faceSizingQ.require(); }
void ExtrinsicGeometryInterface::unrequireFaceSizing() { faceSizingQ.unrequire(); }

void ExtrinsicGeometryInterface::computeVertexSizing() {

  faceSizingQ.ensureHave();
  faceAreasQ.ensureHave();
  vertexDualAreasQ.ensureHave();
  vertexSizing = VertexData<double>(mesh, 0.0);
  double dual_area;
  double sizing;
  for (Vertex v : mesh.vertices()) {
    sizing = 0.;
    dual_area = vertexDualAreas[v];

    for (Face f : v.adjacentFaces()) {
      // std::cout << "THe face sizing is " << faceSizing[f] << "\n";
      // std::cout << " THe face area is " << faceAreas[f] << "\n";
      sizing += faceAreas[f] * faceSizing[f] / 3;
    }
    vertexSizing[v] = sizing / dual_area;
  }
}
void ExtrinsicGeometryInterface::requireVertexSizing() { vertexSizingQ.require(); }
void ExtrinsicGeometryInterface::unrequireVertexSizing() { vertexSizingQ.unrequire(); }
} // namespace surface
} // namespace geometrycentral
