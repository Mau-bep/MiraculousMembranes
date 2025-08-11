#pragma once

#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/surface_mesh.h"
#include <unistd.h>
#include <Eigen/SparseCore>

namespace geometrycentral {
namespace surface {


class VertexPositionGeometry : public EmbeddedGeometryInterface {

  public:
    // Construct empty -- all positions initially set to the origin
    VertexPositionGeometry(SurfaceMesh& mesh_);

    // Construct from positions
    VertexPositionGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions);

    // Construct from positions (stored in an Eigen matrix)
    template <typename T>
    VertexPositionGeometry(SurfaceMesh& mesh_, const Eigen::MatrixBase<T>& vertexPositions);

    // Boring destructor
    virtual ~VertexPositionGeometry() {}

    // Construct a new geometry which is exactly the same as this one, on the
    // same mesh. This is a deep copy, no quantites are shared, etc. Require
    // counts/computed quantities are not copied.
    std::unique_ptr<VertexPositionGeometry> copy();

    // Construct a new geometry which is exactly the same as this one, on
    // another mesh. This is a deep copy, no quantites are shared, etc. Require
    // counts/computed quantities are not copied. The meshes must be in
    // correspondence (have the same connectivity).
    std::unique_ptr<VertexPositionGeometry> reinterpretTo(SurfaceMesh& targetMesh);


    // == Members

    // The actual input data which defines the geometry
    VertexData<Vector3> inputVertexPositions;

    // == Immediates
    double edgeLength(Edge e) const;
    double faceArea(Face f) const;
    double vertexDualArea(Vertex v) const;
    double cornerAngle(Corner c) const;
    double halfedgeCotanWeight(Halfedge he) const;
    double edgeCotanWeight(Edge e) const;
    Vector3 faceNormal(Face f) const;
    Vector3 halfedgeVector(Halfedge he) const;
    double edgeDihedralAngle(Edge e) const;
    double vertexMeanCurvature(Vertex v) const;
    double vertexGaussianCurvature(Vertex v) const;
    double vertexMinPrincipalCurvature(Vertex v) const;
    double vertexMaxPrincipalCurvature(Vertex v) const;

    // CHANGED: for DDG
    int eulerCharacteristic() const;
    double meanEdgeLength() const;
    
    double totalArea() const;
    double totalVolume() const;
    double cotan(Halfedge he) const;
    double barycentricDualArea(Vertex v) const;
    double angle(Corner c) const;
    double dihedralAngle(Halfedge he) const;
    Vector3 vertexNormalEquallyWeighted(Vertex v) const;
    Vector3 vertexNormalAngleWeighted(Vertex v) const;
    Vector3 vertexNormalSphereInscribed(Vertex v) const;
    Vector3 vertexNormalAreaWeighted(Vertex v) const;
    Vector3 vertexNormalGaussianCurvature(Vertex v) const;
    Vector3 vertexNormalMeanCurvature(Vertex v) const;
    
    // Courtesy of Mem3DG
    Vector3 computeHalfedgeMeanCurvatureVector(Halfedge he) const ;
    Vector3 computeHalfedgeGaussianCurvatureVector(Halfedge he) const ;
    Vector3 dihedralAngleGradient(Halfedge he, Vertex v) const ;



    Vector3 Angle_grad(Corner c, Vertex j, Vector3 Normal);
    double angleDefect(Vertex v) const;
    double totalAngleDefect() const;
    double scalarMeanCurvature(Vertex v) const;
    double circumcentricDualArea(Vertex v) const;
    std::pair<double, double> principalCurvatures(Vertex v) const;


    
    SparseMatrix<double> laplaceMatrix() const;
    SparseMatrix<double> laplaceMatrix2() const;
    SparseMatrix<double> laplaceMatrix3d() const ;

    SparseMatrix<double> uniformlaplacianMatrix() const;

    SparseMatrix<double> NormalFlowMat() const;
    SparseMatrix<double> NormalFlowMatx() const;
    SparseMatrix<double> NormalFlowMaty() const;
    SparseMatrix<double> NormalFlowMatz() const;
    
    SparseMatrix<double> GaussFlowMat() const;
    SparseMatrix<double> GaussFlowMat3d(VertexData<double> Scalar_MC) const;
    SparseMatrix<double> MeanHFlowMat3d(VertexData<double> Scalar_MC,Vector<Vector3> face_Normals) const;
    SparseMatrix<double> Schlafi3d(VertexData<double> Scalar_MC,Vector<Vector3> face_Normals) const ;

    Vector3 dihedralAngleGradient(Halfedge he, Vertex v,Vector<Vector3> face_Normals) const;

    SparseMatrix<double> LocalQWillmore(const Halfedge he, const Vector3 x0, const Vector3 x1, const Vector3 x2, const Vector3 x3) const;

    SparseMatrix<double> FullQWillmore() const;
    
    SparseMatrix<double> FullKWillmore() const;
    





    SparseMatrix<double> massMatrix() const;
    SparseMatrix<double> massMatrix3d() const;
    
    SparseMatrix<std::complex<double>> complexLaplaceMatrix() const;
    Vector3 centerOfMass() const;
    void normalize(const Vector3& origin = {0, 0, 0}, bool rescale = false);
    void rescale( double scale_factor);
    SparseMatrix<double> buildExteriorDerivative1Form() const;
    SparseMatrix<double> buildHodgeStar0Form() const;
    SparseMatrix<double> buildHodgeStar1Form() const;
    SparseMatrix<double> buildHodgeStar2Form() const;
    SparseMatrix<double> buildExteriorDerivative0Form() const;
    SparseMatrix<double> buildExteriorDerivative0Form3N() const;


    // We will be adding some cute functions for the Hessian calculations 

    Eigen::Matrix3d Cross_product_matrix(Eigen::Vector3d v) const;
    double Triangle_area(Eigen::Vector<double,9> Positions) const;
    Eigen::Vector<double,9> gradient_triangle_area( Eigen::Vector<double,9> Positions) const; 
    Eigen::Matrix<double, 9, 9> hessian_triangle_area( Eigen::Vector<double,9> Positions) const;


    double Edge_length(Eigen::Vector<double,6> Positions) const;
    Eigen::Vector<double, 6> gradient_edge_length(Eigen::Vector<double,6> Positions) const;
    Eigen::Matrix<double,6,6> hessian_edge_length(Eigen::Vector<double,6> Positions) const;

    double Ej_edge_regular(Eigen::Vector<double,12> Positions) const;
    double Ej_edge_regular(Eigen::Vector<double,12> Positions, Eigen::Vector<double,5> Edge_lengths) const;
    Eigen::Vector<double,12> gradient_edge_regular(Eigen::Vector<double,12> Positions) const;
    Eigen::Vector<double,12> gradient_edge_regular(Eigen::Vector<double,12> Positions, Eigen::Vector<double,5> Edge_lengths) const;
    
    Eigen::Vector<double,12> gradient_edge_regular_I(Eigen::Vector<double,12> Positions, Eigen::Vector<double,5> Edge_lengths) const;

    Eigen::Matrix<double,12,12> hessian_edge_regular(Eigen::Vector<double,12> Positions) const;
    Eigen::Matrix<double,12,12> hessian_edge_regular(Eigen::Vector<double,12> Positions, Eigen::Vector<double,5> Edge_lengths) const;

    double Dihedral_angle(Eigen::Vector<double, 12> Positions) const;
    Eigen::Vector<double, 12> gradient_dihedral_angle(Eigen::Vector<double, 12> Positions) const;
    Eigen::Matrix<double, 12, 12> hessian_dihedral_angle(Eigen::Vector<double, 12> Positions) const;


    double Cotan(Eigen::Vector<double,9> Positions) const;

    double Angle(Eigen::Vector<double,9> Positions) const;
    Eigen::Vector<double, 9 > gradient_angle(Eigen::Vector<double, 9> Positions) const;

    double Cotan_weight(Eigen::Vector<double,12> Positions) const;
    Eigen::Vector<double, 12> gradient_cotan_weight(Eigen::Vector<double, 12> Positions) const;
    Eigen::Matrix<double, 12, 12> hessian_cotan_weight(Eigen::Vector<double, 12> Positions) const;

    double Volume(Eigen::Vector<double, 9> Positions) const;
    Eigen::Vector<double, 9> gradient_volume(Eigen::Vector<double, 9> Positions) const;
    Eigen::Matrix<double, 9, 9> hessian_volume(Eigen::Vector<double, 9> Positions) const;



  protected:
    // Override the compute vertex positions method for embedded geometry
    virtual void computeVertexPositions() override;

    double vertexPrincipalCurvature(int whichCurvature, Vertex v) const;

  private:
};

} // namespace surface
} // namespace geometrycentral

#include "geometrycentral/surface/vertex_position_geometry.ipp"
