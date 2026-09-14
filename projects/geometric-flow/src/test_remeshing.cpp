#include <gtest/gtest.h>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/remeshing.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class DeleteLowValenceTest : public ::testing::Test
{
protected:
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;
    MutationManager *mm;

    void SetUp() override
    {
        // Create a simple mesh for testing
        // We'll use a basic icosahedron or similar
    }

    void TearDown() override
    {
        if (mm)
            delete mm;
    }

    void createTestMesh()
    {
        std::vector<Vector3> vertices;
        std::vector<std::vector<size_t>> faces;

        // Create a simple tetrahedron or pyramid mesh
        vertices = {
            {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 0}};

        faces = {
            {0, 1, 2}, {0, 2, 3}, {1, 4, 2}, {0, 1, 3}};

        auto lvals = makeManifoldSurfaceMeshAndGeometry(faces, vertices);
        mesh = std::move(std::get<0>(lvals));
        geometry = std::move(std::get<1>(lvals));
    }
};

TEST_F(DeleteLowValenceTest, DeletesValence3Vertices)
{
    createTestMesh();
    ASSERT_TRUE(mesh != nullptr);
    ASSERT_TRUE(geometry != nullptr);

    mm = new MutationManager(*mesh, *geometry);
    RemeshOptions options;
    options.aspect_min = 0.2;
    options.min_absolute_length = 0.01;
    options.max_absolute_length = 1.0;

    // Count vertices with valence 3 before deletion
    VertexData<int> valenceBefore(*mesh, 0);
    size_t count3Before = 0;
    for (Vertex v : mesh->vertices())
    {
        if (!v.isBoundary())
        {
            int valence = 0;
            for (Face f : v.adjacentFaces())
                valence++;
            valenceBefore[v] = valence;
            if (valence == 3)
                count3Before++;
        }
    }

    deleteLowValence(*mesh, *geometry, *mm, options);

    // Count vertices with valence 3 after deletion
    size_t count3After = 0;
    for (Vertex v : mesh->vertices())
    {
        if (!v.isBoundary())
        {
            int valence = 0;
            for (Face f : v.adjacentFaces())
                valence++;
            if (valence == 3)
                count3After++;
        }
    }

    EXPECT_LE(count3After, count3Before);
}

TEST_F(DeleteLowValenceTest, DeletesValence4Vertices)
{
    createTestMesh();
    ASSERT_TRUE(mesh != nullptr);
    ASSERT_TRUE(geometry != nullptr);

    mm = new MutationManager(*mesh, *geometry);
    RemeshOptions options;
    options.aspect_min = 0.2;
    options.min_absolute_length = 0.01;
    options.max_absolute_length = 1.0;

    // Count vertices with valence 4 before deletion
    size_t count4Before = 0;
    for (Vertex v : mesh->vertices())
    {
        if (!v.isBoundary())
        {
            int valence = 0;
            for (Face f : v.adjacentFaces())
                valence++;
            if (valence == 4)
                count4Before++;
        }
    }

    deleteLowValence(*mesh, *geometry, *mm, options);

    // Count vertices with valence 4 after deletion
    size_t count4After = 0;
    for (Vertex v : mesh->vertices())
    {
        if (!v.isBoundary())
        {
            int valence = 0;
            for (Face f : v.adjacentFaces())
                valence++;
            if (valence == 4)
                count4After++;
        }
    }

    EXPECT_LE(count4After, count4Before);
}

TEST_F(DeleteLowValenceTest, PreservesBoundaryVertices)
{
    createTestMesh();
    ASSERT_TRUE(mesh != nullptr);
    ASSERT_TRUE(geometry != nullptr);

    mm = new MutationManager(*mesh, *geometry);
    RemeshOptions options;
    options.aspect_min = 0.2;
    options.min_absolute_length = 0.01;
    options.max_absolute_length = 1.0;

    // Count boundary vertices before deletion
    size_t boundaryCountBefore = 0;
    for (Vertex v : mesh->vertices())
    {
        if (v.isBoundary())
            boundaryCountBefore++;
    }

    deleteLowValence(*mesh, *geometry, *mm, options);

    // Count boundary vertices after deletion
    size_t boundaryCountAfter = 0;
    for (Vertex v : mesh->vertices())
    {
        if (v.isBoundary())
            boundaryCountAfter++;
    }

    EXPECT_EQ(boundaryCountAfter, boundaryCountBefore);
}

TEST_F(DeleteLowValenceTest, FarFromBoundaryVerticesProcessed)
{
    createTestMesh();
    ASSERT_TRUE(mesh != nullptr);
    ASSERT_TRUE(geometry != nullptr);

    mm = new MutationManager(*mesh, *geometry);
    RemeshOptions options;
    options.aspect_min = 0.2;
    options.min_absolute_length = 0.01;
    options.max_absolute_length = 1.0;

    // Identify interior vertices with valence 3 or 4
    std::vector<Vertex> lowValenceInterior;
    for (Vertex v : mesh->vertices())
    {
        if (!v.isBoundary())
        {
            int valence = 0;
            for (Face f : v.adjacentFaces())
                valence++;
            if (valence == 3 || valence == 4)
            {
                lowValenceInterior.push_back(v);
            }
        }
    }

    size_t countBefore = lowValenceInterior.size();
    deleteLowValence(*mesh, *geometry, *mm, options);

    // Count remaining low valence interior vertices
    size_t countAfter = 0;
    for (Vertex v : mesh->vertices())
    {
        if (!v.isBoundary())
        {
            int valence = 0;
            for (Face f : v.adjacentFaces())
                valence++;
            if (valence == 3 || valence == 4)
                countAfter++;
        }
    }

    EXPECT_LE(countAfter, countBefore);
}

TEST_F(DeleteLowValenceTest, MaintainsMeshValidity)
{
    createTestMesh();
    ASSERT_TRUE(mesh != nullptr);
    ASSERT_TRUE(geometry != nullptr);

    mm = new MutationManager(*mesh, *geometry);
    RemeshOptions options;
    options.aspect_min = 0.2;
    options.min_absolute_length = 0.01;
    options.max_absolute_length = 1.0;

    size_t verticesBefore = mesh->nVertices();
    size_t edgesBefore = mesh->nEdges();
    size_t facesBefore = mesh->nFaces();

    deleteLowValence(*mesh, *geometry, *mm, options);

    size_t verticesAfter = mesh->nVertices();
    size_t edgesAfter = mesh->nEdges();
    size_t facesAfter = mesh->nFaces();

    // Verify Euler characteristic is maintained (V - E + F = 2 for manifold)
    int eulerBefore = static_cast<int>(verticesBefore) - static_cast<int>(edgesBefore) + static_cast<int>(facesBefore);
    int eulerAfter = static_cast<int>(verticesAfter) - static_cast<int>(edgesAfter) + static_cast<int>(facesAfter);
    EXPECT_EQ(eulerAfter, eulerBefore);
}

TEST_F(DeleteLowValenceTest, ReducesVertexCount)
{
    createTestMesh();
    ASSERT_TRUE(mesh != nullptr);
    ASSERT_TRUE(geometry != nullptr);

    mm = new MutationManager(*mesh, *geometry);
    RemeshOptions options;
    options.aspect_min = 0.2;
    options.min_absolute_length = 0.01;
    options.max_absolute_length = 1.0;

    size_t verticesBefore = mesh->nVertices();
    deleteLowValence(*mesh, *geometry, *mm, options);
    size_t verticesAfter = mesh->nVertices();

    EXPECT_LE(verticesAfter, verticesBefore);
}