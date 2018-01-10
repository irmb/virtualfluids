#include <vector>
#include <memory>

#include "gmock/gmock.h"

#include <GridGenerator/geometries/Triangle/Triangle.cuh>
#include "TriangleNeighborFinder.h"
#include <GridGenerator/utilities/Transformator/Transformator.h>

using namespace testing;


class TriangleNeighborFinderTest : public Test {
public:
    std::shared_ptr<TriangleNeighborFinder> finder;
    std::vector<Triangle> triangles;
    IntegerPtr2D indices;

    void SetUp() {
		triangles.push_back(Triangle(Vertex((real)0, (real)1, (real)1), Vertex((real)0, (real)1, (real)0), Vertex((real)0, (real)0, (real)1)));
		triangles.push_back(Triangle(Vertex((real)0, (real)1, (real)0), Vertex((real)0, (real)0, (real)0), Vertex((real)0, (real)0, (real)1)));
		triangles.push_back(Triangle(Vertex((real)1, (real)1, (real)1), Vertex((real)1, (real)1, (real)0), Vertex((real)0, (real)1, (real)1)));
		triangles.push_back(Triangle(Vertex((real)1, (real)1, (real)0), Vertex((real)0, (real)1, (real)0), Vertex((real)0, (real)1, (real)1)));
		triangles.push_back(Triangle(Vertex((real)1, (real)0, (real)1), Vertex((real)1, (real)0, (real)0), Vertex((real)1, (real)1, (real)1)));
		triangles.push_back(Triangle(Vertex((real)1, (real)0, (real)0), Vertex((real)1, (real)1, (real)0), Vertex((real)1, (real)1, (real)1)));
		triangles.push_back(Triangle(Vertex((real)0, (real)0, (real)1), Vertex((real)0, (real)0, (real)0), Vertex((real)1, (real)0, (real)1)));
		triangles.push_back(Triangle(Vertex((real)0, (real)0, (real)0), Vertex((real)1, (real)0, (real)0), Vertex((real)1, (real)0, (real)1)));
		triangles.push_back(Triangle(Vertex((real)1, (real)1, (real)1), Vertex((real)0, (real)1, (real)1), Vertex((real)1, (real)0, (real)1)));
		triangles.push_back(Triangle(Vertex((real)0, (real)1, (real)1), Vertex((real)0, (real)0, (real)1), Vertex((real)1, (real)0, (real)1)));
		triangles.push_back(Triangle(Vertex((real)0, (real)1, (real)0), Vertex((real)1, (real)1, (real)0), Vertex((real)0, (real)0, (real)0)));
		triangles.push_back(Triangle(Vertex((real)1, (real)1, (real)0), Vertex((real)1, (real)0, (real)0), Vertex((real)0, (real)0, (real)0)));

        finder = std::make_shared<TriangleNeighborFinder>(&triangles[0], (int)triangles.size());
        indices.size = (int) triangles.size();
        indices.DIM = DIMENSION;
        indices.ptr = new int[indices.DIM * triangles.size()];
    }

    void TearDown() {
        delete[] indices.ptr;
    }
};

TEST_F(TriangleNeighborFinderTest, getNeighboursFromTriangle) {
    unsigned int triangleID = 0;
    finder->fillWithNeighborIndices(&indices, &triangles[0]);

    EXPECT_THAT(indices.ptr[triangleID * indices.DIM + 0], Eq(3));
    EXPECT_THAT(indices.ptr[triangleID * indices.DIM + 1], Eq(1));
    EXPECT_THAT(indices.ptr[triangleID * indices.DIM + 2], Eq(9));
}

