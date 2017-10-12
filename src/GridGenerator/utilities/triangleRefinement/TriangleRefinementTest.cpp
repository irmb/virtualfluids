#include <vector>
#include <memory>

#include "gmock/gmock.h"
#include <GridGenerator/geometries/Vertex/Vertex.cuh>
#include <GridGenerator/geometries/Triangle/Triangle.cuh>
#include "TriangleRefinement.h"
#include <GridGenerator/utilities/triangleNeighborFinder/TriangleNeighborFinder.h>
#include <GridGenerator/io/GridVTKWriter/GridVTKWriter.h>
#include <GridGenerator/utilities/Transformator/Transformator.h>

using namespace testing;



class TriangleRefinementTest : public Test {
public:
    std::shared_ptr<TriangleRefinement> refiner;
    std::vector<Triangle> triangles;

    void SetUp() {
		triangles.push_back(Triangle(Vertex((doubflo)0, (doubflo)1, (doubflo)1), Vertex((doubflo)0, (doubflo)1, (doubflo)0), Vertex((doubflo)0, (doubflo)0, (doubflo)1)));
		triangles.push_back(Triangle(Vertex((doubflo)0, (doubflo)1, (doubflo)0), Vertex((doubflo)0, (doubflo)0, (doubflo)0), Vertex((doubflo)0, (doubflo)0, (doubflo)1)));
		triangles.push_back(Triangle(Vertex((doubflo)1, (doubflo)1, (doubflo)1),Vertex((doubflo)1, (doubflo)1, (doubflo)0),Vertex((doubflo)0, (doubflo)1, (doubflo)1)));
		triangles.push_back(Triangle(Vertex((doubflo)1, (doubflo)1, (doubflo)0),Vertex((doubflo)0, (doubflo)1, (doubflo)0),Vertex((doubflo)0, (doubflo)1, (doubflo)1)));
		triangles.push_back(Triangle(Vertex((doubflo)1, (doubflo)0, (doubflo)1),Vertex((doubflo)1, (doubflo)0, (doubflo)0),Vertex((doubflo)1, (doubflo)1, (doubflo)1)));
		triangles.push_back(Triangle(Vertex((doubflo)1, (doubflo)0, (doubflo)0),Vertex((doubflo)1, (doubflo)1, (doubflo)0),Vertex((doubflo)1, (doubflo)1, (doubflo)1)));
		triangles.push_back(Triangle(Vertex((doubflo)0, (doubflo)0, (doubflo)1),Vertex((doubflo)0, (doubflo)0, (doubflo)0),Vertex((doubflo)1, (doubflo)0, (doubflo)1)));
		triangles.push_back(Triangle(Vertex((doubflo)0, (doubflo)0, (doubflo)0),Vertex((doubflo)1, (doubflo)0, (doubflo)0),Vertex((doubflo)1, (doubflo)0, (doubflo)1)));
		triangles.push_back(Triangle(Vertex((doubflo)1, (doubflo)1, (doubflo)1),Vertex((doubflo)0, (doubflo)1, (doubflo)1),Vertex((doubflo)1, (doubflo)0, (doubflo)1)));
		triangles.push_back(Triangle(Vertex((doubflo)0, (doubflo)1, (doubflo)1),Vertex((doubflo)0, (doubflo)0, (doubflo)1),Vertex((doubflo)1, (doubflo)0, (doubflo)1)));
		triangles.push_back(Triangle(Vertex((doubflo)0, (doubflo)1, (doubflo)0),Vertex((doubflo)1, (doubflo)1, (doubflo)0),Vertex((doubflo)0, (doubflo)0, (doubflo)0)));
		triangles.push_back(Triangle(Vertex((doubflo)1, (doubflo)1, (doubflo)0),Vertex((doubflo)1, (doubflo)0, (doubflo)0),Vertex((doubflo)0, (doubflo)0, (doubflo)0)));

        refiner = std::make_shared<TriangleRefinement>(&triangles);
    }
};

TEST_F(TriangleRefinementTest, refineOneTriangle_shouldCreateTwoMoreTriangle) {
    int oldSize = (int)triangles.size();
    refiner->refine(0);

    ASSERT_THAT(triangles.size(), Eq(oldSize + 2));
}


TEST_F(TriangleRefinementTest, refineOneTriangle_triangleShouldBeDeleteAfterRefine) {
    int refine = 0;
    Triangle oldTriangle = triangles[refine];
    refiner->refine(refine);

    ASSERT_FALSE(oldTriangle.getNumberOfCommonEdge(triangles[refine]) == 3);
}

TEST_F(TriangleRefinementTest, refineOneTriangle_newTrianglesAtTheEndMustHaveTwoCommonEdgesWithOld) {
    int refine = 0;
    Triangle oldTriangle = triangles[refine];
    refiner->refine(refine);

    Triangle firstNewTriangle = triangles[triangles.size() - 3];
    Triangle secondNewTriangle = triangles[triangles.size() - 4];
    EXPECT_THAT(oldTriangle.getNumberOfCommonEdge(firstNewTriangle), Eq(2));
    EXPECT_THAT(oldTriangle.getNumberOfCommonEdge(secondNewTriangle), Eq(2));
}

TEST_F(TriangleRefinementTest, refineQuadar_shouldDoubleTheTriangles_PlusTwoCauseOfFirstLoop) {
    int oldSize = (int)triangles.size();
    refiner->refineUntilMinDistance(1);

    ASSERT_THAT(triangles.size(), Eq(oldSize * 2 + 2));
}


TEST_F(TriangleRefinementTest, getHalfVertex){
    Vertex v1(10.0, 5.0, -2.0);
    Vertex v2(5.0, 0.0, 10);

    Vertex half = TriangleRefinement::getHalfVertex(v1, v2);


    EXPECT_THAT((double)half.x, DoubleEq(7.5));
    EXPECT_THAT((double)half.y, DoubleEq(2.5));
    EXPECT_THAT((double)half.z, DoubleEq(4));
}


TEST_F(TriangleRefinementTest, getLongestDistance){

    Triangle t = Triangle(Vertex(0, 0, 0), Vertex(10, 0, 0), Vertex(1, 1, 0), Vertex(0, 0, 0));

    doubflo d = TriangleRefinement::getLongestEdgeDistance(t);
    ASSERT_THAT((double)d, DoubleEq(10.0));
}

TEST_F(TriangleRefinementTest, getLongestDistanceEdge_ExpectEdgeZero){

    Triangle t = Triangle(Vertex(0, 0, 0), Vertex(10, 0, 0), Vertex(4, 4, 0), Vertex(0, 0, 0));

    int edge = TriangleRefinement::getEdgeWithLongestDistance(t);
    ASSERT_THAT(edge, Eq(0));
}

TEST_F(TriangleRefinementTest, getLongestDistanceEdge_ExpectEdgeOne){

    Triangle t = Triangle(Vertex(4, 4, 0), Vertex(10, 0, 0), Vertex(0, 0, 0), Vertex(0, 0, 0));

    int edge = TriangleRefinement::getEdgeWithLongestDistance(t);
    ASSERT_THAT(edge, Eq(1));
}

TEST_F(TriangleRefinementTest, getLongestDistanceEdge_ExpectEdgeTwo){

    Triangle t = Triangle(Vertex(0, 0, 0), Vertex(4, 4, 0), Vertex(10, 0, 0), Vertex(0, 0, 0));

    int edge = TriangleRefinement::getEdgeWithLongestDistance(t);
    ASSERT_THAT(edge, Eq(2));
}
