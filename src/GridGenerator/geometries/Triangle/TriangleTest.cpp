#include "gmock/gmock.h"
#include <GridGenerator/utilities/Transformator/Transformator.h>
#include "TriangleException.h"

#include "Triangle.h"


using namespace testing;

TEST(TriangleTest, callingStandardConstructor_createsCorrectTriangle)
{
	Triangle t;

    bool isV1Null = (t.v1 == Vertex(0, 0, 0));
    bool isV2Null = (t.v2 == Vertex(0, 0, 0));
    bool isV3Null = (t.v3 == Vertex(0, 0, 0));
    bool isNormalNull = (t.normal == Vertex(0, 0, 0));

    EXPECT_TRUE(isV1Null);
    EXPECT_TRUE(isV2Null);
    EXPECT_TRUE(isV3Null);
    EXPECT_TRUE(isNormalNull);
}

TEST(TriangleTest, changeTriangleWithSetMethod)
{
	Triangle t;

	real v1x = 3.0f;
	real v1y = 3.0f;
	real v1z = 3.0f;

	real v2x = 1.0f;
	real v2y = 1.0f;
	real v2z = 3.0f;

	real v3x = -1.0f;
	real v3y = 1.0f;
	real v3z = 3.0f;

	t.set(Vertex(v1x,v1y,v1z), Vertex(v2x, v2y, v2z), Vertex(v3x, v3y, v3z));


    bool isV1Set = (t.v1 == Vertex(v1x, v1y, v1z));
    bool isV2Set = (t.v2 == Vertex(v2x, v2y, v2z));
    bool isV3Set = (t.v3 == Vertex(v3x, v3y, v3z));
    bool isNormalSet = (t.normal == Vertex(0, 0, -1.0));

    EXPECT_TRUE(isV1Set);
    EXPECT_TRUE(isV2Set);
    EXPECT_TRUE(isV3Set);
    EXPECT_TRUE(isNormalSet);
}


TEST(TriangleTest, getClosestPointsOnEdgesFromTriangle)
{
    Triangle  t = Triangle(Vertex(0, 0, 0), Vertex(10, 0, 0), Vertex(0, 10, 0), Vertex(0, 0, 1));
    Vertex P4 = Vertex(0, 1, 5);

    Vertex p[3];
    t.getClosestPointsOnEdges(p, P4);

    bool isP1OnEdge = (p[0] == Vertex(0,0,0));
    bool isP2OnEdge = (p[1] == Vertex(4.5,5.5,0.0));
    bool isP3OnEdge = (p[2] == Vertex(0.0,1.0,0.0));

    EXPECT_TRUE(isP1OnEdge);
    EXPECT_TRUE(isP2OnEdge);
    EXPECT_TRUE(isP3OnEdge);
}

//TEST(TriangleTest, PointIsNotWithinTriangle)
//{
//    Triangle t = Triangle(Vertex(0, 0, 0), Vertex(10, 0, 0), Vertex(0, 10, 0), Vertex(0, 0, 1));
//    Vertex P4 = Vertex(1, 1, -2);
//
//    EXPECT_FALSE(t.isNotNextToFace(P4));
//}

TEST(TriangleTest, calculatePerpendicularPoint)
{
    Triangle t = Triangle(Vertex(0, 0, 0), Vertex(10, 0, 0), Vertex(0, 10, 0), Vertex(0, 0, 1));
    Vertex  P4 = Vertex(8, 5, -3);

    Vertex distPbances = t.getPerpedicularPointFrom(P4);

    bool isVertexEqual = distPbances == Vertex(8.0, 5.0, 0.0);
    EXPECT_TRUE(isVertexEqual);
}

TEST(TriangleTest, pointintersectTriangle_directionIntoTriangle_ExpectReturn1)
{
    Vertex normal = Vertex(-1.0f, 0.0f, 0.0f);
    Triangle t = Triangle(Vertex(3.5, 0, 0), Vertex(3.5, 0, 9), Vertex(3.5, 9, 9), normal);
    Vertex v = Vertex(3, 5, 5);
    Vertex direction = Vertex(1.0f, 0, 0);
    Vertex intersect;

    real q;
    int err = t.getTriangleIntersection(v, direction, intersect, q);
    ASSERT_THAT(err, Eq(0));
}


TEST(TriangleTest, pointintersectTriangle_directionAgainstTriangle_ExpectReturn0)
{
    Vertex normal = Vertex(-1.0f, 0.0f, 0.0f);
    Triangle t = Triangle(Vertex(3.5, 0, 0), Vertex(3.5, 0, 9), Vertex(3.5, 9, 9), normal);
    Vertex v = Vertex(3, 5, 5);
    Vertex direction = Vertex(-1.0f, 0, 0);
    Vertex intersect;

    real q;
    int err = t.getTriangleIntersection(v, direction, intersect, q);
    ASSERT_THAT(err, Eq(1));
}


TEST(TriangleTest, getHalfAngleBetweenTwoEqualTriangles_ExpectNullAngle)
{
    Triangle t1 = Triangle(Vertex(0, 0, 0), Vertex(10, 0, 0), Vertex(0, 10, 0), Vertex(0, 0, 1));
    Triangle t2 = Triangle(Vertex(0, 0, 0), Vertex(10, 0, 0), Vertex(0, 10, 0), Vertex(0, 0, 1));

	real alpha = t1.getHalfAngleBetweenToAdjacentTriangle(t2);
    ASSERT_THAT(alpha, RealEq(0.0));
}

TEST(TriangleTest, checkSTLwith90degreeOutpointingNormal) 
{
	Triangle t1 = Triangle(Vertex(40.0f, 20.0f, 20.0f), Vertex(40.0f, 20.0f, 0.0f), Vertex(60.0f, 20.0f, 20.0f));
	Triangle t2 = Triangle(Vertex(60.0f, 20.0f, 20.0f), Vertex(40.0f, 40.0f, 20.0f), Vertex(40.0f, 20.0f, 20.0f));

	real alpha = t1.getHalfAngleBetweenToAdjacentTriangle(t2);
    ASSERT_THAT(alpha, RealEq(90.0 / 2.0));
}

TEST(TriangleTest, checkSTLwith90degreeInpointingNormal) 
{
	Triangle t1 = Triangle(Vertex(40.0f, 20.0f, 20.0f), Vertex(40.0f, 20.0f, 40.0f), Vertex(60.0f, 20.0f, 20.0f));
	Triangle t2 = Triangle(Vertex(60.0f, 20.0f, 20.0f), Vertex(40.0f, 40.0f, 20.0f), Vertex(40.0f, 20.0f, 20.0f));

	real alpha = t1.getHalfAngleBetweenToAdjacentTriangle(t2);
    ASSERT_THAT(alpha, RealEq(270.0 / 2.0));
}

TEST(TriangleTest, checkSTLwith180degreeOutpointingNormal)
{
	Triangle t1 = Triangle(Vertex(40.0f, 20.0f, 20.0f), Vertex(40.0f, 0.0f, 20.0f), Vertex(60.0f, 20.0f, 20.0f));
	Triangle t2 = Triangle(Vertex(60.0f, 20.0f, 20.0f), Vertex(40.0f, 40.0f, 20.0f), Vertex(40.0f, 20.0f, 20.0f));

	real alpha = t1.getHalfAngleBetweenToAdjacentTriangle(t2);
    ASSERT_THAT(alpha, RealEq(180.0 / 2.0));
}

TEST(TriangleTest, checkSTLwithSmallDegreeOutpointingNormal) 
{
	Triangle t2 = Triangle(Vertex(40.0f, 40.0f, 18.0f), Vertex(60.0f, 20.0f, 20.0f), Vertex(40.0f, 20.0f, 20.0f));
	Triangle t1 = Triangle(Vertex(60.0f, 20.0f, 20.0f), Vertex(40.0f, 40.0f, 20.0f), Vertex(40.0f, 20.0f, 20.0f));

	real alpha = t1.getHalfAngleBetweenToAdjacentTriangle(t2);
    ASSERT_TRUE(alpha < 30.0);
}

TEST(TriangleTest, checkSTLwithBigDegreeInpointingNormal) 
{
	Triangle t1 = Triangle(Vertex(40.0f, 40.0f, 20.0f), Vertex(60.0f, 20.0f, 20.0f), Vertex(40.0f, 20.0f, 20.0f));
	Triangle t2 = Triangle(Vertex(60.0f, 20.0f, 20.0f), Vertex(40.0f, 40.0f, 18.0f), Vertex(40.0f, 20.0f, 20.0f));

    real alpha = t1.getHalfAngleBetweenToAdjacentTriangle(t2);
    ASSERT_TRUE(alpha > 330.0 / 2);
}

TEST(TriangleTest, getCommonEdgesFromTheSameTriangle_ExpectThreeEdges)
{
    Triangle t1 = Triangle(Vertex(0.0f, 0.0f, 0.0f), Vertex(10.0f, 0.0f, 0.0f), Vertex(0.0f, 10.0f, 0.0f));
    Triangle t2 = Triangle(Vertex(0.0f, 0.0f, 0.0f), Vertex(10.0f, 0.0f, 0.0f), Vertex(0.0f, 10.0f, 0.0f));
    int num = t1.getNumberOfCommonEdge(t2);
    ASSERT_THAT(num, Eq(3));
}

TEST(TriangleTest, getCommonEdgeFromDifferentTriangles_ShouldReturnMinusOne)
{
	Triangle t1 = Triangle(Vertex(0.0f, 0.0f, 0.0f), Vertex(10.0f, 0.0f, 0.0f), Vertex(0.0f, 10.0f, 0.0f));
	Triangle t2 = Triangle(Vertex(10.0f, 100.0f, 100.0f), Vertex(102.0f, 122.0f, 100.0f), Vertex(3000.0f, 1000.0f, 100.0f));


	int num = t1.getCommonEdge(t2);
	ASSERT_THAT(num, -1);
}
