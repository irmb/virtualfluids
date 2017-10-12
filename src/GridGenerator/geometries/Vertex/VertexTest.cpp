#include "gmock/gmock.h"
#include "Vertex.cuh"

using namespace testing;

class VertexTest : public Test 
{
public:
    Vertex vec1;
    Vertex vec2;

    void SetUp() 
	{
        vec1.x = vec1.y = vec1.z = 4.0f;
        vec2.x = 0.0f; vec2.y = 2.4f; vec2.x = -1.3f;
    }
};

TEST_F(VertexTest, overloadMinusOperator)
{
    Vertex vec3;
    vec3.x = vec2.x - vec1.x;
    vec3.y = vec2.y - vec1.y;
    vec3.z = vec2.z - vec1.z;

    Vertex v4 = vec2 - vec1;
    ASSERT_THAT((double)v4.x, DoubleEq(vec3.x));
    ASSERT_THAT((double)v4.y, DoubleEq(vec3.y));
    ASSERT_THAT((double)v4.z, DoubleEq(vec3.z));
}

TEST_F(VertexTest, overloadPlusOperator)
{
    Vertex vec3;
    vec3.x = vec2.x + vec1.x;
    vec3.y = vec2.y + vec1.y;
    vec3.z = vec2.z + vec1.z;

    Vertex v4 = vec2 + vec1;
    ASSERT_THAT((double)v4.x, DoubleEq(vec3.x));
    ASSERT_THAT((double)v4.y, DoubleEq(vec3.y));
    ASSERT_THAT((double)v4.z, DoubleEq(vec3.z));
}

TEST_F(VertexTest, overloadTimesOperatorWithSkalarProduct)
{
    doubflo skalar = vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
    ASSERT_THAT((double)(vec1 * vec2), DoubleEq(skalar));
}

TEST_F(VertexTest, overloadTimesOperatorWithSkalarMultiplication)
{
    doubflo skalar = 1.0f / 3.0f;
    Vertex vec3;
    vec3.x = skalar * vec1.x;
    vec3.y = skalar * vec1.y;
    vec3.z = skalar * vec1.z;

    Vertex v4 = vec1 * skalar;

    ASSERT_THAT((double)v4.x, DoubleEq(vec3.x));
    ASSERT_THAT((double)v4.y, DoubleEq(vec3.y));
    ASSERT_THAT((double)v4.z, DoubleEq(vec3.z));
}

TEST_F(VertexTest, getMagnitudeFromVector)
{
    Vertex v;
    v.x = 4.0;
    v.y = 3.0;
    v.z = -1.0;

    float expected = (float)(std::sqrt(16.0 + 9.0 + 1.0));
    ASSERT_FLOAT_EQ(v.getMagnitude(), expected);
}

TEST_F(VertexTest, compareTwoVectors)
{
    Vertex v;
    v.x = vec1.x;
    v.y = vec1.y;
    v.z = vec1.z;
    ASSERT_THAT(v.isEqual(vec1), Eq(1));
}

TEST_F(VertexTest, checkEuclideanDistance)
{
    Vertex v = Vertex(3, 3, 3);

    ASSERT_FLOAT_EQ(v.getEuclideanDistanceTo(vec1), (float)sqrt(3));
}

TEST_F(VertexTest, checkEuclideanDistanceWithNullVector_ExpectNull)
{
    Vertex v1 = Vertex(0.0, 0.0, 0.0);
    Vertex v2 = Vertex(0.0, 0.0, 0.0);

    ASSERT_THAT((double)v1.getEuclideanDistanceTo(v2), DoubleEq(0.0));
}

TEST(VertexAngleTest, checkInnerAngleBetweenToVectors_ExpectRightAngle)
{
    Vertex v1 = Vertex(1.0, 4.0, -2.0);
    Vertex v2 = Vertex(-3.0, 3.0, 1);

    ASSERT_THAT((int)floor(v1.getInnerAngle(v2)), Eq(69));
}

TEST(VertexAngleTest, checkInnerAngleBetweenSameVectors_ExpectNull)
{
    Vertex v1 = Vertex(1.0, 4.0, -2.0);
    Vertex v2 = Vertex(1.0, 4.0, -2.0);

    ASSERT_THAT((int)floor(v1.getInnerAngle(v2)), Eq(0.0));
}

TEST(VertexAngleTest, checkInnerAngleBetweenNullVectors_ExpectNull)
{
    Vertex v1 = Vertex(0.0, 0.0, 0.0);
    Vertex v2 = Vertex(0.0, 0.0, 0.0);

    ASSERT_THAT((int)floor(v1.getInnerAngle(v2)), Eq(0.0));
}


TEST(VertexAngleTest, checkInnerAngleBetweenSecondNullVectors_ExpectNull)
{
    Vertex v1 = Vertex(1.0, 0.0, 0.0);
    Vertex v2 = Vertex(0.0, 0.0, 0.0);

    ASSERT_THAT((int)floor(v1.getInnerAngle(v2)), Eq(0.0));
}

TEST(VertexAngleTest, checkInnerAngleBetweenFirstNullVectors_ExpectNull)
{
    Vertex v1 = Vertex(0.0, 0.0, 0.0);
    Vertex v2 = Vertex(2.0, 0.0, 0.0);

    ASSERT_THAT((int)floor(v1.getInnerAngle(v2)), Eq(0.0));
}


TEST_F(VertexTest, crossProductBetweenTwoVectors)
{
    Vertex v1 = Vertex(-5.0, -5.0, 0.0);
    Vertex v2 = Vertex(5.0, 0.0, 10);

    Vertex crossProd = Vertex(-50.0, 50.0, 25.0);
    Vertex testCrossProduct = v1.crossProduct(v2);

    EXPECT_THAT((double)testCrossProduct.x, DoubleEq(crossProd.x));
    EXPECT_THAT((double)testCrossProduct.y, DoubleEq(crossProd.y));
    EXPECT_THAT((double)testCrossProduct.z, DoubleEq(crossProd.z));
}
