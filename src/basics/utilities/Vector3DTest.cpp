#include "Vector3D.h"
#include "gmock/gmock.h"

#include <cmath>

using namespace testing;

class Vector3DTest : public Test
{
public:
    Vector3D vec1;
    Vector3D vec2;

    void SetUp() override
    {
        vec1[0] = vec1[1] = vec1[2] = 4.0f;
        vec2[0]                     = 0.0f;
        vec2[0]                     = 2.4f;
        vec2[0]                     = -1.3f;
    }
};

TEST_F(Vector3DTest, overloadMinusOperator)
{
    Vector3D vec3;
    vec3[0] = vec2[0] - vec1[0];
    vec3[1] = vec2[1] - vec1[1];
    vec3[2] = vec2[2] - vec1[2];

    Vector3D v4 = vec2 - vec1;
    ASSERT_THAT((double)v4[0], DoubleEq(vec3[0]));
    ASSERT_THAT((double)v4[1], DoubleEq(vec3[1]));
    ASSERT_THAT((double)v4[2], DoubleEq(vec3[2]));
}

TEST_F(Vector3DTest, overloadPlusOperator)
{
    Vector3D vec3;
    vec3[0] = vec2[0] + vec1[0];
    vec3[1] = vec2[1] + vec1[1];
    vec3[2] = vec2[2] + vec1[2];

    Vector3D v4 = vec2 + vec1;
    ASSERT_THAT((double)v4[0], DoubleEq(vec3[0]));
    ASSERT_THAT((double)v4[1], DoubleEq(vec3[1]));
    ASSERT_THAT((double)v4[2], DoubleEq(vec3[2]));
}

TEST_F(Vector3DTest, overloadTimesOperatorWithSkalarProduct)
{
    double skalar = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
    ASSERT_THAT(vec1.Dot(vec2), DoubleEq(skalar));
}

TEST_F(Vector3DTest, overloadTimesOperatorWithSkalarMultiplication)
{
    double skalar = 1.0f / 3.0f;
    Vector3D vec3;
    vec3[0] = skalar * vec1[0];
    vec3[1] = skalar * vec1[1];
    vec3[2] = skalar * vec1[2];

    Vector3D v4 = vec1 * skalar;

    ASSERT_THAT((double)v4[0], DoubleEq(vec3[0]));
    ASSERT_THAT((double)v4[1], DoubleEq(vec3[1]));
    ASSERT_THAT((double)v4[2], DoubleEq(vec3[2]));
}

TEST_F(Vector3DTest, getLengthFromVector)
{
    Vector3D v;
    v[0] = 4.0;
    v[1] = 3.0;
    v[2] = -1.0;

    double expected = std::sqrt(16.0 + 9.0 + 1.0);
    ASSERT_THAT(v.Length(), testing::DoubleEq(expected));
}

TEST_F(Vector3DTest, compareTwoVectors)
{
    Vector3D v;
    v[0] = vec1[0];
    v[1] = vec1[1];
    v[2] = vec1[2];
    ASSERT_TRUE(v == vec1);
}
//
// TEST_F(Vector3DTest, checkEuclideanDistance)
//{
//    Vector3D v = Vector3D(3, 3, 3);
//
//    ASSERT_FLOAT_EQ(v.getEuclideanDistanceTo(vec1), (float)sqrt(3));
//}
//
// TEST_F(Vector3DTest, checkEuclideanDistanceWithNullVector_ExpectNull)
//{
//    Vector3D v1 = Vector3D(0.0, 0.0, 0.0);
//    Vector3D v2 = Vector3D(0.0, 0.0, 0.0);
//
//    ASSERT_THAT((double)v1.getEuclideanDistanceTo(v2), DoubleEq(0.0));
//}
//
// TEST(Vector3DAngleTest, checkInnerAngleBetweenToVectors_ExpectRightAngle)
//{
//    Vector3D v1 = Vector3D(1.0, 4.0, -2.0);
//    Vector3D v2 = Vector3D(-3.0, 3.0, 1);
//
//    ASSERT_THAT((int)floor(v1.getInnerAngle(v2)), Eq(69));
//}
//
// TEST(Vector3DAngleTest, checkInnerAngleBetweenSameVectors_ExpectNull)
//{
//    Vector3D v1 = Vector3D(1.0, 4.0, -2.0);
//    Vector3D v2 = Vector3D(1.0, 4.0, -2.0);
//
//    ASSERT_THAT((int)floor(v1.getInnerAngle(v2)), Eq(0.0));
//}
//
// TEST(Vector3DAngleTest, checkInnerAngleBetweenNullVectors_ExpectNull)
//{
//    Vector3D v1 = Vector3D(0.0, 0.0, 0.0);
//    Vector3D v2 = Vector3D(0.0, 0.0, 0.0);
//
//    ASSERT_THAT((int)floor(v1.getInnerAngle(v2)), Eq(0.0));
//}
//
//
// TEST(Vector3DAngleTest, checkInnerAngleBetweenSecondNullVectors_ExpectNull)
//{
//    Vector3D v1 = Vector3D(1.0, 0.0, 0.0);
//    Vector3D v2 = Vector3D(0.0, 0.0, 0.0);
//
//    ASSERT_THAT((int)floor(v1.getInnerAngle(v2)), Eq(0.0));
//}
//
// TEST(Vector3DAngleTest, checkInnerAngleBetweenFirstNullVectors_ExpectNull)
//{
//    Vector3D v1 = Vector3D(0.0, 0.0, 0.0);
//    Vector3D v2 = Vector3D(2.0, 0.0, 0.0);
//
//    ASSERT_THAT((int)floor(v1.getInnerAngle(v2)), Eq(0.0));
//}

TEST_F(Vector3DTest, crossProductBetweenTwoVectors)
{
    Vector3D v1 = Vector3D(-5.0, -5.0, 0.0);
    Vector3D v2 = Vector3D(5.0, 0.0, 10);

    Vector3D crossProd        = Vector3D(-50.0, 50.0, 25.0);
    Vector3D testCrossProduct = v1.Cross(v2);

    EXPECT_THAT(testCrossProduct[0], DoubleEq(crossProd[0]));
    EXPECT_THAT(testCrossProduct[1], DoubleEq(crossProd[1]));
    EXPECT_THAT(testCrossProduct[2], DoubleEq(crossProd[2]));
}
