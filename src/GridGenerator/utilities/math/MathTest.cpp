#include "gmock/gmock.h"
#include "Math.h"

using namespace testing;

class MathTest : public Test {
public:
	void SetUp() {

	}
};

TEST_F(MathTest, getDecimalPart) {
    ASSERT_THAT(vf::Math::getDecimalPart(2.23232), RealEq(0.23232));
}

TEST_F(MathTest, compareBigValues) {
	ASSERT_TRUE(vf::Math::equal(10000000.0f, 10000001.0f));
	ASSERT_TRUE(vf::Math::equal(10000001.0f, 10000000.0f));
	ASSERT_FALSE(vf::Math::equal(100000.0f, 100001.0f));
	ASSERT_FALSE(vf::Math::equal(100001.0f, 100000.0f));
}

TEST_F(MathTest, compareSmallValues) {
	ASSERT_TRUE(vf::Math::equal(-10000000.0f, -10000001.0f));
	ASSERT_TRUE(vf::Math::equal(-10000001.0f, -10000000.0f));
	ASSERT_FALSE(vf::Math::equal(-100000.0f, -100001.0f));
	ASSERT_FALSE(vf::Math::equal(-100001.0f, -100000.0f));
}

TEST_F(MathTest, compareNumbersAroundOne) {
	ASSERT_TRUE(vf::Math::equal(1.0000001f, 1.0000002f));
	ASSERT_TRUE(vf::Math::equal(1.0000002f, 1.0000001f));
	ASSERT_FALSE(vf::Math::equal(1.0002f, 1.0001f));
	ASSERT_FALSE(vf::Math::equal(1.0001f, 1.0002f));
}

TEST_F(MathTest, compareNumbersAroundMinusOne) {
	ASSERT_TRUE(vf::Math::equal(-1.0000001f, -1.0000002f));
	ASSERT_TRUE(vf::Math::equal(-1.0000002f, -1.0000001f));
	ASSERT_FALSE(vf::Math::equal(-1.0002f, -1.0001f));
	ASSERT_FALSE(vf::Math::equal(-1.0001f, -1.0002f));
}

TEST_F(MathTest, compareNumbersBetweenOneAndZero) {
	ASSERT_TRUE(vf::Math::equal(0.0000000010000001f, 0.0000000010000002f));
	ASSERT_TRUE(vf::Math::equal(0.0000000010000002f, 0.0000000010000001f));
	ASSERT_FALSE(vf::Math::equal(0.000000000001002f, 0.000000000001001f));
	ASSERT_FALSE(vf::Math::equal(0.000000000001001f, 0.000000000001002f));
}

TEST_F(MathTest, compareNumbersBetweenMinusOneAndZero) {
	ASSERT_TRUE(vf::Math::equal(-0.0000000010000001f, -0.0000000010000002f));
	ASSERT_TRUE(vf::Math::equal(-0.0000000010000002f, -0.0000000010000001f));
	ASSERT_FALSE(vf::Math::equal(-0.000000000001002f, -0.000000000001001f));
	ASSERT_FALSE(vf::Math::equal(-0.000000000001001f, -0.000000000001002f));
}

TEST_F(MathTest, compareNumbersZero) {
	ASSERT_TRUE(vf::Math::equal(0.0f, 0.0f));
	ASSERT_TRUE(vf::Math::equal(0.0f, -0.0f));
	ASSERT_TRUE(vf::Math::equal(-0.0f, -0.0f));
	ASSERT_FALSE(vf::Math::equal(0.000000001f, 0.0f));
	ASSERT_FALSE(vf::Math::equal(0.0f, 0.000000001f));
	ASSERT_FALSE(vf::Math::equal(-0.000000001f, 0.0f));
	ASSERT_FALSE(vf::Math::equal(0.0f, -0.000000001f));

	ASSERT_FALSE(vf::Math::equal(0.0f, 1e-40f, 0.01f));
	ASSERT_FALSE(vf::Math::equal(1e-40f, 0.0f, 0.01f));
}

TEST_F(MathTest, compareNumbersExtremeMaxValues) {
	ASSERT_TRUE(vf::Math::equal(FLT_MAX, FLT_MAX));
	ASSERT_FALSE(vf::Math::equal(FLT_MAX, -FLT_MAX));
	ASSERT_FALSE(vf::Math::equal(-FLT_MAX, FLT_MAX));
}

TEST_F(MathTest, compareNumbersExtremeMinValues) {
	ASSERT_TRUE(vf::Math::equal(FLT_MIN, FLT_MIN));

	ASSERT_FALSE(vf::Math::equal(FLT_MIN, 0.0f));
	ASSERT_FALSE(vf::Math::equal(0.0f, FLT_MIN));
	ASSERT_FALSE(vf::Math::equal(FLT_MIN, 0.0000000000001f));
	ASSERT_FALSE(vf::Math::equal(0.0000000000001f, FLT_MIN));

	ASSERT_FALSE(vf::Math::equal(FLT_MIN, -FLT_MIN));
	ASSERT_FALSE(vf::Math::equal(-FLT_MIN, FLT_MIN));
}

TEST_F(MathTest, compareNumbersNotANumber) {
	ASSERT_FALSE(vf::Math::equal(NAN, NAN));
	ASSERT_FALSE(vf::Math::equal(NAN, 0.0f));
	ASSERT_FALSE(vf::Math::equal(0.0f, NAN));
}


TEST_F(MathTest, lessEqual) {
	ASSERT_TRUE(vf::Math::lessEqual(10000001.0f, 10000000.0f));
	ASSERT_FALSE(vf::Math::lessEqual(10000001.0f, 1000000.0f));
	ASSERT_TRUE(vf::Math::lessEqual(100000.0f, 100000000.0f));
}

TEST_F(MathTest, greaterEqual) {
	ASSERT_TRUE(vf::Math::greaterEqual(10000001.0f, 10000000.0f));
	ASSERT_FALSE(vf::Math::greaterEqual(10000.0f, 100000001.0f));
	ASSERT_TRUE(vf::Math::greaterEqual(100000000.0f, 100000.0f));
}
