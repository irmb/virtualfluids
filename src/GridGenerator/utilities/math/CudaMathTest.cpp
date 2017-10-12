#include "gmock/gmock.h"
#include "CudaMath.cuh"

using namespace testing;

class CudaMathTest : public Test {
public:
	void SetUp() {

	}
};

TEST_F(CudaMathTest, compareBigValues) {
	ASSERT_TRUE(CudaMath::equal(10000000.0f, 10000001.0f));
	ASSERT_TRUE(CudaMath::equal(10000001.0f, 10000000.0f));
	ASSERT_FALSE(CudaMath::equal(100000.0f, 100001.0f));
	ASSERT_FALSE(CudaMath::equal(100001.0f, 100000.0f));
}

TEST_F(CudaMathTest, compareSmallValues) {
	ASSERT_TRUE(CudaMath::equal(-10000000.0f, -10000001.0f));
	ASSERT_TRUE(CudaMath::equal(-10000001.0f, -10000000.0f));
	ASSERT_FALSE(CudaMath::equal(-100000.0f, -100001.0f));
	ASSERT_FALSE(CudaMath::equal(-100001.0f, -100000.0f));
}

TEST_F(CudaMathTest, compareNumbersAroundOne) {
	ASSERT_TRUE(CudaMath::equal(1.0000001f, 1.0000002f));
	ASSERT_TRUE(CudaMath::equal(1.0000002f, 1.0000001f));
	ASSERT_FALSE(CudaMath::equal(1.0002f, 1.0001f));
	ASSERT_FALSE(CudaMath::equal(1.0001f, 1.0002f));
}

TEST_F(CudaMathTest, compareNumbersAroundMinusOne) {
	ASSERT_TRUE(CudaMath::equal(-1.0000001f, -1.0000002f));
	ASSERT_TRUE(CudaMath::equal(-1.0000002f, -1.0000001f));
	ASSERT_FALSE(CudaMath::equal(-1.0002f, -1.0001f));
	ASSERT_FALSE(CudaMath::equal(-1.0001f, -1.0002f));
}

TEST_F(CudaMathTest, compareNumbersBetweenOneAndZero) {
	ASSERT_TRUE(CudaMath::equal(0.0000000010000001f, 0.0000000010000002f));
	ASSERT_TRUE(CudaMath::equal(0.0000000010000002f, 0.0000000010000001f));
	ASSERT_FALSE(CudaMath::equal(0.000000000001002f, 0.000000000001001f));
	ASSERT_FALSE(CudaMath::equal(0.000000000001001f, 0.000000000001002f));
}

TEST_F(CudaMathTest, compareNumbersBetweenMinusOneAndZero) {
	ASSERT_TRUE(CudaMath::equal(-0.0000000010000001f, -0.0000000010000002f));
	ASSERT_TRUE(CudaMath::equal(-0.0000000010000002f, -0.0000000010000001f));
	ASSERT_FALSE(CudaMath::equal(-0.000000000001002f, -0.000000000001001f));
	ASSERT_FALSE(CudaMath::equal(-0.000000000001001f, -0.000000000001002f));
}

TEST_F(CudaMathTest, compareNumbersZero) {
	ASSERT_TRUE(CudaMath::equal(0.0f, 0.0f));
	ASSERT_TRUE(CudaMath::equal(0.0f, -0.0f));
	ASSERT_TRUE(CudaMath::equal(-0.0f, -0.0f));
	ASSERT_FALSE(CudaMath::equal(0.000000001f, 0.0f));
	ASSERT_FALSE(CudaMath::equal(0.0f, 0.000000001f));
	ASSERT_FALSE(CudaMath::equal(-0.000000001f, 0.0f));
	ASSERT_FALSE(CudaMath::equal(0.0f, -0.000000001f));

	ASSERT_FALSE(CudaMath::equal(0.0f, 1e-40f, 0.01f));
	ASSERT_FALSE(CudaMath::equal(1e-40f, 0.0f, 0.01f));
}

TEST_F(CudaMathTest, compareNumbersExtremeMaxValues) {
	ASSERT_TRUE(CudaMath::equal(FLT_MAX, FLT_MAX));
	ASSERT_FALSE(CudaMath::equal(FLT_MAX, -FLT_MAX));
	ASSERT_FALSE(CudaMath::equal(-FLT_MAX, FLT_MAX));
}

TEST_F(CudaMathTest, compareNumbersExtremeMinValues) {
	ASSERT_TRUE(CudaMath::equal(FLT_MIN, FLT_MIN));

	ASSERT_FALSE(CudaMath::equal(FLT_MIN, 0.0f));
	ASSERT_FALSE(CudaMath::equal(0.0f, FLT_MIN));
	ASSERT_FALSE(CudaMath::equal(FLT_MIN, 0.0000000000001f));
	ASSERT_FALSE(CudaMath::equal(0.0000000000001f, FLT_MIN));

	ASSERT_FALSE(CudaMath::equal(FLT_MIN, -FLT_MIN));
	ASSERT_FALSE(CudaMath::equal(-FLT_MIN, FLT_MIN));
}

TEST_F(CudaMathTest, compareNumbersNotANumber) {
	ASSERT_FALSE(CudaMath::equal(NAN, NAN));
	ASSERT_FALSE(CudaMath::equal(NAN, 0.0f));
	ASSERT_FALSE(CudaMath::equal(0.0f, NAN));
}


TEST_F(CudaMathTest, lessEqual) {
	ASSERT_TRUE(CudaMath::lessEqual(10000001.0f, 10000000.0f));
	ASSERT_FALSE(CudaMath::lessEqual(10000001.0f, 1000000.0f));
	ASSERT_TRUE(CudaMath::lessEqual(100000.0f, 100000000.0f));
}

TEST_F(CudaMathTest, greaterEqual) {
	ASSERT_TRUE(CudaMath::greaterEqual(10000001.0f, 10000000.0f));
	ASSERT_FALSE(CudaMath::greaterEqual(10000.0f, 100000001.0f));
	ASSERT_TRUE(CudaMath::greaterEqual(100000000.0f, 100000.0f));
}
