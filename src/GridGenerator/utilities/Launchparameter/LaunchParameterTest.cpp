#include "LaunchParameter.cuh"

#include "gmock/gmock.h"


using namespace testing;


TEST(LaunchParameterTest, get1d1dLaunchParameterWithThreadsEqualToSize)
{
	int threads = 256;
	int size = 256;

	LaunchParameter p = LaunchParameter::make_1D1D_launchParameter(size, threads);

	EXPECT_THAT(p.threads.x, Eq(256));
	EXPECT_THAT(p.threads.y, Eq(1));
	EXPECT_THAT(p.threads.z, Eq(1));

	EXPECT_THAT(p.blocks.x, Eq(1));
	EXPECT_THAT(p.blocks.y, Eq(1));
	EXPECT_THAT(p.blocks.z, Eq(1));
}


TEST(LaunchParameterTest, get1d1dLaunchParameterWithThreadsSmallerThanSize)
{
	int threads = 256;
	int size = 257;

	LaunchParameter p = LaunchParameter::make_1D1D_launchParameter(size, threads);

	EXPECT_THAT(p.threads.x, Eq(256));
	EXPECT_THAT(p.threads.y, Eq(1));
	EXPECT_THAT(p.threads.z, Eq(1));

	EXPECT_THAT(p.blocks.x, Eq(2));
	EXPECT_THAT(p.blocks.y, Eq(1));
	EXPECT_THAT(p.blocks.z, Eq(1));
}


TEST(LaunchParameterTest, get2d1dLaunchParameterWithThreadsSmallerThanSize)
{
	int threads = 256;
	int size = 65535 * 257;

	LaunchParameter p = LaunchParameter::make_2D1D_launchParameter(size, threads);

	EXPECT_THAT(p.threads.x, Eq(256));
	EXPECT_THAT(p.threads.y, Eq(1));
	EXPECT_THAT(p.threads.z, Eq(1));

	EXPECT_THAT(p.blocks.x, Eq(256));
	EXPECT_THAT(p.blocks.y, Eq(256));
	EXPECT_THAT(p.blocks.z, Eq(1));
}