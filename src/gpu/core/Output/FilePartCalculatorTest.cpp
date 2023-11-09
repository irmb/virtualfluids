#include "FilePartCalculator.h"
#include "gpu/core/Utilities/testUtilitiesGPU.h"
#include <gmock/gmock.h>

TEST(FilePartCalculatorTest, calculateNumberOfParts)
{
    EXPECT_THAT(FilePartCalculator::calculateNumberOfParts(FilePartCalculator::limitOfNodesForVTK * 2), testing::Eq(3));

    EXPECT_THAT(FilePartCalculator::calculateNumberOfParts(FilePartCalculator::limitOfNodesForVTK), testing::Eq(2));

    EXPECT_THAT(FilePartCalculator::calculateNumberOfParts( FilePartCalculator::limitOfNodesForVTK - 1), testing::Eq(1));

    EXPECT_THAT(FilePartCalculator::calculateNumberOfParts(1), testing::Eq(1));
}

TEST(FilePartCalculatorTest, calculateNumberOfNodesInPart)
{
    uint part = 0;
    EXPECT_THAT(FilePartCalculator::calculateNumberOfNodesInPart(FilePartCalculator::limitOfNodesForVTK + 13, part),
                FilePartCalculator::limitOfNodesForVTK);

    part = 1;
    EXPECT_THAT(FilePartCalculator::calculateNumberOfNodesInPart(FilePartCalculator::limitOfNodesForVTK + 13, part), 13);

    part = 2;
    EXPECT_THROW(FilePartCalculator::calculateNumberOfNodesInPart(FilePartCalculator::limitOfNodesForVTK + 13, part), std::runtime_error);
}

TEST(FilePartCalculatorTest, getStartingPostionOfPart)
{
    uint part = 0;
    EXPECT_THAT(FilePartCalculator::calculateStartingPostionOfPart(part), 0);

    part = 1;
    EXPECT_THAT(FilePartCalculator::calculateStartingPostionOfPart(part), FilePartCalculator::limitOfNodesForVTK);
}
