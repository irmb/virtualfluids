#include "FilePartCalculator.h"
#include "gpu/core/Utilities/testUtilitiesGPU.h"
#include <gmock/gmock.h>

TEST(FilePartCalculatorTest, calculateNumberOfParts)
{
    const uint level = 0;
    auto parameter = testingVF::createParameterForLevel(level);
    parameter->setLimitOfNodesForVTK(10);

    parameter->getParH(level)->numberOfNodes = 23;
    EXPECT_THAT(FilePartCalculator::calculateNumberOfParts(parameter.get(), level), testing::Eq(3));

    parameter->getParH(level)->numberOfNodes = 13;
    EXPECT_THAT(FilePartCalculator::calculateNumberOfParts(parameter.get(), level), testing::Eq(2));

    parameter->getParH(level)->numberOfNodes = 3;
    EXPECT_THAT(FilePartCalculator::calculateNumberOfParts(parameter.get(), level), testing::Eq(1));
}

TEST(FilePartCalculatorTest, calculateNumberOfNodesInPart)
{
    const uint level = 0;
    auto parameter = testingVF::createParameterForLevel(level);
    parameter->getParH(level)->numberOfNodes = 13;
    parameter->setLimitOfNodesForVTK(10);

    uint part = 0;
    EXPECT_THAT(FilePartCalculator::calculateNumberOfNodesInPart(parameter.get(), level, part), testing::Eq(10));

    part = 1;
    EXPECT_THAT(FilePartCalculator::calculateNumberOfNodesInPart(parameter.get(), level, part), testing::Eq(3));

    part = 2;
    EXPECT_THROW(FilePartCalculator::calculateNumberOfNodesInPart(parameter.get(), level, part), std::runtime_error);
}
