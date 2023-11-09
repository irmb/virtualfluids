#include "FilePartCalculator.h"
#include "gpu/core/Utilities/testUtilitiesGPU.h"
#include <gmock/gmock.h>

TEST(FilePartCalculatorTest, calculateNumberOfParts)
{
    const uint level = 0;
    auto parameter = testingVF::createParameterForLevel(level);

    parameter->getParH(level)->numberOfNodes = FilePartCalculator::getLimitOfNodesForVTK() * 2;
    EXPECT_THAT(FilePartCalculator::calculateNumberOfParts(parameter.get(), level), testing::Eq(3));

    parameter->getParH(level)->numberOfNodes = FilePartCalculator::getLimitOfNodesForVTK();
    EXPECT_THAT(FilePartCalculator::calculateNumberOfParts(parameter.get(), level), testing::Eq(2));

    parameter->getParH(level)->numberOfNodes = FilePartCalculator::getLimitOfNodesForVTK() - 1;
    EXPECT_THAT(FilePartCalculator::calculateNumberOfParts(parameter.get(), level), testing::Eq(1));

    parameter->getParH(level)->numberOfNodes = 1;
    EXPECT_THAT(FilePartCalculator::calculateNumberOfParts(parameter.get(), level), testing::Eq(1));
}

TEST(FilePartCalculatorTest, calculateNumberOfNodesInPart)
{
    const uint level = 0;
    auto parameter = testingVF::createParameterForLevel(level);
    parameter->getParH(level)->numberOfNodes = FilePartCalculator::getLimitOfNodesForVTK() + 13;

    uint part = 0;
    EXPECT_THAT(FilePartCalculator::calculateNumberOfNodesInPart(parameter.get(), level, part),
                FilePartCalculator::getLimitOfNodesForVTK());

    part = 1;
    EXPECT_THAT(FilePartCalculator::calculateNumberOfNodesInPart(parameter.get(), level, part), 13);

    part = 2;
    EXPECT_THROW(FilePartCalculator::calculateNumberOfNodesInPart(parameter.get(), level, part), std::runtime_error);
}

TEST(FilePartCalculatorTest, getStartingPostionOfPart)
{
    const uint level = 0;
    auto parameter = testingVF::createParameterForLevel(level);
    parameter->getParH(level)->numberOfNodes = FilePartCalculator::getLimitOfNodesForVTK() + 13;

    uint part = 0;
    EXPECT_THAT(FilePartCalculator::calculateStartingPostionOfPart(part),
                0);

    part = 1;
    EXPECT_THAT(FilePartCalculator::calculateStartingPostionOfPart(part), FilePartCalculator::getLimitOfNodesForVTK());

    part = 2;
    EXPECT_THROW(FilePartCalculator::calculateNumberOfNodesInPart(parameter.get(), level, part), std::runtime_error);
}
