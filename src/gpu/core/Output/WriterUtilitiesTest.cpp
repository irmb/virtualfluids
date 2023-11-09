
#include "WriterUtilities.h"
#include "gpu/VirtualFluids_GPU/Utilities/testUtilitiesGPU.h"
#include <gmock/gmock.h>
#include <numeric>

TEST(WriterUtilitiesTest, calculateNumberOfParts)
{
    const uint level = 0;
    auto parameter = testingVF::createParameterForLevel(level);
    parameter->setLimitOfNodesForVTK(10);

    parameter->getParH(level)->numberOfNodes = 23;
    EXPECT_THAT(WriterUtilities::calculateNumberOfParts(parameter.get(), level), testing::Eq(3));

    parameter->getParH(level)->numberOfNodes = 13;
    EXPECT_THAT(WriterUtilities::calculateNumberOfParts(parameter.get(), level), testing::Eq(2));

    parameter->getParH(level)->numberOfNodes = 3;
    EXPECT_THAT(WriterUtilities::calculateNumberOfParts(parameter.get(), level), testing::Eq(1));
}

TEST(WriterUtilitiesTest, calculateNumberOfNodesInPart)
{
    const uint level = 0;
    auto parameter = testingVF::createParameterForLevel(level);
    parameter->getParH(level)->numberOfNodes = 13;
    parameter->setLimitOfNodesForVTK(10);

    uint part = 0;
    EXPECT_THAT(WriterUtilities::calculateNumberOfNodesInPart(parameter.get(), level, part), testing::Eq(10));

    part = 1;
    EXPECT_THAT(WriterUtilities::calculateNumberOfNodesInPart(parameter.get(), level, part), testing::Eq(3));

    part = 2;
    EXPECT_THROW(WriterUtilities::calculateNumberOfNodesInPart(parameter.get(), level, part), std::runtime_error);
}

class WriterUtilitiesPeriodicCellTest : public testing::Test
{
protected:
    std::unique_ptr<LBMSimulationParameter> parH = std::make_unique<LBMSimulationParameter>();
    const uint level = 0;
    const uint baseNodeIndex = 0;
    const uint otherNodeIndex = 1;
    std::array<real, 2> coordinates = {0.0, 1.0};

    void SetUp() override
    {
        // create a domain with only three layers of nodes
        // nodes are at the coordinates 0.0, 1.0 and 2.0
        parH->gridSpacing = 1.0;

        parH->coordinateX = new real[2];
        parH->coordinateY = new real[2];
        parH->coordinateZ = new real[2];

        parH->coordinateX[baseNodeIndex] = coordinates[baseNodeIndex];
        parH->coordinateY[baseNodeIndex] = coordinates[baseNodeIndex];
        parH->coordinateZ[baseNodeIndex] = coordinates[baseNodeIndex];
        parH->coordinateX[otherNodeIndex] = coordinates[otherNodeIndex];
        parH->coordinateY[otherNodeIndex] = coordinates[otherNodeIndex];
        parH->coordinateZ[otherNodeIndex] = coordinates[otherNodeIndex];
    }

    void TearDown() override
    {
        delete[] parH->coordinateX;
        delete[] parH->coordinateY;
        delete[] parH->coordinateZ;
    }
};

TEST_F(WriterUtilitiesPeriodicCellTest, cellIsNotPeriodic)
{
    EXPECT_FALSE(WriterUtilities::isPeriodicCell(parH.get(), baseNodeIndex, otherNodeIndex));
    EXPECT_FALSE(WriterUtilities::isPeriodicCell(parH.get(), otherNodeIndex, baseNodeIndex));
}

TEST_F(WriterUtilitiesPeriodicCellTest, cellIsPeriodicInX)
{
    parH->coordinateX[1] = 2.0;
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH.get(), baseNodeIndex, otherNodeIndex));
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH.get(), otherNodeIndex, baseNodeIndex));
}

TEST_F(WriterUtilitiesPeriodicCellTest, cellIsPeriodicInY)
{
    parH->coordinateY[1] = 2.0;
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH.get(), baseNodeIndex, otherNodeIndex));
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH.get(), otherNodeIndex, baseNodeIndex));
}

TEST_F(WriterUtilitiesPeriodicCellTest, cellIsPeriodicInZ)
{
    parH->coordinateZ[1] = 2.0;
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH.get(), baseNodeIndex, otherNodeIndex));
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH.get(), otherNodeIndex, baseNodeIndex));
}
TEST_F(WriterUtilitiesPeriodicCellTest, cellIsPeriodicInXY)
{
    parH->coordinateX[1] = 2.0;
    parH->coordinateY[1] = 2.0;
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH.get(), baseNodeIndex, otherNodeIndex));
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH.get(), otherNodeIndex, baseNodeIndex));
}

TEST_F(WriterUtilitiesPeriodicCellTest, cellIsPeriodicInXZ)
{
    parH->coordinateX[1] = 2.0;
    parH->coordinateZ[1] = 2.0;
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH.get(), baseNodeIndex, otherNodeIndex));
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH.get(), otherNodeIndex, baseNodeIndex));
}

TEST_F(WriterUtilitiesPeriodicCellTest, cellIsPeriodicInYZ)
{
    parH->coordinateY[1] = 2.0;
    parH->coordinateZ[1] = 2.0;
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH.get(), baseNodeIndex, otherNodeIndex));
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH.get(), otherNodeIndex, baseNodeIndex));
}

TEST_F(WriterUtilitiesPeriodicCellTest, cellIsPeriodicInXYZ)
{
    parH->coordinateX[1] = 2.0;
    parH->coordinateY[1] = 2.0;
    parH->coordinateZ[1] = 2.0;
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH.get(), baseNodeIndex, otherNodeIndex));
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH.get(), otherNodeIndex, baseNodeIndex));
}

class WriterUtilitiesNeighborOctTest : public testing::Test
{
    static void setUpNeighborsNeighborsForOct(LBMSimulationParameter* parH, const std::array<uint, 8>& nodeIndices)
    {
        // node indices: MMM, PMM, PPM, MPM,
        //               MMP, PMP, PPP, MPP

        for (uint i = 0; i < (uint)nodeIndices.size(); i++) {
            const uint currentNodeIndex = nodeIndices[i];
            if (i < 4)
                parH->neighborZ[currentNodeIndex] = nodeIndices[i + 4];
            else
                parH->neighborZ[currentNodeIndex] = 99;

            if (i == 0 || i == 4)
                parH->neighborY[currentNodeIndex] = nodeIndices[i + 3];
            else if (i == 1 || i == 5)
                parH->neighborY[currentNodeIndex] = nodeIndices[i + 1];
            else
                parH->neighborY[currentNodeIndex] = 999;

            if (i == 0 || i == 4)
                parH->neighborX[currentNodeIndex] = nodeIndices[i + 1];
            else if (i == 3 || i == 7)
                parH->neighborX[currentNodeIndex] = nodeIndices[i - 1];
            else
                parH->neighborX[currentNodeIndex] = 9999;
        }
    }

public:
    std::unique_ptr<LBMSimulationParameter> parH = std::make_unique<LBMSimulationParameter>();
    std::array<uint, 8> nodeIndices;

    void SetUp() override
    {
        // set up some node indices from 0 to 7
        std::iota(nodeIndices.begin(), nodeIndices.end(), 0);
        std::reverse(nodeIndices.begin(), nodeIndices.end());

        parH->neighborX = new uint[8];
        parH->neighborY = new uint[8];
        parH->neighborZ = new uint[8];
        setUpNeighborsNeighborsForOct(parH.get(), nodeIndices);
    }

    void TearDown() override
    {
        delete[] parH->neighborX;
        delete[] parH->neighborY;
        delete[] parH->neighborZ;
    }
};

TEST_F(WriterUtilitiesNeighborOctTest, getIndicesOfAllNodesInOct)
{
    std::array<uint, 8> resultingNodeIndices;
    WriterUtilities::getIndicesOfAllNodesInOct(resultingNodeIndices, nodeIndices[0], parH.get());
    for (uint i = 0; i < 8; i++)
        EXPECT_THAT(resultingNodeIndices[i], testing::Eq(nodeIndices[i])) << "for index i = " << i << " in nodeIndices";
}

TEST(WriterUtilitiesTest, calculateRelativeNodeIndexInPart)
{
    std::array<uint, 8> indicesOfOct = { 10, 12, 14, 16, 20, 22, 24, 26 };
    uint startPositionOfPart = 10;
    std::array<uint, 8> expected = { 0, 2, 4, 6, 10, 12, 14, 16 };

    std::array<uint, 8> result;
    WriterUtilities::calculateRelativeNodeIndexInPart(result, indicesOfOct, startPositionOfPart);
    EXPECT_THAT(result, testing::Eq(expected));
}

class WriterUtilitiesTestNodeValidity : public testing::Test
{
protected:
    std::unique_ptr<LBMSimulationParameter> parH = std::make_unique<LBMSimulationParameter>();
    std::array<uint, 8> nodeIndices;
    std::array<uint, 8> typeOfGridNode;

    void SetUp() override
    {
        // set up node indices from 0 to 7
        std::iota(nodeIndices.begin(), nodeIndices.end(), 0);

        std::fill(typeOfGridNode.begin(), typeOfGridNode.end(), GEO_FLUID);
        parH->typeOfGridNode = typeOfGridNode.data();
    }
};

TEST_F(WriterUtilitiesTestNodeValidity, allNodesInOctValidForWriting)
{
    uint endPositionOfPart = 7;
    EXPECT_TRUE(WriterUtilities::areAllNodesInOctValidForWriting(nodeIndices, parH.get(), endPositionOfPart));
}

TEST_F(WriterUtilitiesTestNodeValidity, areAllNodesInOctValidForWriting_NodeOutOfPart)
{
    uint endPositionOfPart = 6;
    EXPECT_FALSE(WriterUtilities::areAllNodesInOctValidForWriting(nodeIndices, parH.get(), endPositionOfPart));
}

TEST_F(WriterUtilitiesTestNodeValidity, areAllNodesInOctValidForWriting_NonFluidNode)
{
    uint endPositionOfPart = 7;
    typeOfGridNode[0] = GEO_SOLID;
    EXPECT_FALSE(WriterUtilities::areAllNodesInOctValidForWriting(nodeIndices, parH.get(), endPositionOfPart));
}

TEST_F(WriterUtilitiesTestNodeValidity, areAllNodesInOctValidForWriting_NonFluidNodeAtEnd)
{
    uint endPositionOfPart = 7;
    typeOfGridNode[7] = GEO_SOLID;
    EXPECT_FALSE(WriterUtilities::areAllNodesInOctValidForWriting(nodeIndices, parH.get(), endPositionOfPart));
}
