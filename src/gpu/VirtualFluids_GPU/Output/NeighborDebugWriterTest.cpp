#include <gmock/gmock.h>
#include "NeighborDebugWriter.hpp"
#include "gpu/VirtualFluids_GPU/Utilities/testUtilitiesGPU.h"

class WbWriterSpy : public WbWriter
{
public:
    std::string writeLines(const std::string & /*filename*/, std::vector<UbTupleFloat3> &nodes,
                           std::vector<UbTupleInt2> &lines) override
    {
        this->nodes = nodes;
        this->lines = lines;
        return "";
    }
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleInt2> lines;

    std::string getFileExtension() override { return ""; }
};

class NeighborDebugWriterTest : public testing::Test
{
protected:
    void SetUp() override
    {
        typeOfGridNode = std::vector<uint>(numberOfNodes, GEO_FLUID);
        neighbors = std::vector<uint>(numberOfNodes, 2);
        coordinates = std::vector<real>(numberOfNodes, 1.0);
        coordinates[2] = 3.0;

        para = testingVF::createParameterForLevel(0);
        para->getParH(level)->numberOfNodes = numberOfNodes;
        para->getParH(level)->coordinateX = coordinates.data();
        para->getParH(level)->coordinateY = coordinates.data();
        para->getParH(level)->coordinateZ = coordinates.data();
        para->getParH(level)->neighborX = neighbors.data();
        para->getParH(level)->typeOfGridNode = typeOfGridNode.data();
    }

    const int level = 0;
    const unsigned long long numberOfNodes = 3;
    const uint direction = vf::lbm::dir::DIR_P00; // x
    std::shared_ptr<Parameter> para;
    WbWriterSpy writerSpy;
    std::vector<uint> typeOfGridNode;
    std::vector<uint> neighbors;
    std::vector<real> coordinates;
};

TEST_F(NeighborDebugWriterTest, writeNeighborLinkLinesOnlyFLuidNodes)
{
    UbTupleFloat3 oneCoord(1.0, 1.0, 1.0);
    UbTupleFloat3 threeCoord(3.0, 3.0, 3.0);
    std::vector<UbTupleFloat3> expectedNodes = { oneCoord, threeCoord, oneCoord, threeCoord, threeCoord, threeCoord };
    std::vector<UbTupleInt2> expectedLines = { UbTupleInt2(0, 1), UbTupleInt2(2, 3), UbTupleInt2(4, 5) };

    NeighborDebugWriter::writeNeighborLinkLines(para->getParH(level).get(), direction, "name", &writerSpy);

    EXPECT_THAT(writerSpy.nodes.size(), testing::Eq(numberOfNodes * 2));
    EXPECT_THAT(writerSpy.lines.size(), testing::Eq(numberOfNodes));
    EXPECT_THAT(writerSpy.nodes, testing::Eq(expectedNodes));
    EXPECT_THAT(writerSpy.lines, testing::Eq(expectedLines));
}

TEST_F(NeighborDebugWriterTest, writeNeighborLinkLinesOneSolidNode)
{
    typeOfGridNode[2] = GEO_SOLID;
    
    UbTupleFloat3 oneCoord(1.0, 1.0, 1.0);
    UbTupleFloat3 threeCoord(3.0, 3.0, 3.0);
    std::vector<UbTupleFloat3> expectedNodes = { oneCoord, threeCoord, oneCoord, threeCoord};
    std::vector<UbTupleInt2> expectedLines = { UbTupleInt2(0, 1), UbTupleInt2(2, 3)};

    NeighborDebugWriter::writeNeighborLinkLines(para->getParH(level).get(), direction, "name", &writerSpy);

    EXPECT_THAT(writerSpy.nodes.size(), testing::Eq((numberOfNodes-1) * 2));
    EXPECT_THAT(writerSpy.lines.size(), testing::Eq(numberOfNodes-1));
    EXPECT_THAT(writerSpy.nodes, testing::Eq(expectedNodes));
    EXPECT_THAT(writerSpy.lines, testing::Eq(expectedLines));
}
