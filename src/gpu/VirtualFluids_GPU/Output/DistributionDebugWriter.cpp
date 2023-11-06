#include <basics/writer/WbWriterVtkXmlBinary.h>
#include <lbm/constants/D3Q27.h>
#include "DistributionDebugWriter.h"
#include "Parameter/Parameter.h"
#include "WriterUtilities.h"

using namespace vf::lbm::dir;

void DistributionDebugWriter::writeDistributions(const Parameter* para, uint timestep)
{
    for (int level = para->getCoarse(); level <= para->getFine(); level++) {
        DistributionDebugWriter::writeDistributionsForLevel(para, level, timestep);
    }
}

void DistributionDebugWriter::writeDistributionsForLevel(const Parameter* para, uint level, uint timestep)
{
    const uint numberOfParts = WriterUtilities::calculateNumberOfParts(para, level);

    std::vector<std::string> fileNames;
    for (uint i = 1; i <= numberOfParts; i++) {
        fileNames.push_back(para->getFName() + "_bin_distributions" +
                            WriterUtilities::makePartFileNameEnding(level, para->getMyProcessID(), i, timestep));
    }

    std::vector<std::string> nodeDataNames(NUMBER_Of_DIRECTIONS);

    for (uint dir = STARTDIR; dir <= ENDDIR; dir++) {
        const size_t minLenghtOfNumberString = 2; // the number is padded with zeros to this length
        const auto numberString = std::to_string(dir);
        nodeDataNames[dir] =
            "f_" + std::string(minLenghtOfNumberString - std::min(minLenghtOfNumberString, numberString.length()), '0') +
            numberString;
    }

    uint sizeOfNodes;
    uint startPosition;
    uint endPosition;
    std::array<uint, 8> indicesOfOct;
    std::array<uint, 8> relativePosInPart;
    uint relPosInPart;

    const LBMSimulationParameter* parH = para->getParHConst(level).get();
    Distributions27 distributions = parH->distributions;

    for (unsigned int part = 0; part < (uint)fileNames.size(); part++) {
        sizeOfNodes = WriterUtilities::calculateNumberOfNodesInPart(para, level, part);
        startPosition = part * para->getLimitOfNodesForVTK();
        endPosition = startPosition + sizeOfNodes;

        std::vector<UbTupleFloat3> nodes(sizeOfNodes);
        std::vector<UbTupleUInt8> cells;
        std::vector<std::vector<double>> nodeData(nodeDataNames.size());
        for (uint i = 0; i < (uint)nodeDataNames.size(); i++)
            nodeData[i].resize(sizeOfNodes);

        for (unsigned int pos = startPosition; pos < endPosition; pos++) {

            if (parH->typeOfGridNode[pos] != GEO_FLUID)
                continue;

            relPosInPart = pos - startPosition;

            // node
            double x1 = parH->coordinateX[pos];
            double x2 = parH->coordinateY[pos];
            double x3 = parH->coordinateZ[pos];
            nodes[relPosInPart] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));

            // node data
            for (uint dir = STARTDIR; dir <= ENDDIR; dir++) {
                nodeData[dir][relPosInPart] = distributions.f[0][dir*parH->numberOfNodes + pos];
            }

            WriterUtilities::getIndicesOfAllNodesInOct(indicesOfOct, pos, parH);
            if (WriterUtilities::isPeriodicCell(parH, indicesOfOct[0], indicesOfOct[6])) {
                continue;
            }

            if (WriterUtilities::areAllNodesInOctValidForWriting(indicesOfOct, parH, endPosition)) {
                WriterUtilities::calculateRelativeNodeIndexInPart(relativePosInPart, indicesOfOct, startPosition);
                cells.push_back(makeUbTupleFromArray(relativePosInPart));
            }
        }

        std::string fileName = WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fileNames[part], nodes, cells,
                                                                                          nodeDataNames, nodeData);
        VF_LOG_DEBUG("DistributionDebugWriter wrote to {} ", fileName);
    }
}
