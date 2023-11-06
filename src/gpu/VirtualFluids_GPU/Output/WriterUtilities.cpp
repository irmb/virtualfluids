#include "WriterUtilities.h"
#include "Parameter/Parameter.h"
#include <basics/StringUtilities/StringUtil.h>

std::string WriterUtilities::makePartFileNameEnding(uint level, int ID, int part, int timestep)
{
    return "_lev_" + StringUtil::toString<int>(level) + "_ID_" + StringUtil::toString<int>(ID) + "_Part_" +
           StringUtil::toString<int>(part) + "_t_" + StringUtil::toString<int>(timestep) + ".vtk";
}

uint WriterUtilities::calculateNumberOfParts(const Parameter* parameter, uint level)
{
    return (uint)parameter->getParHConst(level)->numberOfNodes / parameter->getLimitOfNodesForVTK() + 1;
}

bool WriterUtilities::isPeriodicCell(const Parameter* para, int level, unsigned int baseNodeOfCell,
                                     unsigned int otherNodeInCell)
{
    real distance = sqrt(
        pow(para->getParHConst(level)->coordinateX[otherNodeInCell] - para->getParHConst(level)->coordinateX[baseNodeOfCell], 2.) +
        pow(para->getParHConst(level)->coordinateY[otherNodeInCell] - para->getParHConst(level)->coordinateY[baseNodeOfCell], 2.) +
        pow(para->getParHConst(level)->coordinateZ[otherNodeInCell] - para->getParHConst(level)->coordinateZ[baseNodeOfCell], 2.));
    return distance > 1.01 * sqrt(3 * pow(para->getParHConst(level)->gridSpacing, 2.));
}

uint WriterUtilities::calculateNumberOfNodesInPart(const Parameter* para, uint level, uint part)
{
    if (part >= WriterUtilities::calculateNumberOfParts(para, level))
        throw std::runtime_error("The number of nodes for a non-existing part can not be calculated");
    if (((part + 1) * para->getLimitOfNodesForVTK()) > (uint)para->getParHConst(level)->numberOfNodes)
        return (uint)para->getParHConst(level)->numberOfNodes - (part * para->getLimitOfNodesForVTK());
    return para->getLimitOfNodesForVTK();
}

void WriterUtilities::getIndicesOfAllNodesInOct(std::array<uint, 8>& nodeIndices, uint baseNodeOfOct,
                                                const LBMSimulationParameter* parH)
{
    nodeIndices[0] = baseNodeOfOct;
    nodeIndices[1] = parH->neighborX[nodeIndices[0]];
    nodeIndices[2] = parH->neighborY[nodeIndices[1]];
    nodeIndices[3] = parH->neighborY[nodeIndices[0]];
    nodeIndices[4] = parH->neighborZ[nodeIndices[0]];
    nodeIndices[5] = parH->neighborZ[nodeIndices[1]];
    nodeIndices[6] = parH->neighborZ[nodeIndices[2]];
    nodeIndices[7] = parH->neighborZ[nodeIndices[3]];
}

void WriterUtilities::calculateRelativeNodeIndexInPart(std::array<uint, 8>& relativePositionInPart,
                                                      const std::array<uint, 8>& indicesOfOct, uint startPositionOfPart)
{
    for (size_t i = 0; i < relativePositionInPart.size(); i++) {
        relativePositionInPart[i] = indicesOfOct[i] - startPositionOfPart;
    }
}

bool WriterUtilities::areAllNodesInOctValidForWriting(const std::array<uint, 8>& indicesOfOct,
                                                      const LBMSimulationParameter* parH, uint endPositionOfPart)
{
    bool neighborsAreFluid = std::all_of(indicesOfOct.begin(), indicesOfOct.end(),
                                         [&](uint index) { return parH->typeOfGridNode[index] == GEO_FLUID; });

    bool neighborIsOutOfPart =
        (std::any_of(indicesOfOct.begin(), indicesOfOct.end(), [&](uint index) { return index > endPositionOfPart; }));

    return neighborsAreFluid && !neighborIsOutOfPart;
}