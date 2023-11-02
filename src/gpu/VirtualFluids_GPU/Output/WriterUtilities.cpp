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

bool WriterUtilities::isPeriodicCell(const Parameter* para, int level, unsigned int number1, unsigned int number7)
{
    real distance =
        sqrt(pow(para->getParHConst(level)->coordinateX[number7] - para->getParHConst(level)->coordinateX[number1], 2.) +
             pow(para->getParHConst(level)->coordinateY[number7] - para->getParHConst(level)->coordinateY[number1], 2.) +
             pow(para->getParHConst(level)->coordinateZ[number7] - para->getParHConst(level)->coordinateZ[number1], 2.));
    return distance > 1.01 * sqrt(3 * para->getParHConst(level)->gridSpacing);
}

uint WriterUtilities::calculateNumberOfNodesInPart(const Parameter* para, uint level, uint part)
{
    if (((part + 1) * para->getLimitOfNodesForVTK()) > (uint)para->getParHConst(level)->numberOfNodes)
        return (uint)para->getParHConst(level)->numberOfNodes - (part * para->getLimitOfNodesForVTK());
    return para->getLimitOfNodesForVTK();
}

void WriterUtilities::getIndicesOfAllNodesInOct(std::array<uint, 8>& indices, uint baseNodeOfOct,
                                                const LBMSimulationParameter* parH)
{
    indices[0] = baseNodeOfOct;
    indices[1] = parH->neighborX[indices[0]];
    indices[2] = parH->neighborY[indices[1]];
    indices[3] = parH->neighborY[indices[0]];
    indices[4] = parH->neighborZ[indices[0]];
    indices[5] = parH->neighborZ[indices[1]];
    indices[6] = parH->neighborZ[indices[2]];
    indices[7] = parH->neighborZ[indices[3]];
}

void WriterUtilities::calculateRelativePositionInPart(std::array<uint, 8>& relativePositionInPart,
                                                      const std::array<uint, 8>& indicesOfOct, uint startPosition)
{
    for (size_t i = 0; i < relativePositionInPart.size(); i++) {
        relativePositionInPart[i] = indicesOfOct[i] - startPosition;
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