#include "FilePartCalculator.h"
#include <stdexcept>

uint FilePartCalculator::calculateNumberOfNodesInPart(uint numberOfNodes, uint indexOfFilePart)
{
    uint indexOfLastFilePart = FilePartCalculator::calculateNumberOfParts(numberOfNodes) - 1;

    if (indexOfFilePart > indexOfLastFilePart)
        throw std::runtime_error("The number of nodes for a non-existing part can not be calculated");
    if (indexOfFilePart == indexOfLastFilePart)
        return numberOfNodes - (indexOfFilePart * FilePartCalculator::limitOfNodesForVTK);
    return FilePartCalculator::limitOfNodesForVTK;
}

uint FilePartCalculator::calculateNumberOfParts(uint numberOfNodes)
{
    return numberOfNodes / FilePartCalculator::limitOfNodesForVTK + 1;
}

uint FilePartCalculator::calculateStartingPostionOfPart(uint indexOfPart)
{
    return indexOfPart * FilePartCalculator::limitOfNodesForVTK;
}

const uint FilePartCalculator::limitOfNodesForVTK;
