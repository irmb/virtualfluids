#include "FilePartCalculator.h"
#include "Parameter/Parameter.h"

uint FilePartCalculator::calculateNumberOfNodesInPart(const Parameter* parameter, uint level, uint indexOfFilePart)
{
    uint indexOfLastFilePart = FilePartCalculator::calculateNumberOfParts(parameter, level) - 1;

    if (indexOfFilePart > indexOfLastFilePart)
        throw std::runtime_error("The number of nodes for a non-existing part can not be calculated");
    if (indexOfFilePart == indexOfLastFilePart)
        return (uint)parameter->getParHostAsReference(level).numberOfNodes -
               (indexOfFilePart * FilePartCalculator::limitOfNodesForVTK);
    return FilePartCalculator::limitOfNodesForVTK;
}

uint FilePartCalculator::calculateNumberOfParts(const Parameter* parameter, uint level)
{
    return (uint)parameter->getParHostAsReference(level).numberOfNodes / FilePartCalculator::limitOfNodesForVTK + 1;
}

uint FilePartCalculator::getLimitOfNodesForVTK()
{
    return FilePartCalculator::limitOfNodesForVTK;
}

uint FilePartCalculator::calculateStartingPostionOfPart(uint indexOfPart)
{
    return indexOfPart * FilePartCalculator::limitOfNodesForVTK;
}
