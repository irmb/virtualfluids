#include "FilePartCalculator.h"
#include "Parameter/Parameter.h"

uint FilePartCalculator::calculateNumberOfNodesInPart(const Parameter* para, uint level, uint indexOfFilePart)
{
    uint indexOfLastFilePart = FilePartCalculator::calculateNumberOfParts(para, level) - 1;

    if (indexOfFilePart > indexOfLastFilePart)
        throw std::runtime_error("The number of nodes for a non-existing part can not be calculated");
    if (indexOfFilePart == indexOfLastFilePart)
        return (uint)para->getParHostAsReference(level).numberOfNodes - (indexOfFilePart * para->getLimitOfNodesForVTK());
    return para->getLimitOfNodesForVTK();
}

uint FilePartCalculator::calculateNumberOfParts(const Parameter* parameter, uint level)
{
    return (uint)parameter->getParHostAsReference(level).numberOfNodes / parameter->getLimitOfNodesForVTK() + 1;
}
