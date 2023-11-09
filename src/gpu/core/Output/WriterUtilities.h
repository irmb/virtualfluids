#ifndef WRITER_UTILITIES
#define WRITER_UTILITIES

#include <array>

#include <basics/DataTypes.h>

class Parameter;
struct LBMSimulationParameter;

class WriterUtilities
{
public:
    static uint calculateNumberOfParts(const Parameter* para, uint level);
    static uint calculateNumberOfNodesInPart(const Parameter* para, uint level, uint part);
    static bool isPeriodicCell(const LBMSimulationParameter* parH, unsigned int baseNodeOfCell,
                               unsigned int otherNodeInCell);
    static void getIndicesOfAllNodesInOct(std::array<uint, 8>& nodeIndices, uint baseNodeOfOct,
                                          const LBMSimulationParameter* parH);
    static void calculateRelativeNodeIndexInPart(std::array<uint, 8>& relativePositionInPart,
                                                 const std::array<uint, 8>& indicesOfOct, uint startPositionOfPart);
    static bool areAllNodesInOctValidForWriting(const std::array<uint, 8>& indicesOfOct, const LBMSimulationParameter* parH,
                                                uint endPositionOfPart);

    static std::string makePartFileNameEnding(uint level, int processID, int part, int timestep);
};

#endif
