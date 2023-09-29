#ifndef WRITER_UTILITIES
#define WRITER_UTILITIES

#include <array>

#include <basics/DataTypes.h>

class Parameter;
class LBMSimulationParameter;

class WriterUtilities
{
public:
    static uint calculateNumberOfParts(const Parameter* para, uint level);
    static uint calculateNumberOfNodesInPart(const Parameter* para, uint level, uint part);
    static std::string makePartFileNameEnding(uint level, int processID, int part, int timestep);
    static void getIndicesOfAllNodesInOct(std::array<uint, 8>& indices, uint baseNodeOfOct,
                                          const LBMSimulationParameter* parH);
    static void calculateRelativePositionInPart(std::array<uint, 8>& relativePositionInPart,
                                                const std::array<uint, 8>& indicesOfOct, uint startPosition);
    static bool areAllNodesInOctValidForWriting(const std::array<uint, 8>& indicesOfOct, const LBMSimulationParameter* parH,
                                                uint endPositionOfPart);
    static bool isPeriodicCell(const Parameter* para, int level, unsigned int number1, unsigned int number7);
};

#endif
