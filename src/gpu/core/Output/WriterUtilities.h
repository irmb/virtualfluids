#ifndef WRITER_UTILITIES
#define WRITER_UTILITIES

#include <array>

#include <basics/DataTypes.h>

class Parameter;
struct LBMSimulationParameter;

class WriterUtilities
{
public:
    //! \brief check whether a grid cell is part of a periodic boundary condition
    //! \param baseNodeOfCell is the index of one node of the grid cell
    //! \param otherNodeInCell is the index of the node which is not on the same face of the grid cell as the base node (i.e.
    //! it is on the other end of the space diagonal)
    static bool isPeriodicCell(const LBMSimulationParameter& parH, unsigned int baseNodeOfCell,
                               unsigned int otherNodeInCell);

    //! \brief use the neighbor relations to find the indices of all nodes in an oct cell
    static void getIndicesOfAllNodesInOct(std::array<uint, 8>& nodeIndices, uint baseNodeOfOct,
                                          const LBMSimulationParameter& parH);

    //! \brief calculate the node index relative to the start position of the part
    static void calculateRelativeNodeIndexInPart(std::array<uint, 8>& relativePositionInPart,
                                                 const std::array<uint, 8>& indicesOfOct, uint startPositionOfPart);

    //! \brief check if all nodes in an oct are valid to be written into an output file
    //! \details to be valid the nodes need to be: 1. have the type GEO_FLUID, 2. not be outside the current file part
    //! \param endPositionOfPart specifies the index of the last node in the current file part
    static bool areAllNodesInOctValidForWriting(const std::array<uint, 8>& indicesOfOct, const LBMSimulationParameter& parH,
                                                uint endPositionOfPart);

    //! \brief create the ending of the file name for a file part
    static std::string makePartFileNameEnding(uint level, int processID, int part, int timestep);
};

#endif
