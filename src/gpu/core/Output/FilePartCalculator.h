#ifndef FILE_PART_CALCULATOR
#define FILE_PART_CALCULATOR

#include <basics/DataTypes.h>

//! \brief calculations needed for writing the simulation's output into multiple file parts
//! \details To prevent output files from becoming to large, they are split into multiple part files. For this process some
//! calculations are needed.
class FilePartCalculator
{
public:
    //! \brief calculate how many output vtk-files are created for one timestep of the given grid level
    static uint calculateNumberOfParts(uint numberOfNodes);
    //! \brief calculate how many grid nodes are written to the file with the given index
    static uint calculateNumberOfNodesInPart(uint numberOfNodes, uint indexOfFilePart);
    //! \returns index of the first node in this file part
    static uint calculateStartingPostionOfPart(uint indexOfPart);

    //! \brief limits how many grid nodes are written into a single vtk file
    static const uint limitOfNodesForVTK = 30e6; // max 30 million nodes per VTK file
};

#endif
