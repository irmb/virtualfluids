#ifndef FILE_PART_CALCULATOR
#define FILE_PART_CALCULATOR

#include <basics/DataTypes.h>

class Parameter;

//! \brief calculations needed for writing the simulation's output into multiple file parts
//! \details To prevent output files from becoming to large, they are split into multiple part files. For this process some
//! calculations are needed.
class FilePartCalculator
{
private:
    //! \brief limits how many grid nodes are written into a single vtk file
    static const uint limitOfNodesForVTK = 30000000; // max 30 Million nodes per VTK file

public:
    //! \brief calculate how many output vtk-files are created for one timestep of the given grid level
    static uint calculateNumberOfParts(const Parameter* parameter, uint level);
    //! \brief calculate how many grid nodes are written to the file with the given index
    static uint calculateNumberOfNodesInPart(const Parameter* parameter, uint level, uint indexOfFilePart);
    //! \returns index of the first node in this file part
    static uint calculateStartingPostionOfPart(uint indexOfPart);
    static uint getLimitOfNodesForVTK();
};

#endif
