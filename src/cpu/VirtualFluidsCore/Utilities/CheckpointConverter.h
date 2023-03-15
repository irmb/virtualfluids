#ifndef _UTILITACONVERTOR_H_
#define _UTILITACONVERTOR_H_

#include "MPIIODataStructures.h"
#include <PointerDefinitions.h>
#include <mpi.h>
#include <string>
#include <vector>

class Grid3D;
namespace vf::mpi {class Communicator;}

//! \class UtilConvertor
//! \brief Converts timestep data from MPIIORestartSimulationObserver format into MPIIOMigrationSimulationObserver format
class CheckpointConverter
{
public:
    CheckpointConverter(SPtr<Grid3D> grid, const std::string &path, std::shared_ptr<vf::mpi::Communicator> comm);
    virtual ~CheckpointConverter();

    void convert(int step, int procCount);
    void convertBlocks(int step, int procCount);
    void convertDataSet(int step, int procCount);
    void convertBC(int step, int procCount);
    void convert___Array(int step, int procCount, std::string filenameR, std::string filenameW);

protected:
    std::string path;
    std::shared_ptr<vf::mpi::Communicator> comm;
    SPtr<Grid3D> grid;

private:
    MPI_Datatype gridParamType, block3dType;
    MPI_Datatype dataSetParamType, dataSetTypeRead, dataSetTypeWrite;
    MPI_Datatype boundCondType, boundCondType1000;

    MPIIODataStructures::boundCondParam boundCondParamStr;
};

#endif
