#ifndef _MPIIOSimulationObserver_H_
#define _MPIIOSimulationObserver_H_

#include "SimulationObserver.h"
#include <PointerDefinitions.h>
#include <mpi.h>
#include <string>

class Grid3D;
class UbScheduler;
namespace vf::mpi {class Communicator;}

//! \class MPIWriteBlocksBESimulationObserver
//! \brief Writes the grid each timestep into the files and reads the grip from the files before regenerating
class MPIIOSimulationObserver : public SimulationObserver
{
public:
    MPIIOSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path, std::shared_ptr<vf::mpi::Communicator> comm);
    ~MPIIOSimulationObserver() override;

    //! Each timestep writes the grid into the files
    void process(real step) override = 0;

    //! Writes the blocks of the grid into the file cpBlocks.bin
    void writeBlocks(int step);

    //! Reads the blocks of the grid from the file cpBlocks.bin
    void readBlocks(int step);

    //! The function truncates the data files
    void clearAllFiles(int step);

    //! The function write a time step of last check point
    void writeCpTimeStep(int step);
    //! The function read a time step of last check point
    int readCpTimeStep();

protected:
    std::string path;
    std::shared_ptr<vf::mpi::Communicator> comm;
    MPI_Datatype gridParamType, block3dType, dataSetParamType, boundCondType, arrayPresenceType;
};
#endif // ! _MPIIOSimulationObserver_H_
#define _MPIIOSimulationObserver_H_
