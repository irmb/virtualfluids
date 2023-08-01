#ifndef _MPIIORestartSimulationObserver_H_
#define _MPIIORestartSimulationObserver_H_

#include <mpi.h>
//#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "MPIIOSimulationObserver.h"
#include "MPIIODataStructures.h"

class Grid3D;
class UbScheduler;
namespace vf::parallel {class Communicator;}
class BCSet;
class LBMKernel;

//! \class MPIIORestartSimulationObserver
//! \brief Writes the grid each timestep into the files and reads the grip from the files before regenerating
class MPIIORestartSimulationObserver : public MPIIOSimulationObserver
{
public:
    enum Arrays {
        AverageDensity = 1,
        AverageVelocity = 2,
        AverageFluktuations = 3,
        AverageTriple = 4,
        ShearStressVal = 5,
        RelaxationFactor = 6,
        PhaseField1 = 7,
        PhaseField2 = 8,
        PressureField = 9
    };

    MPIIORestartSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path, std::shared_ptr<vf::parallel::Communicator> comm);
    ~MPIIORestartSimulationObserver() override;
    //! Each timestep writes the grid into the files
    void update(real step) override;
    //! Reads the grid from the files before grid reconstruction
    void restart(int step);
    //! Writes the blocks of the grid into the file cpBlocks.bin
    void writeBlocks(int step);
    //! Writes the datasets of the blocks into the file cpDataSet.bin
    void writeDataSet(int step);
    void write4DArray(int step, Arrays arrType, std::string fname);
    void write3DArray(int step, Arrays arrType, std::string fname);
    //void writeAverageDensityArray(int step);
    //void writeAverageVelocityArray(int step);
    //void writeAverageFluktuationsArray(int step);
    //void writeAverageTripleArray(int step);
    //void writeShearStressValArray(int step);
    //void writeRelaxationFactor(int step);
    //void writePhaseField(int step, int num);
    //void writePressureField(int step);
    //! Writes the boundary conditions of the blocks into the file cpBC.bin
    void writeBoundaryConds(int step);

    //! Reads the blocks of the grid from the file cpBlocks.bin
    void readBlocks(int step);
    //! Reads the datasets of the blocks from the file cpDataSet.bin
    void readDataSet(int step);
    void readArray(int step, Arrays arrType, std::string fname);

    //void readAverageDensityArray(int step);
    //void readAverageVelocityArray(int step);
    //void readAverageFluktuationsArray(int step);
    //void readAverageTripleArray(int step);
    //void readShearStressValArray(int step);
    //void readRelaxationFactor(int step);
    //void readPhaseField(int step, int num);
    //void readPressureField(int step);
    // 
   //! Reads the boundary conditions of the blocks from the file cpBC.bin
    void readBoundaryConds(int step);
    //! The function sets LBMKernel
    void setLBMKernel(SPtr<LBMKernel> kernel);
    //! The function sets BCSet
    void setBCSet(SPtr<BCSet> BCSet);
    //! The function truncates the data files
    void clearAllFiles(int step);

private:
    // MPI_Datatype gridParamType, block3dType;
    MPI_Datatype dataSetType, dataSetSmallType, dataSetDoubleType;
    MPI_Datatype boundCondParamType, boundCondType1000, boundCondTypeAdd, bcindexmatrixType;

    MPIIODataStructures::boundCondParam boundCondParamStr;
    SPtr<LBMKernel> lbmKernel;
    SPtr<BCSet> bcSet;

    //std::vector<double> doubleValuesArrayRW;
};

#endif
