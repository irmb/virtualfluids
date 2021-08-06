#ifndef _MPIIOMigrationCoProcessor_H_
#define _MPIIOMigrationCoProcessor_H_

#include <mpi.h>
#include <string>

#include "MPIIOCoProcessor.h"
#include "MPIIODataStructures.h"

class Grid3D;
class UbScheduler;
class Communicator;
class BCProcessor;
class LBMKernel;
class Grid3DVisitor;

//! \class MPIWriteBlocksCoProcessor
//! \brief Writes the grid each timestep into the files and reads the grip from the files before regenerating
class MPIIOMigrationCoProcessor : public MPIIOCoProcessor
{
public:
    enum Arrays {
        AverageDensity      = 1,
        AverageVelocity     = 2,
        AverageFluktuations = 3,
        AverageTriple       = 4,
        ShearStressVal      = 5,
        RelaxationFactor = 6,
        PhaseField1 = 7,
        PhaseField2 = 8
    };

    MPIIOMigrationCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<Grid3DVisitor> mV, const std::string &path, SPtr<Communicator> comm);
    ~MPIIOMigrationCoProcessor() override;
    //! Each timestep writes the grid into the files
    void process(double step) override;
    //! Reads the grid from the files before grid reconstruction
    void restart(int step);
    //! Writes the blocks of the grid into the file cpBlocks.bin
    void writeBlocks(int step);
    //! Writes the datasets of the blocks into the file cpDataSet.bin
    void writeDataSet(int step);
    void write4DArray(int step, Arrays arrType, std::string fname);
    void write3DArray(int step, Arrays arrType, std::string fname);
    //   void writeAverageDensityArray(int step);
    //   void writeAverageVelocityArray(int step);
    //   void writeAverageFluktuationsArray(int step);
    //   void writeAverageTripleArray(int step);
    //   void writeShearStressValArray(int step);
    //   void writeRelaxationFactor(int step);
    //! Writes the boundary conditions of the blocks into the file cpBC.bin
    void writeBoundaryConds(int step);

    //! Reads the blocks of the grid from the file cpBlocks.bin
    void readBlocks(int step);
    //! Reads the datasets of the blocks from the file cpDataSet.bin
    void readDataSet(int step);
    void readArray(int step, Arrays arrType, std::string fname);
    //   void readAverageDensityArray(int step);
    //   void readAverageVelocityArray(int step);
    //   void readAverageFluktuationsArray(int step);
    //   void readAverageTripleArray(int step);
    //   void readShearStressValArray(int step);
    //   void readRelaxationFactor(int step);
    //! Reads the boundary conditions of the blocks from the file cpBC.bin
    void readBoundaryConds(int step);
    //! The function sets LBMKernel
    void setLBMKernel(SPtr<LBMKernel> kernel);
    //! The function sets BCProcessor
    void setBCProcessor(SPtr<BCProcessor> bcProcessor);
    //! The function truncates the data files
    void clearAllFiles(int step);
    // void setNu(double nu);

protected:
    // std::string path;
    // SPtr<Communicator> comm;

private:
    // MPI_Datatype gridParamType, block3dType;
    MPI_Datatype dataSetType, dataSetSmallType, dataSetDoubleType;
    MPI_Datatype boundCondParamType, boundCondTypeAdd, bcindexmatrixType;

    MPIIODataStructures::boundCondParam boundCondParamStr;
    SPtr<LBMKernel> lbmKernel;
    SPtr<BCProcessor> bcProcessor;
    // double nue;
    SPtr<Grid3DVisitor> metisVisitor;
};

#endif
