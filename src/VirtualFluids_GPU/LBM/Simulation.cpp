#include "Simulation.h"

#include <stdio.h>
#include <vector>

#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>

#include "LBM/LB.h"
#include "Communication/Communicator.h"
#include "Parameter/Parameter.h"
#include "GPU/GPU_Interface.h"
#include "DataStructureInitializer/GridProvider.h"


#include "Output/WriteData.h"
#include "Init/InitLattice.h"

#include <utilities/logger/Logger.h>
#include <utilities/StringUtil/StringUtil.h>



Simulation::Simulation()
{

}

Simulation::~Simulation()
{

}

void Simulation::init(std::shared_ptr<Parameter> para, std::shared_ptr<GridProvider> gridProvider)
{
    this->gridProvider = gridProvider;
    comm = Communicator::getInstance();
    this->para = para;

    para->setMyID(comm->getPID());
    para->setNumprocs(comm->getNummberOfProcess());
    devCheck(comm->mapCudaDevice(para->getMyID(), para->getNumprocs(), para->getDevices(), para->getMaxDev()));

    gridProvider->allocAndCopyForcing();
    gridProvider->setDimensions();
    gridProvider->setBoundingBox();

    para->initParameter();

    para->setRe(para->getVelocity() * (doubflo)1.0 / para->getViscosity());

    para->setlimitOfNodesForVTK(30000000); //max 30 Million nodes per VTK file

    *logging::out << logging::Logger::INTERMEDIATE << "LB_Modell: D3Q" << para->getD3Qxx() << "\n";
    *logging::out << logging::Logger::INTERMEDIATE << "Re:           " << para->getRe() << "\n";
    *logging::out << logging::Logger::INTERMEDIATE << "vis_ratio:    " << para->getViscosityRatio() << "\n";
    *logging::out << logging::Logger::INTERMEDIATE << "u0_ratio:     " << para->getVelocityRatio() << "\n";
    *logging::out << logging::Logger::INTERMEDIATE << "delta_rho:    " << para->getDensityRatio() << "\n";

    para->setMemsizeGPU(0, true);
    
    gridProvider->allocArrays_CoordNeighborGeo();
    gridProvider->allocArrays_BoundaryValues();
    gridProvider->allocArrays_BoundaryQs();

    *logging::out << logging::Logger::INTERMEDIATE << "init lattice...";
    initLattice(para);
    *logging::out << logging::Logger::INTERMEDIATE << "done.\n";
    
    *logging::out << "Print files Init...";
    gridProvider->cudaCopyDataToHost(0);
    DataWriter::writeInit(para);
    *logging::out << logging::Logger::INTERMEDIATE << "done.\n";
    *logging::out << logging::Logger::INTERMEDIATE << "used Device Memory: " << para->getMemsizeGPU() / 1000000.0 << " MB\n";
}

void Simulation::run()
{
    this->logOutputHeading();
  
    CudaTimer cudaTimer;
    this->createAndStartTimer(cudaTimer);

    unsigned int timestep;
    for (timestep = para->getTStart(); timestep <= para->getTEnd(); timestep++)
        this->calculateTimestep(timestep, cudaTimer);

    this->stopLogAndDeleteTimer(cudaTimer, timestep);

    this->gridProvider->freeMemoryOnHost();
}


void Simulation::stopLogAndDeleteTimer(CudaTimer &cudaTimer, unsigned int timestep)
{
    double totalTimeSdkTimer, totalTimeEventTimer;
    float timeSizeLastTimeStep;

    cudaTimer.stopSdkTimer(timeSizeLastTimeStep, totalTimeSdkTimer);
    cudaTimer.stopEventTimer(timeSizeLastTimeStep, totalTimeEventTimer);

    this->logTotalSimulationCharacteristics(timestep, totalTimeSdkTimer);
    this->logTotalSimulationCharacteristics(timestep, totalTimeEventTimer);

    cudaTimer.deleteSdkTimer();
    cudaTimer.deleteEventTimer();
}

void Simulation::createAndStartTimer(CudaTimer &cudaTimer)
{
    cudaTimer.createSdkTimer();
    cudaTimer.createEventTimer();

    cudaTimer.startSdkTimer();
    cudaTimer.startEventTimer();
}

void Simulation::logOutputHeading()
{
    *logging::out << logging::Logger::INTERMEDIATE << "Processing time (ms) \t Nups in Mio \t throughput in GB/sec\n";
    *logging::out << logging::Logger::INTERMEDIATE << "getMaxLevel = " << para->getMaxLevel() << "\n";
}

void Simulation::calculateTimestep(unsigned int timestep, CudaTimer &cudaTimer)
{
    ////////////////////////////////////////////////////////////////////////////////
        // Collision and Propagation
        ////////////////////////////////////////////////////////////////////////////////      
        //comp
    KernelBGKSPSimple27(para->getParD(0)->numberofthreads,
        para->getParD(0)->omega,
        para->getParD(0)->geoSP,
        para->getParD(0)->neighborX_SP,
        para->getParD(0)->neighborY_SP,
        para->getParD(0)->neighborZ_SP,
        para->getParD(0)->d0SP.f[0],
        para->getParD(0)->size_Mat_SP,
        para->getParD(0)->evenOrOdd);
    getLastCudaError("KernelCasSP27 execution failed");
    //////////////////////////////////////////////////////////////////////////////
    if (para->getNumprocs() > 1)
    {
        // ...
    }
    ////////////////////////////////////////////////////////////////////////////////
    /*QVelDevComp27(para->getParD(0)->numberofthreads, para->getParD(0)->nx, para->getParD(0)->ny,
        para->getParD(0)->Qinflow.Vx, para->getParD(0)->Qinflow.Vy, para->getParD(0)->Qinflow.Vz,
        para->getParD(0)->d0SP.f[0], para->getParD(0)->Qinflow.k, para->getParD(0)->Qinflow.q27[0],
        para->getParD(0)->kInflowQ, para->getParD(0)->kInflowQ, para->getParD(0)->omega,
        para->getParD(0)->neighborX_SP, para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
        para->getParD(0)->size_Mat_SP, para->getParD(0)->evenOrOdd);
    getLastCudaError("QVelDevComp27 execution failed");*/
    ////////////////////////////////////////////////////////////////////////////////
    //QDevComp27(para->getParD(0)->numberofthreads,       para->getParD(0)->nx,           para->getParD(0)->ny,
    //		     para->getParD(0)->d0SP.f[0],             para->getParD(0)->QWall.k,		para->getParD(0)->QWall.q27[0], 
    //		     para->getParD(0)->kQ,                    para->getParD(0)->kQ,           para->getParD(0)->omega,
    //		     para->getParD(0)->neighborX_SP,          para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
    //		     para->getParD(0)->size_Mat_SP,           para->getParD(0)->evenOrOdd);
    //getLastCudaError("QDevComp27 (Wall) execution failed");
    ////////////////////////////////////////////////////////////////////////////////
    /*QDevComp27(para->getParD(0)->numberofthreads, para->getParD(0)->nx, para->getParD(0)->ny,
        para->getParD(0)->d0SP.f[0], para->getParD(0)->QGeom.k, para->getParD(0)->QGeom.q27[0],
        para->getParD(0)->QGeom.kQ, para->getParD(0)->QGeom.kQ, para->getParD(0)->omega,
        para->getParD(0)->neighborX_SP, para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
        para->getParD(0)->size_Mat_SP, para->getParD(0)->evenOrOdd);
    getLastCudaError("QDevComp27 (Geom) execution failed");*/
    ////////////////////////////////////////////////////////////////////////////////
    /*QPressDevOld27(para->getParD(0)->numberofthreads, para->getParD(0)->QPress.RhoBC,
        para->getParD(0)->d0SP.f[0], para->getParD(0)->QPress.k,
        para->getParD(0)->QPress.kN, para->getParD(0)->QPress.kQ, para->getParD(0)->omega,
        para->getParD(0)->neighborX_SP, para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
        para->getParD(0)->size_Mat_SP, para->getParD(0)->evenOrOdd);
    getLastCudaError("QPressDev27 execution failed");*/
    ////////////////////////////////////////////////////////////////////////////////
    if (para->getParD(0)->evenOrOdd == true)  para->getParD(0)->evenOrOdd = false;
    else                                    para->getParD(0)->evenOrOdd = true;
    ////////////////////////////////////////////////////////////////////////////////
    if (para->getNumprocs() > 1)
    {
        // ...
    }
    ////////////////////////////////////////////////////////////////////////////////

    this->logAndWriteResults(timestep, cudaTimer);
}

void Simulation::logAndWriteResults(unsigned int timestep, CudaTimer &cudaTimer)
{
    if (para->getTOut() > 0 && timestep % para->getTOut() == 0)
    {

        checkCudaErrors(cudaDeviceSynchronize());

        double totalTime;
        float timeSizeLastTimeStep;
        cudaTimer.stopSdkTimer(timeSizeLastTimeStep, totalTime);
        logTimeStepValues(timestep, totalTime, timeSizeLastTimeStep);

        cudaTimer.stopEventTimer(timeSizeLastTimeStep, totalTime);
        logTimeStepValues(timestep, totalTime, timeSizeLastTimeStep);

        //////////////////////////////////////////////////////////////////////////
        //exchange data for valid post process
        if (para->getNumprocs() > 1)
        {
            // ...
        }
        //////////////////////////////////////////////////////////////////////////

        if (para->getPrintFiles())
        {
            *logging::out << logging::Logger::INTERMEDIATE << "write t=" << (int)timestep << "...";
            ////////////////////////////////////////////////////////////////////////
            CalcMacCompSP27(para->getParD(0)->vx_SP,
                para->getParD(0)->vy_SP,
                para->getParD(0)->vz_SP,
                para->getParD(0)->rho_SP,
                para->getParD(0)->press_SP,
                para->getParD(0)->geoSP,
                para->getParD(0)->neighborX_SP,
                para->getParD(0)->neighborY_SP,
                para->getParD(0)->neighborZ_SP,
                para->getParD(0)->size_Mat_SP,
                para->getParD(0)->numberofthreads,
                para->getParD(0)->d0SP.f[0],
                para->getParD(0)->evenOrOdd);
            getLastCudaError("CalcMacSP27 execution failed");
            ////////////////////////////////////////////////////////////////////////
            gridProvider->cudaCopyDataToHost(0);
            ////////////////////////////////////////////////////////////////////////
            DataWriter::writeTimestep(para, timestep);
            ////////////////////////////////////////////////////////////////////////
            *logging::out << logging::Logger::INTERMEDIATE << "done.\n";
        }

        cudaTimer.startSdkTimer();
        cudaTimer.startEventTimer();
    }
}

void Simulation::logTimeStepValues(unsigned int timestep, double totalTime, float timeSizeLastTimeStep)
{
    double  fnups, throughput;
    fnups = 0.0;
    throughput = 0.0;
    fnups += 1000.0 * (timestep - para->getTStart()) * para->getParH(0)->size_Mat_SP / (totalTime*1.0E6);
    throughput += (27.0 + 1.0) * 4.0 * 1000.0 * (timestep - para->getTStart()) * para->getParH(0)->size_Mat_SP / (totalTime*1.0E9);
    *logging::out << logging::Logger::INTERMEDIATE << timeSizeLastTimeStep << " / " << totalTime << " \t " << fnups << " \t " << throughput << "\n";
}

void Simulation::logTotalSimulationCharacteristics(unsigned int timestep, double totalTime)
{
    double fnups, throughput;
    fnups = 0.0;
    throughput = 0.0;
    fnups += 1000.0 * (timestep - para->getTStart()) * para->getParH(0)->size_Mat_SP / (totalTime*1.0E6);
    throughput += (27.0 + 1.0) * 4.0 * 1000.0 * (timestep - para->getTStart()) * para->getParH(0)->size_Mat_SP / (totalTime*1.0E9);
    *logging::out << logging::Logger::INTERMEDIATE << "Processing time: " << totalTime << "(ms)\n";
    *logging::out << logging::Logger::INTERMEDIATE << "Nups in Mio: " << fnups << "\n";
    *logging::out << logging::Logger::INTERMEDIATE << "throughput in GB/sec: " << throughput << "\n";
}

