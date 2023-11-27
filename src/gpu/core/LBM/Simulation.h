#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <memory>
#include <vector>

#include <PointerDefinitions.h>

#include "LBM/LB.h"
#include "Utilities/Buffer2D.hpp"

namespace vf::parallel
{
class Communicator;
}

class CudaMemoryManager;
class Parameter;
class GridProvider;
class RestartObject;
class ForceCalculations;
class DataWriter;
class Kernel;
class AdvectionDiffusionKernel;
class KernelFactory;
class PreProcessor;
class PreProcessorFactory;
class TrafficMovementFactory;
class UpdateGrid27;
class KineticEnergyAnalyzer;
class EnstrophyAnalyzer;
class BoundaryConditionFactory;
class GridScalingFactory;
class TurbulenceModelFactory;
class Timer;

class Simulation
{
public:
    Simulation(std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> memoryManager,
               vf::parallel::Communicator &communicator, GridProvider &gridProvider, BoundaryConditionFactory* bcFactory, GridScalingFactory* scalingFactory = nullptr);    
    Simulation(std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> memoryManager,
               vf::parallel::Communicator &communicator, GridProvider &gridProvider, BoundaryConditionFactory* bcFactory, SPtr<TurbulenceModelFactory> tmFactory, GridScalingFactory* scalingFactory = nullptr);

    ~Simulation();
    void run();

    void setFactories(std::unique_ptr<KernelFactory> &&kernelFactory,
                      std::unique_ptr<PreProcessorFactory> &&preProcessorFactory);
    void setDataWriter(std::shared_ptr<DataWriter> dataWriter);
    void addKineticEnergyAnalyzer(uint tAnalyse);
    void addEnstrophyAnalyzer(uint tAnalyse);

    //! \brief can be used as an alternative to run(), if the simulation needs to be controlled from the outside (e. g. for fluid structure interaction FSI)
    void calculateTimestep(uint timestep);
    //! \brief needed to initialize the simulation timers if calculateTimestep is used instead of run()
    void initTimers();

private:
    void init(GridProvider &gridProvider, BoundaryConditionFactory *bcFactory, SPtr<TurbulenceModelFactory> tmFactory, GridScalingFactory *scalingFactory);
    void allocNeighborsOffsetsScalesAndBoundaries(GridProvider& gridProvider);
    void readAndWriteFiles(uint timestep);

    std::unique_ptr<KernelFactory> kernelFactory;
    std::shared_ptr<PreProcessorFactory> preProcessorFactory;

    Buffer2D <real> sbuf_t;
    Buffer2D <real> rbuf_t;
    Buffer2D <real> sbuf_b;
    Buffer2D <real> rbuf_b;

    Buffer2D <int> geo_sbuf_t;
    Buffer2D <int> geo_rbuf_t;
    Buffer2D <int> geo_sbuf_b;
    Buffer2D <int> geo_rbuf_b;

    vf::parallel::Communicator& communicator;
    SPtr<Parameter> para;
    std::shared_ptr<DataWriter> dataWriter;
    std::shared_ptr<CudaMemoryManager> cudaMemoryManager;
    std::vector < SPtr< Kernel>> kernels;
    std::vector < SPtr< AdvectionDiffusionKernel>> adKernels;
    std::shared_ptr<PreProcessor> preProcessor;
    std::shared_ptr<PreProcessor> preProcessorAD;
    SPtr<TurbulenceModelFactory> tmFactory;

    SPtr<RestartObject> restart_object;

    // Timer
    std::unique_ptr<Timer> averageTimer;
    uint previousTimestepForAveraging;
    uint previousTimestepForTurbulenceIntensityCalculation;
    uint timestepForMeasuringPoints;

    //Forcing Calculation
    std::shared_ptr<ForceCalculations> forceCalculator;

    std::unique_ptr<KineticEnergyAnalyzer> kineticEnergyAnalyzer;
    std::unique_ptr<EnstrophyAnalyzer> enstrophyAnalyzer;
    std::unique_ptr<UpdateGrid27> updateGrid27;
};
#endif
