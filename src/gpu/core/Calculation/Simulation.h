//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \author Martin Schoenherr
//=======================================================================================
#ifndef GPU_SIMULATION_H_
#define GPU_SIMULATION_H_

#include <memory>
#include <vector>

#include "PointerDefinitions.h"
#include "Calculation/Calculation.h"

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
