#ifndef _SIMULATION_H_
#define _SIMULATION_H_
#include "VirtualFluids_GPU_EXPORT.h"

#include <memory>
#include <string>


#include "CudaTimer.h"


class GridProvider;
class Parameter;
class Communicator;

class Simulation
{
public:
    VirtualFluids_GPU_EXPORT Simulation();
    VirtualFluids_GPU_EXPORT ~Simulation();
    VirtualFluids_GPU_EXPORT void run();
    VirtualFluids_GPU_EXPORT void init(std::shared_ptr<Parameter> para, std::shared_ptr<GridProvider> gridReader);

protected:
    std::shared_ptr<Communicator> comm;
    std::shared_ptr<Parameter> para;
    std::shared_ptr<GridProvider> gridProvider;

private:
    void logOutputHeading();
    void createAndStartTimer(CudaTimer &cudaTimer);
    void calculateTimestep(unsigned int timestep, CudaTimer &cudaTimer);

    void logAndWriteResults(unsigned int timestep, CudaTimer &cudaTimer);

    void logTimeStepValues(unsigned int timestep, double totalTime, float timeSizeLastTimeStep);

    void logTotalSimulationCharacteristics(unsigned int timestep, double totalTime);
    void stopLogAndDeleteTimer(CudaTimer &cudaTimer, unsigned int timestep);

};
#endif
