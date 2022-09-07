#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <memory>
#include <vector>

#include <PointerDefinitions.h>

#include "Utilities/Buffer2D.hpp"
#include "LBM/LB.h"


namespace vf::gpu { class Communicator; }

class CudaMemoryManager;
class Parameter;
class GridProvider;
class PorousMedia;
class RestartObject;
class ForceCalculations;
class DataWriter;
class Kernel;
class ADKernel;
class KernelFactory;
class PreProcessor;
class PreProcessorFactory;
class TrafficMovementFactory;
class UpdateGrid27;
class KineticEnergyAnalyzer;
class EnstrophyAnalyzer;
class BoundaryConditionFactory;

class Simulation
{
public:
    Simulation(std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> memoryManager,
               vf::gpu::Communicator &communicator, GridProvider &gridProvider, BoundaryConditionFactory* bcFactory);

    ~Simulation();
    void run();

    void setFactories(std::unique_ptr<KernelFactory> &&kernelFactory,
               std::unique_ptr<PreProcessorFactory> &&preProcessorFactory);
    void setDataWriter(std::shared_ptr<DataWriter> dataWriter);
    void addKineticEnergyAnalyzer(uint tAnalyse);
    void addEnstrophyAnalyzer(uint tAnalyse);

private:
	void init(GridProvider &gridProvider, BoundaryConditionFactory *bcFactory);
    void allocNeighborsOffsetsScalesAndBoundaries(GridProvider& gridProvider);
    void porousMedia();
    void definePMarea(std::shared_ptr<PorousMedia>& pm);

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


	vf::gpu::Communicator& communicator;
    SPtr<Parameter> para;
    std::shared_ptr<DataWriter> dataWriter;
	std::shared_ptr<CudaMemoryManager> cudaMemoryManager;
	std::vector < SPtr< Kernel>> kernels;
	std::vector < SPtr< ADKernel>> adKernels;
	std::shared_ptr<PreProcessor> preProcessor;

    SPtr<RestartObject> restart_object;

	//Forcing Calculation
	std::shared_ptr<ForceCalculations> forceCalculator;

	//Porous Media
	std::vector<std::shared_ptr<PorousMedia>> pm;
	//PorousMedia* pm0;
	//PorousMedia* pm1;
	//PorousMedia* pm2;

    // TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/29
	//KQ - Schlaff
	// unsigned int            kNQ, kSQ, kEQ, kWQ;
	// QforBoundaryConditions  QnH, QnD;
	// QforBoundaryConditions  QsH, QsD;
	// QforBoundaryConditions  QeH, QeD;
	// QforBoundaryConditions  QwH, QwD;
	// real *VxNH,          *VyNH,       *VzNH,       *deltaVNH;
	// real *VxND,          *VyND,       *VzND,       *deltaVND;
	// real *VxSH,          *VySH,       *VzSH,       *deltaVSH;
	// real *VxSD,          *VySD,       *VzSD,       *deltaVSD;
	// real *VxEH,          *VyEH,       *VzEH,       *deltaVEH;
	// real *VxED,          *VyED,       *VzED,       *deltaVED;
	// real *VxWH,          *VyWH,       *VzWH,       *deltaVWH;
	// real *VxWD,          *VyWD,       *VzWD,       *deltaVWD;


    std::unique_ptr<KineticEnergyAnalyzer> kineticEnergyAnalyzer;
    std::unique_ptr<EnstrophyAnalyzer> enstrophyAnalyzer;
    std::unique_ptr<UpdateGrid27> updateGrid27;
};
#endif
