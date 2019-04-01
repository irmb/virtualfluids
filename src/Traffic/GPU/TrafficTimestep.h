#ifndef  TrafficTimestep_H
#define  TrafficTimestep_H

#include <vector>
#include <memory>
#include <random>
#include <thrust/device_vector.h>

#include <VirtualFluidsDefinitions.h>
#include "Utilities/RandomHelper.h"
#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include <curand_kernel.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

struct RoadNetworkData;
class Sink;
class Junction;
class Source;

class VF_PUBLIC TrafficTimestep
{
private:
	bool timestepIsEven = true;

	uint size_roads;
	uint maxVelocity;
	uint safetyDistance;
	real dawdlePossibility;
	bool useSlowToStart;
	real slowStartPossibility;
	uint maxAcceleration;

	uint size_junctions;
	uint size_juncInCells;
	uint size_juncOutCells;
	uint size_sources;
	uint size_sinks;

	dim3 gridRoad;
	dim3 threadsRoads;
	dim3 gridJunctions;
	dim3 threadsJunctions;
	dim3 gridSources;
	dim3 threadsSources;

	thrust::device_vector<int> neighbors;

	thrust::device_vector<int> roadCurrent;
	thrust::device_vector<int> roadNext;

	thrust::device_vector<uint> juncInCellIndices;
	thrust::device_vector<bool> juncCarCanEnter;
	thrust::device_vector<int> juncCarsOnJunction;	
	thrust::device_vector<uint> juncAlreadyMoved;
	thrust::device_vector<uint> juncOldSpeeds;

	thrust::device_vector<int> juncOutCellIndices;
	thrust::device_vector<int> juncCarCanNotEnterThisOutCell; //no such inCell: -2
	thrust::device_vector<bool> juncOutCellIsOpen;
	thrust::device_vector<uint> juncStartInIncells;
	thrust::device_vector<uint> juncStartInOutcells;

	thrust::device_vector<real> sinkCarBlockedPossibilities;;

	thrust::device_vector<float> sourcePossibilities;
	thrust::device_vector<uint> sourceIndices;

	curandState *statesJunctions;
	curandState *statesSources;
	curandState *statesRoad;

	real * pConcArray;

public:

	TrafficTimestep(std::shared_ptr<RoadNetworkData> road, real * pConcArray);
	void run(std::shared_ptr<RoadNetworkData> road);
	void cleanUp();
	uint getNumCarsOnJunctions(); //only used for debugging
	void copyCurrentDeviceToHost(std::shared_ptr<RoadNetworkData> road);

private:
;
	//timestep
	void switchCurrentNext();
	void resetOutCellIsOpen();

	void callTrafficTimestepKernel();
	void callSourceKernel();
	void callJunctionKernel();

	//init
	void calculateTrafficTimestepKernelDimensions();
	void calculateJunctionKernelDimensions();
	void calculateSourceKernelDimensions();

	void combineJuncInCellIndices(std::vector<std::shared_ptr<Junction> > &junctions);
	void combineJuncOutCellIndices(std::vector<std::shared_ptr<Junction> > &junctions);
	void combineJuncCarCanNotEnterThisOutCell(std::vector<std::shared_ptr<Junction> > &junctions);
	void initjuncOutCellIsOpen();
	void initJuncCarCanEnter();
	void initJuncCarsOnJunction();
	void initJuncAlreadyMoved();
	void initJuncOldSpeeds();

	void combineSinkBlockedPossibilities(std::vector<std::shared_ptr<Sink>>& sinks);

	void combineSourcePossibilities(std::vector<std::shared_ptr<Source> > &sources);
	void combineSourceIndices(std::vector<std::shared_ptr<Source> > &sources);


	
};

#endif
