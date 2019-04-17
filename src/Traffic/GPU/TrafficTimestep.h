#ifndef  TrafficTimestep_H
#define  TrafficTimestep_H

#include <vector>
#include <memory>

#include <thrust/device_vector.h>
#include <curand_kernel.h>

#include <VirtualFluidsDefinitions.h>
#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"

struct RoadNetworkData;
class Sink;
class Junction;
class Source;

class VF_PUBLIC TrafficTimestep
{
private:

	bool timestepIsEven = true;
	uint numTimestep = 0;


	uint maxVelocity;
	uint safetyDistance;
	real dawdlePossibility;
	bool useSlowToStart;
	real slowStartPossibility;
	uint maxAcceleration;


	//sizes
	uint size_roads;
	uint size_junctions;
	uint size_juncInCells;
	uint size_juncOutCells;
	uint size_sources;
	uint size_sinks;


	//grids
	dim3 gridRoad;
	dim3 threadsRoads;
	dim3 gridJunctions;
	dim3 threadsJunctions;
	dim3 gridSources;
	dim3 threadsSources;


	//road
	thrust::device_vector<int> neighbors;
	thrust::device_vector<int> roadCurrent;
	thrust::device_vector<int> roadNext;


	//junctions
	thrust::device_vector<uint> juncTrafficLightSwitchTime; //no TrafficLight: 0

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


	//sinks
	thrust::device_vector<real> sinkCarBlockedPossibilities;;


	//sources
	thrust::device_vector<float> sourcePossibilities;
	thrust::device_vector<uint> sourceIndices;


	//concentrations
	real * pConcArray;


	//curandStates
	curandState *statesJunctions;
	curandState *statesSources;
	curandState *statesRoad;

	int* pRoadCurrent;
	int* pRoadNext;

public:

	TrafficTimestep(std::shared_ptr<RoadNetworkData> road, real * pConcArray);
	void calculateTimestep(std::shared_ptr<RoadNetworkData> road);
	void cleanUp();
	uint getNumCarsOnJunctions(); //only used for debugging
	void copyCurrentDeviceToHost(std::shared_ptr<RoadNetworkData> road);

private:
	
	//timestep
	void switchCurrentNext();
	void resetOutCellIsOpen();
	void resetNext();

	//kernel calls
	void callTrafficTimestepKernel();
	void callSourceTimestepKernel();
	void callJunctionTimestepKernel();

	//init grids
	void calculateTrafficTimestepKernelDimensions();
	void calculateJunctionKernelDimensions();
	void calculateSourceKernelDimensions();

	//init junctions
	void combineJuncInCellIndices(std::vector<std::shared_ptr<Junction> > &junctions);
	void combineJuncOutCellIndices(std::vector<std::shared_ptr<Junction> > &junctions);
	void combineJuncCarCanNotEnterThisOutCell(std::vector<std::shared_ptr<Junction> > &junctions);
	void combineUseTrafficLights(std::vector<std::shared_ptr<Junction>>& junctions);
	void initjuncOutCellIsOpen();
	void initJuncCarCanEnter();
	void initJuncCarsOnJunction();
	void initJuncAlreadyMoved();
	void initJuncOldSpeeds();
	
	//init sinks
	void combineSinkBlockedPossibilities(std::vector<std::shared_ptr<Sink>>& sinks);

	//init sources
	void combineSourcePossibilities(std::vector<std::shared_ptr<Source> > &sources);
	void combineSourceIndices(std::vector<std::shared_ptr<Source> > &sources);
};

#endif
