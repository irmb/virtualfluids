#ifndef  TrafficTimestep_H
#define  TrafficTimestep_H

#include <vector>
//#include <string>
#include <memory>
#include <random>
#include <thrust/device_vector.h>

#include "VirtualFluidsDefinitions.h"
#include "Utilities/RandomHelper.h"
#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"

//#include "TrafficMovement.h"

struct RoadNetworkData;
//class TrafficMovement;
class Sink;
class Junction;
class Source;

class VF_PUBLIC TrafficTimestep
{
private:

	//std::vector<real> kineticEnergyTimeSeries;
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

	thrust::device_vector<int> neighbors;

	thrust::device_vector<int> roadCurrent;
	thrust::device_vector<int> roadNext;
	thrust::device_vector<int> oldSpeeds;

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

	thrust::device_vector<bool> sinkCarCanEnter;

	thrust::device_vector<float> sourcePossibilities;
	thrust::device_vector<uint> sourceIndices;

public:

	TrafficTimestep(std::shared_ptr<RoadNetworkData> road);
	void run(std::shared_ptr<RoadNetworkData> road);
	uint getNumCarsOnJunctions(); //only used for debugging

private:
	void callTrafficMovementKernel();
	void callSourceKernel();
	void callJunctionKernel();

	void combineJuncInCellIndices(std::vector<std::shared_ptr<Junction> > &junctions);
	void combineJuncOutCellIndices(std::vector<std::shared_ptr<Junction> > &junctions);
	void combineJuncCarCanNotEnterThisOutCell(std::vector<std::shared_ptr<Junction> > &junctions);
	void initjuncOutCellIsOpen();
	void initJuncCarCanEnter();
	void initJuncCarsOnJunction();
	void initJuncAlreadyMoved();
	void initJuncOldSpeeds();

	void combineSinkCarCanEnterSink(std::vector<std::shared_ptr<Sink> > &sinks);

	void combineSourcePossibilities(std::vector<std::shared_ptr<Source> > &sources);
	void combineSourceIndices(std::vector<std::shared_ptr<Source> > &sources);

	void resetOutCellIsOpen();
};

#endif
