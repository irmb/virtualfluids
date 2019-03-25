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
	int size_roads;
	uint maxVelocity;
	uint safetyDistance;

	real dawdlePossibility;
	bool useSlowToStart;
	real slowStartPossibility;
	uint maxAcceleration;

	int size_junctions;
	int size_sources;
	int size_sinks;

	thrust::device_vector<int> neighbors;

	thrust::device_vector<int> roadCurrent;
	thrust::device_vector<int> roadNext;
	thrust::device_vector<int> oldSpeeds;

	thrust::device_vector<uint> junctionInCellIndices;
	thrust::device_vector<bool> junctionCarCanEnter;
	thrust::device_vector<int> junctionCarsOnJunction;	
	thrust::device_vector<uint> junctionAlreadyMoved;
	thrust::device_vector<uint> JunctionOldSpeeds;

	thrust::device_vector<bool> sinkCarCanEnter;

	thrust::device_vector<float> sourcePossibilities;

	std::mt19937 engine = RandomHelper::make_engine();
	std::uniform_real_distribution<float> distFloat{ 0.0, 1.0 };
public:


	TrafficTimestep(std::shared_ptr<RoadNetworkData> road);
	//TrafficTimestep(RoadNetworkData* road, TrafficMovement* traffic);
	void run(std::shared_ptr<RoadNetworkData> road);

	void combineJunctionInCellIndices(std::vector<std::shared_ptr<Junction> > junctions);
	void combineJunctionCarCanEnter(std::vector<std::shared_ptr<Junction> > junctions);
	void combineJunctionCarsOnJunction(std::vector<std::shared_ptr<Junction> > junctions);
	void combineSinkCarCanEnterSink(std::vector<std::shared_ptr<Sink> > sinks);
	void combineSourcePossibilities(std::vector<std::shared_ptr<Source>> sources);
	void combineJunctionAlreadyMoved(std::vector<std::shared_ptr<Junction> > junctions);
	void combineJunctionOldSpeeds(std::vector<std::shared_ptr<Junction> > junctions);

};

#endif
