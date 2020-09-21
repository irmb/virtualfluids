#pragma once
#include <memory>
#include <vector>



#include "Source/Source.h"
#include "Sink/Sink.h"
#include "Junction/Junction.h"

#include "Traffic_export.h"


struct TRAFFIC_EXPORT RoadNetworkData
{
protected:
	friend class TrafficMovement;
	friend class TrafficTimestep;

	uint roadLength;
	uint maxVelocity;
	uint vehicleLength;
	uint safetyDistance;

	std::vector<int> current;
	std::vector<int> currentWithLongVehicles;
	std::vector<int> next;						//for temporary calculations
	std::vector<int> neighbors;

	std::vector<std::shared_ptr<Sink> > sinks; 
	std::vector<std::shared_ptr<Junction> > junctions;
	std::vector<std::shared_ptr<Source> > sources;

	std::vector<int> *pcurrent;
	std::vector<int> *pnext;
	std::vector<int> *pdummy;

	std::vector<int> oldSpeeds;

	real dawdlePossibility;
	bool useSlowToStart = false;
	real slowStartPossibility;
	uint maxAcceleration = 1;

	std::vector<real> conc; //dispConcFromGPU
};

