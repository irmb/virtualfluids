#pragma once
#include <memory>
#include <vector>

#include <VirtualFluidsDefinitions.h>


#include "Source/Source.h"
#include "Sink/Sink.h"
#include "Junction/Junction.h"


struct VF_PUBLIC RoadNetworkData
{
public:
	friend class TrafficMovement;

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
};

