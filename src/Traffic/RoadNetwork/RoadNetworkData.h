#pragma once
#include <memory>
#include <vector>

#include "Utilities/VectorHelper.h"
#include "Utilities/invalidInput_error.h"

#include <VirtualFluidsDefinitions.h>

#include "Source/Source.h"
#include "Sink/Sink.h"
#include "Junction/Junction.h"


struct VF_PUBLIC RoadNetworkData
{
public:
	//friend class TrafficMovement;

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
};

