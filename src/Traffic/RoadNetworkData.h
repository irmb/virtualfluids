#pragma once
#include <memory>
#include <vector>

#include "VectorHelper.h"
#include "invalidInput_error.h"

#include <VirtualFluidsDefinitions.h>

#include "Source.h"
#include "Sink.h"
#include "Junction.h"


struct VF_PUBLIC RoadNetworkData
{
public:
	//friend class OneWayRoad;
	//friend class OneWayRoadSSJ;

	unsigned int roadLength;
	unsigned int maxVelocity;
	unsigned int vehicleLength;
	unsigned int safetyDistance;

	std::vector<int> current;
	std::vector<int> next;					//for temporary calculations
	std::vector<int> neighbors;

	std::vector<std::shared_ptr<Sink> > sinks;
	std::vector<std::shared_ptr<Junction> > junctions;
	std::vector<std::shared_ptr<Source> > sources;
};

