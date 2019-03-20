# pragma once

//#include <VirtualFluidsDefinitions.h>
#include "TrafficMovementFactoryImpl.h"

#include <vector>
#include <memory>

#include "Core/DataTypes.h"
#include "GridGenerator/StreetPointFinder/StreetPointFinder.h"


class VF_PUBLIC TrafficMovementFactoryTest:
public TrafficMovementFactoryImpl{
public:
	TrafficMovementFactoryTest() {};
	~TrafficMovementFactoryTest() {};
	virtual void initTrafficMovement(real * pconcArrayStart = nullptr);
	virtual void calculateTimestep(uint step, uint stepForVTK);
	void loopThroughTimesteps(uint timeSteps);
};