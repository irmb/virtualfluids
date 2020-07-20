# pragma once

#include <VirtualFluidsDefinitions.h>

#include <vector>
#include <memory>

#include "TrafficMovementFactory.h"
#include "Core/DataTypes.h"
#include "GridGenerator/StreetPointFinder/StreetPointFinder.h"


class VIRTUALFLUIDS_GPU_EXPORT TrafficMovementFactoryTest :
	public TrafficMovementFactory {
public:
	TrafficMovementFactoryTest() {};
	~TrafficMovementFactoryTest() {};
	virtual void initTrafficMovement(real * pConcArray = nullptr);
	virtual void calculateTimestep(uint step, uint stepForVTK);
	void loopThroughTimesteps(uint timeSteps);
};