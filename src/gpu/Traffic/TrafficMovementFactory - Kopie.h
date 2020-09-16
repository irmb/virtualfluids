# pragma once


#include "Core/DataTypes.h"

#include <vector>
#include <memory>

#include "TrafficMovementFactory.h"
#include "GridGenerator/StreetPointFinder/StreetPointFinder.h"

#include "Traffic_export.h"


class TRAFFIC_EXPORT TrafficMovementFactoryTest :
	public TrafficMovementFactory {
public:
	TrafficMovementFactoryTest() {};
	~TrafficMovementFactoryTest() {};
	virtual void initTrafficMovement(bool useGPU, real * pConcArray = nullptr);
	virtual void calculateTimestep(uint step, uint stepForVTK);
	void loopThroughTimesteps(uint timeSteps);
};