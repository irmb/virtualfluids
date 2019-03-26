# pragma once
#include "Core/DataTypes.h"

class TrafficMovementFactory {
public:
	virtual void initTrafficMovement(real * pconcArrayStart = nullptr) = 0;
	virtual void calculateTimestep(uint step, uint stepForVTK) = 0;
};
