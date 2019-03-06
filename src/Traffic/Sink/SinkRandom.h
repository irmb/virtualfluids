#pragma once

#include <VirtualFluidsDefinitions.h>

#include <iostream>
#include <random>

#include "Utilities/RandomHelper.h"
#include "Sink.h"
#include "Utilities/invalidInput_error.h"

class VF_PUBLIC SinkRandom:
	public Sink
{
private:
	SinkData data;

	std::mt19937 engine = RandomHelper::make_engine();
	std::uniform_real_distribution<float> distFloat{ 0.0, 1.0 };

public:
	SinkRandom(uint sinkIndex, float sinkBlockedPossibility);
	~SinkRandom() {};

	float getPossibilityBeingBlocked() const;
	bool carCanEnter();
	uint getIndex() const;
};

