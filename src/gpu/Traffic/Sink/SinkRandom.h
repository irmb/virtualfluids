#pragma once


#include <random>

#include "Sink.h"

#include "Utilities/RandomHelper.h"

class VF_PUBLIC SinkRandom:
	public Sink
{
private:
	SinkData data;

	std::mt19937 engine = RandomHelper::make_engine();
	std::uniform_real_distribution<float> distFloat{ 0.0, 1.0 };

public:
	SinkRandom(uint sinkIndex, real sinkBlockedPossibility);
	~SinkRandom() {};

	real getPossibilityBeingBlocked() const;
	bool carCanEnter();
	uint getIndex() const;
};

