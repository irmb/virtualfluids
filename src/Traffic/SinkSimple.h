#pragma once

#include <VirtualFluidsDefinitions.h>

#include <iostream>
#include <random>

#include "RandomHelper.h"
#include "Sink.h"
#include "invalidInput_error.h"

using namespace std;

class VF_PUBLIC SinkSimple:
	public Sink
{
private:
	SinkData data;

	mt19937 engine = RandomHelper::make_engine();
	uniform_real_distribution<float> distFloat{ 0.0, 1.0 };

public:
	SinkSimple(unsigned int sinkIndex, float sinkBlockedPossibility);
	~SinkSimple() {};

	float getPossibilityBeingBlocked() const;
	bool carCanEnter();
	unsigned int getIndex() const;
};

