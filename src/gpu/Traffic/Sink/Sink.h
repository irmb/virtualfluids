#pragma once

#include "SinkData.h"

class VIRTUALFLUIDS_GPU_EXPORT Sink
{
public:
	virtual real getPossibilityBeingBlocked() const = 0;
	virtual bool carCanEnter() = 0;
	virtual uint getIndex() const = 0;
};

