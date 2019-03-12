#pragma once

#include "SinkData.h"
#include <VirtualFluidsDefinitions.h>

class VF_PUBLIC Sink
{
public:
	virtual real getPossibilityBeingBlocked() const = 0;
	virtual bool carCanEnter() = 0;
	virtual uint getIndex() const = 0;
};

