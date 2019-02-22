#pragma once

#include "SinkData.h"
#include <VirtualFluidsDefinitions.h>

class VF_PUBLIC Sink
{
public:
	virtual float getPossibilityBeingBlocked() const = 0;
	virtual bool carCanEnter() = 0;
	virtual unsigned int getIndex() const = 0;
};

