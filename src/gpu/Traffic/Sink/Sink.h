#pragma once

#include "SinkData.h"

#include "Traffic_export.h"

class TRAFFIC_EXPORT Sink
{
public:
	virtual real getPossibilityBeingBlocked() const = 0;
	virtual bool carCanEnter() = 0;
	virtual uint getIndex() const = 0;
};

