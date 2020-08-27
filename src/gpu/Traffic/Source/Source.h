#pragma once
#include <VirtualFluidsDefinitions.h>

#include "SourceData.h"
#include "Traffic_export.h"

class TRAFFIC_EXPORT Source
{
public:
	virtual uint getIndex() const = 0;
	virtual real getPossibility() const = 0;
	virtual uint getSourceCar() = 0;
};

