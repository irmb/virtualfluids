#pragma once
#include <VirtualFluidsDefinitions.h>

#include "SourceData.h"

class VF_PUBLIC Source
{
public:
	virtual uint getIndex() const = 0;
	virtual real getPossibility() const = 0;
	virtual uint getSourceCar() = 0;
};

