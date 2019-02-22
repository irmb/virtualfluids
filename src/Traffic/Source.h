#pragma once
#include <VirtualFluidsDefinitions.h>

#include "SourceData.h"

class VF_PUBLIC Source
{
public:
	virtual unsigned int getIndex() const = 0;
	virtual float getPossibility() const = 0;
	virtual unsigned int getSourceCar() = 0;
};

