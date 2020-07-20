#pragma once
#include <VirtualFluidsDefinitions.h>
#include "Core/DataTypes.h"

#include <memory>
#include <vector>


struct VIRTUALFLUIDS_GPU_EXPORT SourceData {
	uint sourceIndex;
	real sourcePossibility;
	uint maxVelocity;
};