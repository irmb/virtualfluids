#pragma once
#include <VirtualFluidsDefinitions.h>
#include "Core/DataTypes.h"

#include <memory>
#include <vector>


struct VF_PUBLIC SourceData {
	uint sourceIndex;
	float sourcePossibility;
	uint maxVelocity;
};