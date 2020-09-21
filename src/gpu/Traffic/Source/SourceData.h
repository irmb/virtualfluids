#pragma once

#include "Core/DataTypes.h"

#include <memory>
#include <vector>
#include "Traffic_export.h"


struct TRAFFIC_EXPORT SourceData {
	uint sourceIndex;
	real sourcePossibility;
	uint maxVelocity;
};