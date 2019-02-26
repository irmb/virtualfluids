#pragma once
#include <VirtualFluidsDefinitions.h>

#include <memory>
#include <vector>


struct VF_PUBLIC SourceData {
	unsigned int sourceIndex;
	float sourcePossibility;
	unsigned int maxVelocity;
};