#pragma once
#include <VirtualFluidsDefinitions.h>

#include <vector>
#include <memory>

struct VF_PUBLIC JunctionData
{
public:
	std::vector<unsigned int> inCellIndices;
	std::vector<unsigned int> outCellIndices;
};

