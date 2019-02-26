#pragma once
#include <VirtualFluidsDefinitions.h>

#include <vector>
#include <memory>
#include "Utilities\RandomHelper.h"

struct VF_PUBLIC JunctionData
{
public:
	std::vector<unsigned int> inCellIndices;
	std::vector<unsigned int> outCellIndices; 
	std::vector<int> carCanNotEnterThisOutCell; //no such inCell: -2

	std::vector<unsigned int> freeOutCells;

	std::vector<bool> carCanEnter;
	std::vector<int> carsOnJunction;
	std::vector<unsigned int> alreadyMoved;

	std::mt19937 engine = RandomHelper::make_engine();
	uniform_int_distribution<unsigned int> distInt2{ 0, 1 };
	uniform_int_distribution<unsigned int> distInt3{ 0, 2 };
	uniform_int_distribution<unsigned int> distInt4{ 0, 3 };
};

