#pragma once
#include <VirtualFluidsDefinitions.h>

#include <vector>
#include <memory>

#include "Utilities/RandomHelper.h"


struct VF_PUBLIC JunctionData
{
public:
	std::vector<uint> inCellIndices;
	std::vector<uint> outCellIndices; 
	std::vector<int> carCanNotEnterThisOutCell; //no such inCell: -2

	std::vector<uint> possibleOutCells;

	std::vector<bool> carCanEnter;
	std::vector<int> carsOnJunction;
	std::vector<uint> alreadyMoved;
	std::vector<uint> oldSpeeds;

	std::mt19937 engine = RandomHelper::make_engine();
	std::uniform_int_distribution<uint> distInt2{ 0, 1 };
	std::uniform_int_distribution<uint> distInt3{ 0, 2 };
	std::uniform_int_distribution<uint> distInt4{ 0, 3 };
};

