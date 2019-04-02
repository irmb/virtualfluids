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

	uint trafficLightSwitchTime; //no TrafficLight: 0

	std::mt19937 engine = RandomHelper::make_engine();
	std::uniform_real_distribution<float> distFloat{ 0.0, 1.0 };
};

