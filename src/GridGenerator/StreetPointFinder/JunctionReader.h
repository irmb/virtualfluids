#ifndef JUNCTIOREADER_H
#define JUNCTIONREADER_H

#include <vector>

#include "Core/DataTypes.h"
#include "Core/Logger/Logger.h"

#include "StreetPointFinder.h"

#include <VirtualFluidsDefinitions.h>

struct VF_PUBLIC JunctionInReader
{
	std::vector<uint> inCells;
	std::vector<uint> outCells;
	std::vector<int> carCanNotEnterThisOutCell;
	uint trafficLightSwitchTime;

	JunctionInReader(std::vector<uint> inCells, std::vector<uint> outCells, std::vector<int> carCanNotEnterThisOutCell, uint trafficLightSwitchTime);
};


struct VF_PUBLIC JunctionReader
{
	std::vector<JunctionInReader> junctions;
	StreetPointFinder streetPointFinder;

	void readJunctions(std::string filename, StreetPointFinder streetPointFinder);


private:
	unsigned int getCellIndex(unsigned int streetIndex, char startOrEnd);
};
#endif
