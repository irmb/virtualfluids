#ifndef JUNCTIONREADER_H
#define JUNCTIONREADER_H

#include <vector>

#include "GridGenerator_export.h"

#include "Core/DataTypes.h"
#include "Core/Logger/Logger.h"

#include "StreetPointFinder.h"



struct GRIDGENERATOR_EXPORT JunctionReaderData
{
	std::vector<uint> inCells;
	std::vector<uint> outCells;
	std::vector<int> carCanNotEnterThisOutCell;
	uint trafficLightSwitchTime;

	JunctionReaderData(std::vector<uint> inCells, std::vector<uint> outCells, std::vector<int> carCanNotEnterThisOutCell, uint trafficLightSwitchTime);
};


struct GRIDGENERATOR_EXPORT Neighbors
{
	std::vector<int> cells;
	std::vector<int> neighbors;
};



struct GRIDGENERATOR_EXPORT JunctionReader
{
	std::vector<JunctionReaderData> junctions;
	Neighbors specialNeighbors;
	StreetPointFinder* streetPointFinder;

	void readJunctions(std::string filename, StreetPointFinder* streetPointFinder);


private:
	unsigned int getCellIndex(unsigned int streetIndex, char startOrEnd);
};
#endif
