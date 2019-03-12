#ifndef SINKREADER_H
#define  SINKREADER_H

#include <vector>

#include "Core/DataTypes.h"
#include "Core/Logger/Logger.h"

#include "StreetPointFinder.h"

#include <VirtualFluidsDefinitions.h>

struct VF_PUBLIC SinkInReader{
	uint sinkIndex;
	float sinkBlockedPossibility;
	SinkInReader(uint sinkIndex, float sinkBlockedPossibility);
};

struct VF_PUBLIC SinkReader
{
	std::vector<SinkInReader> sinks;
	StreetPointFinder streetPointFinder;

	void readSinks(std::string filename, StreetPointFinder streetPointFinder);

private:
	unsigned int getCellIndexEnd(unsigned int streetIndex);
};


#endif