#ifndef SINKREADER_H
#define  SINKREADER_H

#include <vector>

#include "GridGenerator_export.h"

#include "Core/DataTypes.h"
#include "Core/Logger/Logger.h"

#include "StreetPointFinder.h"

#include <VirtualFluidsDefinitions.h>

struct GRIDGENERATOR_EXPORT SinkReaderData{
	uint sinkIndex;
	float sinkBlockedPossibility;
	SinkReaderData(uint sinkIndex, float sinkBlockedPossibility);
};

struct GRIDGENERATOR_EXPORT SinkReader
{
	std::vector<SinkReaderData> sinks;
	StreetPointFinder* streetPointFinder;

	void readSinks(std::string filename, StreetPointFinder* streetPointFinder);

private:
	unsigned int getCellIndexEnd(unsigned int streetIndex);
};


#endif