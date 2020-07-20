#ifndef SINKREADER_H
#define  SINKREADER_H

#include <vector>

#include "Core/DataTypes.h"
#include "Core/Logger/Logger.h"

#include "StreetPointFinder.h"

#include <VirtualFluidsDefinitions.h>

struct VIRTUALFLUIDS_GPU_EXPORT SinkReaderData{
	uint sinkIndex;
	float sinkBlockedPossibility;
	SinkReaderData(uint sinkIndex, float sinkBlockedPossibility);
};

struct VIRTUALFLUIDS_GPU_EXPORT SinkReader
{
	std::vector<SinkReaderData> sinks;
	StreetPointFinder* streetPointFinder;

	void readSinks(std::string filename, StreetPointFinder* streetPointFinder);

private:
	unsigned int getCellIndexEnd(unsigned int streetIndex);
};


#endif