#ifndef SOURCEREADER_H
#define  SOURCEREADER_H

#include <vector>

#include "Core/DataTypes.h"
#include "Core/Logger/Logger.h"

#include "StreetPointFinder.h"

#include <VirtualFluidsDefinitions.h>

struct VF_PUBLIC SourceInReader {
	unsigned int sourceIndex;
	float sourcePossibility;
	SourceInReader(unsigned int sourceIndex, float sourcePossibility);
};

struct VF_PUBLIC SourceReader
{
	std::vector<SourceInReader> sources;
	StreetPointFinder streetPointFinder;

	void readSources(std::string filename, StreetPointFinder streetPointFinder);

private:
	unsigned int getCellIndex(unsigned int streetIndex, char startOrEnd);
};


#endif