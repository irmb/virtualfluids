#ifndef SOURCEREADER_H
#define  SOURCEREADER_H

#include <vector>

#include "Core/DataTypes.h"
#include "Core/Logger/Logger.h"

#include "StreetPointFinder.h"

#include <VirtualFluidsDefinitions.h>

struct VF_PUBLIC SourceReaderData {
	unsigned int sourceIndex;
	float sourcePossibility;
	SourceReaderData(unsigned int sourceIndex, float sourcePossibility);
};

struct VF_PUBLIC SourceReader
{
	std::vector<SourceReaderData> sources;
	StreetPointFinder streetPointFinder;

	void readSources(std::string filename, StreetPointFinder streetPointFinder);

private:
	unsigned int getCellIndexStart(unsigned int streetIndex);
};


#endif