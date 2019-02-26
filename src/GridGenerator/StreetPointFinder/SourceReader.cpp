#include "SourceReader.h"

#include <fstream>
#include <iostream>

SourceInReader::SourceInReader(unsigned int sourceIndex, float sourcePossibility):
	sourceIndex{sourceIndex}, sourcePossibility{sourcePossibility}
{}


void SourceReader::readSources(std::string filename, StreetPointFinder streetPointFinder)
{
	*logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::readSources( " << filename << " )" << "\n";

	this->streetPointFinder = streetPointFinder;

	std::ifstream file;
	file.open(filename.c_str());
	if (!file.is_open()) std::cerr << "File not found" << std::endl;

	uint numberOfSources;
	file >> numberOfSources;

	uint streetIndex;
	char startOrEnd;
	float sourcePossibility;


	for (uint i = 0; i < numberOfSources; i++) {
		file >> streetIndex >> startOrEnd >> sourcePossibility;
		sources.push_back(SourceInReader(getCellIndex(streetIndex, startOrEnd), sourcePossibility));
	}
}


unsigned int SourceReader::getCellIndex(unsigned int streetIndex, char startOrEnd)
{
	uint i = 0;
	unsigned int cellIndex = 0;
	while (i < streetIndex) {
		cellIndex += streetPointFinder.streets[i].numberOfCells;
		++i;
	}
	if (startOrEnd == 's') {
		return cellIndex;
	}
	return cellIndex + streetPointFinder.streets[streetIndex].numberOfCells - 1;
}
