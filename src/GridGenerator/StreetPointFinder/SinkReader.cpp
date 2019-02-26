#include "SinkReader.h"

#include <fstream>
#include <iostream>

SinkInReader::SinkInReader(uint sinkIndex, float sinkBlockedPossibility) :
	sinkIndex{ sinkIndex }, sinkBlockedPossibility{ sinkBlockedPossibility }
{}

void SinkReader::readSinks(std::string filename, StreetPointFinder streetPointFinder)
{
	*logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::readSinks( " << filename << " )" << "\n";

	this->streetPointFinder = streetPointFinder;

	std::ifstream file;
	file.open(filename.c_str());
	if (!file.is_open()) std::cerr << "File not found" << std::endl;

	uint numberOfSinks;
	file >> numberOfSinks;

	uint streetIndex;
	char startOrEnd;
	float sinkBlockedPossibility;


	for (uint i = 0; i < numberOfSinks; i++) {
		file >> streetIndex >> startOrEnd >> sinkBlockedPossibility;
		sinks.push_back(SinkInReader(getCellIndex(streetIndex, startOrEnd), sinkBlockedPossibility));
	}
}

unsigned int SinkReader::getCellIndex(unsigned int streetIndex, char startOrEnd)
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
