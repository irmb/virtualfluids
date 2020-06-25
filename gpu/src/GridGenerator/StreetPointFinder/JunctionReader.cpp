#include "JunctionReader.h"

#include <fstream>
#include <iostream>
#include <string>


JunctionReaderData::JunctionReaderData(std::vector<uint> inCells, std::vector<uint> outCells, std::vector<int> carCanNotEnterThisOutCell, uint trafficLightSwitchTime = 0) :
	inCells{ inCells }, outCells{ outCells }, carCanNotEnterThisOutCell{ carCanNotEnterThisOutCell }, trafficLightSwitchTime{ trafficLightSwitchTime }
{}

void JunctionReader::readJunctions(std::string filename, StreetPointFinder* streetPointFinder)
{
	*logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::readJunctions( " << filename << " )" << "\n";

	std::ifstream file;
	file.open(filename.c_str());
	if (!file.is_open()) std::cerr << "File not found" << std::endl;
	this->streetPointFinder = streetPointFinder;

	uint numberOfJunctions;
	file >> numberOfJunctions;

	std::string inOutDummy;
	int streetIndex = 0;
	uint trafficLightTime = 0;
	bool onlyNeighbors = false;

	file >> inOutDummy;

	for (uint i = 0; i < numberOfJunctions; i++) {
		std::vector<uint> inCells, outCells;
		std::vector<int> carCanNotEnterThisOutCell;

		//inCells
		file >> inOutDummy;
		while (inOutDummy.compare("out") != 0) {
			streetIndex = std::stoi(inOutDummy);

			if (streetIndex >= 0)
				inCells.push_back(getCellIndex(streetIndex, 'e'));

			file >> inOutDummy;
		}

		//outCells
		file >> inOutDummy;
		while (inOutDummy.compare("in") != 0 && inOutDummy.compare("end") != 0 && inOutDummy.compare("t") != 0 && inOutDummy.compare("c") != 0) {
			streetIndex = std::stoi(inOutDummy);

			if (streetIndex >= 0) {
				outCells.push_back(getCellIndex(streetIndex, 's'));
				if (carCanNotEnterThisOutCell.size() < inCells.size())
					carCanNotEnterThisOutCell.push_back(getCellIndex(streetIndex, 's'));
			}
			else if (streetIndex == -2) //no prohibited outCell
				carCanNotEnterThisOutCell.push_back(-2);

			file >> inOutDummy;
		}

		//trafficLightTime
		if (inOutDummy.compare("t") == 0) {
			file >> inOutDummy;
			trafficLightTime = std::stoi(inOutDummy);
			file >> inOutDummy;
		}
		else
			trafficLightTime = 0;

		// only neighbors (used for curves)
		if (inOutDummy.compare("c") == 0) {
			onlyNeighbors = true;
			file >> inOutDummy;
		}


		//make Junction or neighbors
		if (onlyNeighbors) {
			if (inCells.size() == 2 && outCells.size() == 2) {
				specialNeighbors.cells.insert(specialNeighbors.cells.end(), inCells.begin(), inCells.end());
				specialNeighbors.neighbors.push_back(outCells[1]);     
				specialNeighbors.neighbors.push_back(outCells[0]);

				onlyNeighbors = false;
			}
			else std::cerr << "can't add curve" << std::endl; continue;
		}
		else
			junctions.push_back(JunctionReaderData(inCells, outCells, carCanNotEnterThisOutCell, trafficLightTime));

	}
}


unsigned int JunctionReader::getCellIndex(unsigned int streetIndex, char startOrEnd)
{
	uint i = 0;
	unsigned int cellIndex = 0;
	while (i < streetIndex) {
		cellIndex += streetPointFinder->streets[i].numberOfCells;
		++i;
	}
	if (startOrEnd == 's') 	return cellIndex;
	return cellIndex + streetPointFinder->streets[streetIndex].numberOfCells - 1;
}

