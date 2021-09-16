//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file SinkReader.cpp
//! \ingroup StreetPointFinder
//! \author Stephan Lenz
//=======================================================================================
#include "SinkReader.h"

#include <fstream>
#include <iostream>

SinkReaderData::SinkReaderData(uint sinkIndex, float sinkBlockedPossibility) :
	sinkIndex{ sinkIndex }, sinkBlockedPossibility{ sinkBlockedPossibility }
{}

void SinkReader::readSinks(std::string filename, StreetPointFinder* streetPointFinder)
{
	*logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::readSinks( " << filename << " )" << "\n";

	this->streetPointFinder = streetPointFinder;

	std::ifstream file;
	file.open(filename.c_str());
	if (!file.is_open()) std::cerr << "File not found" << std::endl;

	uint numberOfSinks;
	file >> numberOfSinks;

	uint streetIndex;
	float sinkBlockedPossibility;


	for (uint i = 0; i < numberOfSinks; i++) {
		file >> streetIndex >> sinkBlockedPossibility;
		sinks.push_back(SinkReaderData(getCellIndexEnd(streetIndex), sinkBlockedPossibility));
	}
}

unsigned int SinkReader::getCellIndexEnd(unsigned int streetIndex)
{
	uint i = 0;
	unsigned int cellIndex = 0;
	while (i < streetIndex) {
		cellIndex += streetPointFinder->streets[i].numberOfCells;
		++i;
	}
	
	return cellIndex + streetPointFinder->streets[streetIndex].numberOfCells - 1;
}
