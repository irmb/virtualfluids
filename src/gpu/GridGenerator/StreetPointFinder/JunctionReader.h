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
//! \file JunctionReader.h
//! \ingroup StreetPointFinder
//! \author Stephan Lenz
//=======================================================================================
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
