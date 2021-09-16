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
//! \file StreetPointFinder.h
//! \ingroup StreetPointFinder
//! \author Stephan Lenz
//=======================================================================================
#ifndef StreetPointFinder_H
#define StreetPointFinder_H

#include <vector>
#include <string>

#include "GridGenerator_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"



class Grid;

struct GRIDGENERATOR_EXPORT Street
{
    // The start and end coordinates are stored for cell centers!
    //
    //     |---x---|---x---|---x---|---x---|---x---|---x---|---x---|---x---|---x---|---x---|
    //         |--->                       |<----->|                               <---|
    //         xStart                          dx                                   xEnd
    //
    // dx = (xStart - xEnd) / (numberOfCells - 1)

    uint numberOfCells;
    real xStart, yStart, xEnd, yEnd;

    std::vector<uint> matrixIndicesLB;
    std::vector<uint> sparseIndicesLB;

    // The constructor expect start and end for cells
    Street( real xStartCell, real yStartCell, real xEndCell, real yEndCell, real dx );

    real getCoordinateX( int cellIndex );
    real getCoordinateY( int cellIndex );

	real getVectorX();
	real getVectorY();

    void findIndicesLB( SPtr<Grid> grid, real initialSearchHeight);
};

struct GRIDGENERATOR_EXPORT StreetPointFinder
{
    std::vector<Street> streets;

    std::vector<uint> sparseIndicesLB;
    std::vector<uint> mapNashToConc;

    void prepareSimulationFileData();

    void readStreets(std::string filename);

    void findIndicesLB( SPtr<Grid> grid, real initialSearchHeight );

	void writeVTK(std::string filename, const std::vector<int>& cars = std::vector<int>());

	void writeReducedVTK(std::string filename, const std::vector<int>& cars = std::vector<int>());

	void prepareWriteVTK(std::ofstream& file, uint & numberOfCells);

	void writeCarsVTK(std::ofstream& file, uint numberOfCells, const std::vector<int>& cars);

	void writeLengthsVTK(std::ofstream& file, uint numberOfCells);

	void writeStreetsVTK(std::ofstream& file, uint numberOfCells);

    void writeConnectionVTK(std::string filename, SPtr<Grid> grid);

	void writeSimulationFile(std::string gridPath, real concentration, uint numberOfLevels, uint level);

	void writeStreetVectorFile(std::string gridPath, real concentration, uint numberOfLevels, uint level);

    void writeSimulationFileSorted( std::string gridPath, real concentration, uint numberOfLevels, uint level );

    void writeMappingFile( std::string gridPath );

	//////////////////////////////////////////////////////////////////////////
	// 3D cars writer hacked by Stephan L.

	void write3DVTK(std::string filename, const std::vector<int>& cars = std::vector<int>());

	void prepareWrite3DVTK(std::ofstream& file, uint & numberOfCells, const std::vector<int>& cars);
};


#endif