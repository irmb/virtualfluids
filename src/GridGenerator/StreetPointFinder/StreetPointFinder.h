#ifndef StreetPointFinder_H
#define StreetPointFinder_H

#include <vector>

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"

#include <VirtualFluidsDefinitions.h>

class Grid;

struct VF_PUBLIC Street
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

    void findIndicesLB( SPtr<Grid> grid );
};

struct VF_PUBLIC StreetPointFinder
{
    std::vector<Street> streets;

    std::vector<uint> sparseIndicesLB;
    std::vector<uint> mapNashToConc;

    void prepareSimulationFileData();

    void readStreets(std::string filename);

    void findIndicesLB( SPtr<Grid> grid );

	void writeVTK(std::string filename, const std::vector<int>& cars = std::vector<int>());

	void writeReducedVTK(std::string filename, const std::vector<int>& cars = std::vector<int>());

	void prepareWriteVTK(std::ofstream& file, uint & numberOfCells);

	void writeCarsVTK(std::ofstream& file, uint numberOfCells, const std::vector<int>& cars);

	void writeLengthsVTK(std::ofstream& file, uint numberOfCells);

	void writeStreetsVTK(std::ofstream& file, uint numberOfCells);

    void writeConnectionVTK(std::string filename, SPtr<Grid> grid);

    void writeSimulationFile( std::string gridPath, real concentration, uint numberOfLevels, uint level );

    void writeSimulationFileSorted( std::string gridPath, real concentration, uint numberOfLevels, uint level );

    void writeMappingFile( std::string gridPath );
};


#endif