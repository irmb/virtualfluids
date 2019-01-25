#ifndef StreetPointFinder_H
#define StreetPointFinder_H

#include <vector>

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"

#include "Core/Logger/Logger.h"

#include <VirtualFluidsDefinitions.h>

class Grid;

struct VF_PUBLIC Street
{
    uint numberOfCells;
    real xStart, yStart, xEnd, yEnd;

    std::vector<uint> matrixIndicesLB;
    std::vector<uint> sparseIndicesLB;

    Street( real xStart, real yStart, real xEnd, real yEnd, real dx );

    real getCoordinateX( int cellIndex );
    real getCoordinateY( int cellIndex );

    void findIndicesLB( SPtr<Grid> grid );
};

struct VF_PUBLIC StreetPointFinder
{
    std::vector<Street> streets;

    void readStreets(std::string filename);

    void findIndicesLB( SPtr<Grid> grid );

    void writeVTK( std::string filename );

    void writeConnectionVTK(std::string filename, SPtr<Grid> grid);

    void writeSimulationFile( std::string gridPath, real concentration, uint numberOfLevels, uint level );
};


#endif