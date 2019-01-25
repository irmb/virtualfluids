#include "StreetPointFinder.h"

#include <string>
#include <sstream>
#include <fstream>

#include "grid/Grid.h"
#include "grid/NodeValues.h"

Street::Street(real xStart, real yStart, real xEnd, real yEnd, real dx)
{
    real length = sqrt( (xEnd - xStart)*(xEnd - xStart) + (yEnd - yStart)*(yEnd - yStart) );

    this->xStart = xStart + dx * (xEnd - xStart) / length;
    this->yStart = yStart + dx * (yEnd - yStart) / length;
    
    this->xEnd   = xEnd   - dx * (xEnd - xStart) / length;
    this->yEnd   = yEnd   - dx * (yEnd - yStart) / length;

    this->numberOfCells = lround( length / dx ) + 1;
}

real Street::getCoordinateX(int cellIndex)
{
    return xStart + real(cellIndex) / real(numberOfCells) * (xEnd - xStart);
}

real Street::getCoordinateY(int cellIndex)
{
    return yStart + real(cellIndex) / real(numberOfCells) * (yEnd - yStart);
}

void Street::findIndicesLB( SPtr<Grid> grid )
{
    for( uint i = 0; i < numberOfCells; i++ )
    {
        real x = getCoordinateX(i);
        real y = getCoordinateY(i);
        
        uint matrixIndex = grid->transCoordToIndex(x, y, 10.0);

        real xLB, yLB, zLB;
        grid->transIndexToCoords(matrixIndex, xLB, yLB, zLB);

        while( grid->getFieldEntry(matrixIndex) != BC_SOLID )
        {
            zLB -= grid->getDelta();
            matrixIndex = grid->transCoordToIndex( xLB, yLB, zLB );
        }

        std::stringstream msg;

        
        msg << "( " << x   << ", " << y   <<                " )" << "  ==>  ";
        msg << "( " << xLB << ", " << yLB << ", " << zLB << " ) [" << (int)grid->getFieldEntry(matrixIndex) << "] \n";

        *logging::out << logging::Logger::INFO_LOW << msg.str();

        this->matrixIndicesLB.push_back( matrixIndex );
        this->sparseIndicesLB.push_back( grid->getSparseIndex(matrixIndex) );
    }
}

void StreetPointFinder::readStreets(std::string filename)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::readStreets( " << filename << " )" << "\n";

    uint numberOfStreets;

    std::ifstream file;

    file.open( filename.c_str() );

    file >> numberOfStreets;

    for( uint i = 0; i < numberOfStreets; i++ )
    {
        real xStart, yStart, xEnd, yEnd;

        real dx;

        file >> xStart >> yStart >> xEnd >> yEnd >> dx;

        streets.push_back( Street(xStart, yStart, xEnd, yEnd, dx) );
    }

    file.close();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

void StreetPointFinder::findIndicesLB( SPtr<Grid> grid )
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::findIndicesLB()\n";

    for( auto& street : streets ) street.findIndicesLB( grid );

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

void StreetPointFinder::writeVTK(std::string filename)
{
    uint numberOfCells = 0;
    uint numberOfNodes = 0;

    for( auto& street : streets )
    {
        numberOfCells += street.numberOfCells;
        numberOfNodes += street.numberOfCells + 1;
    }

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::writeVTK( " << filename << " )" << "\n";

    std::ofstream file;

    file.open(filename);

    file << "# vtk DataFile Version 3.0\n";
    file << "by MeshGenerator\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    file << "POINTS " << numberOfNodes << " float" << std::endl;

    for( auto& street : streets )
    {
        for( uint i = 0; i <= street.numberOfCells; i++ )
        {
            file << 0.5 * ( street.getCoordinateX(i-1) + street.getCoordinateX(i) ) << " " 
                 << 0.5 * ( street.getCoordinateY(i-1) + street.getCoordinateY(i) ) << " " << 0.0 << std::endl;
        }
    }

    //////////////////////////////////////////////////////////////////////////

    file << "CELLS " << numberOfCells << " " << 3 * numberOfCells << std::endl;


    uint nodeIndex = 0;
    for( auto& street : streets )
    {
        for( uint i = 0; i < street.numberOfCells; i++ )
        {
            file << "2 " << nodeIndex << " " << nodeIndex + 1 << std::endl;
            nodeIndex++;
        }
        nodeIndex++;
    }

    //////////////////////////////////////////////////////////////////////////

    file << "CELL_TYPES " << numberOfCells << std::endl;

    for ( uint i = 0; i < numberOfCells; i++ ){
        file << "3" << std::endl;
    }
    //////////////////////////////////////////////////////////////////////////

    file << "\nCELL_DATA " << numberOfCells << std::endl;

    file << "FIELD Label " << 2 << std::endl;

    //////////////////////////////////////////////////////////////////////////

    file << "StreetIndex 1 " << numberOfCells << " int" << std::endl;

    uint streetIndex = 0;
    for( auto& street : streets )
    {
        for( uint i = 0; i < street.numberOfCells; i++ )
        {
            file << streetIndex << std::endl;
        }
        streetIndex++;
    }

    //////////////////////////////////////////////////////////////////////////

    file << "length 1 " << numberOfCells << " float" << std::endl;

    for( auto& street : streets )
    {
        for( uint i = 0; i < street.numberOfCells; i++ )
        {
            real length = sqrt( ( street.getCoordinateX(i) - street.getCoordinateX(0) ) * ( street.getCoordinateX(i) - street.getCoordinateX(0) )
                              + ( street.getCoordinateY(i) - street.getCoordinateY(0) ) * ( street.getCoordinateY(i) - street.getCoordinateY(0) ) );

            file << length << std::endl;
        }
    }

    //////////////////////////////////////////////////////////////////////////

    file.close();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}


void StreetPointFinder::writeConnectionVTK(std::string filename, SPtr<Grid> grid)
{
    uint numberOfCells = 0;

    for( auto& street : streets )
    {
        numberOfCells += street.numberOfCells;
    }

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::writeConnectionVTK( " << filename << " )" << "\n";

    std::ofstream file;

    file.open(filename);

    file << "# vtk DataFile Version 3.0\n";
    file << "by MeshGenerator\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    file << "POINTS " << 2 * numberOfCells << " float" << std::endl;

    for( auto& street : streets )
    {
        for( uint i = 0; i < street.numberOfCells; i++ )
        {
            real xLB, yLB, zLB;
            grid->transIndexToCoords(street.matrixIndicesLB[i], xLB, yLB, zLB);

            file << street.getCoordinateX(i) << " " << street.getCoordinateY(i) << " " << 5.0 << std::endl;
            file << xLB << " " << yLB << " " << zLB << std::endl;
        }
    }

    //////////////////////////////////////////////////////////////////////////

    file << "CELLS " << numberOfCells << " " << 3 * numberOfCells << std::endl;


    uint nodeIndex = 0;
    for( auto& street : streets )
    {
        for( uint i = 0; i < street.numberOfCells; i++ )
        {
            file << "2 " << nodeIndex << " " << nodeIndex + 1 << std::endl;
            nodeIndex += 2;
        }
    }

    //////////////////////////////////////////////////////////////////////////

    file << "CELL_TYPES " << numberOfCells << std::endl;

    for ( uint i = 0; i < numberOfCells; i++ ){
        file << "3" << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////

    file.close();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

void StreetPointFinder::writeSimulationFile(std::string gridPath, real concentration, uint numberOfLevels, uint level)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "StreetPointFinder::writeSimulationFile( " << gridPath << "conc.dat )" << "\n";

    std::ofstream file;

    file.open(gridPath + "conc.dat");

    file << "concentration\n";

    file << numberOfLevels - 1 << "\n";

    for( uint currentLevel = 0; currentLevel < numberOfLevels; currentLevel++ )
    {
        if( currentLevel == level )
        {
            uint numberOfCells = 0;
            for( auto& street : streets )
            {
                numberOfCells += street.numberOfCells;
            }

            file << numberOfCells << "\n";

            for( auto& street : streets )
            {
                for( auto& sparseIndexLB : street.sparseIndicesLB )
                {
                    // + 1 for numbering shift between GridGenerator and VF_GPU
                    file << sparseIndexLB + 1 << "\n";
                }
            }
        }
        else
        {
            file << "0\n";
        }
    }

    file.close();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}
