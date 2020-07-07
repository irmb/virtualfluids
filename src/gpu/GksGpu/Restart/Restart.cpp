#include "Restart.h"

#include <iostream>
#include <fstream>

#include "Core/PointerDefinitions.h"
#include "Core/RealConstants.h"
#include "Core/Logger/Logger.h"

#include "DataBase/DataBase.h"

#include "Definitions/MemoryAccessPattern.h"

namespace GksGpu {

void Restart::writeRestart( SPtr<DataBase> dataBase, std::string filename, uint iter )
{
    filename += ".rst";

    *logging::out << logging::Logger::INFO_HIGH << "Writing restart file " << filename << " ... ";
	
    std::ofstream file;

	file.open( filename.c_str(), std::ios::binary );

	if (!file.is_open()) {
		throw std::runtime_error("\nFile cannot be opened.\n\nERROR!\n\n\n");
        return;
	}

    //////////////////////////////////////////////////////////////////////////

    file.write( (char*) &iter, sizeof( uint ) );

    file.write( (char*) &dataBase->numberOfLevels, sizeof( uint ) );
    file.write( (char*) &dataBase->numberOfCells,  sizeof( uint ) );
    file.write( (char*) &dataBase->numberOfFaces,  sizeof( uint ) );

    file.write( (char*) dataBase->dataHost.data(), LENGTH_CELL_DATA * dataBase->numberOfCells * sizeof( real ) );

    file.close();

    *logging::out << logging::Logger::INFO_HIGH << "done!\n";
}

bool Restart::readRestart( SPtr<DataBase> dataBase, std::string filename, uint& iter )
{
    filename += ".rst";

    *logging::out << logging::Logger::INFO_HIGH << "Reading restart file " << filename << " ... ";
	
    std::ifstream file;

	file.open( filename.c_str(), std::ios::binary );

	if (!file.is_open()) {
		throw std::runtime_error("\nFile cannot be opened.\n\nERROR!\n\n\n");
        return false;
	}

    //////////////////////////////////////////////////////////////////////////

    file.read( (char*) &iter, sizeof( uint ) );

    uint numberOfLevelsRead;
    uint numberOfCellsRead;
    uint numberOfFacesRead;
    
    file.read( (char*) &numberOfLevelsRead, sizeof( uint ) );
    file.read( (char*) &numberOfCellsRead,  sizeof( uint ) );
    file.read( (char*) &numberOfFacesRead,  sizeof( uint ) );

    if( numberOfLevelsRead != dataBase->numberOfLevels ||
        numberOfCellsRead  != dataBase->numberOfCells  ||
        numberOfFacesRead  != dataBase->numberOfFaces  ){
    
        *logging::out << logging::Logger::INFO_HIGH << "\n";
        *logging::out << logging::Logger::INFO_HIGH << "Levels: " << numberOfLevelsRead << " vs. " << dataBase->numberOfLevels << "\n";
        *logging::out << logging::Logger::INFO_HIGH << "Cells:  " << numberOfCellsRead  << " vs. " << dataBase->numberOfCells  << "\n";
        *logging::out << logging::Logger::INFO_HIGH << "Faces:  " << numberOfFacesRead  << " vs. " << dataBase->numberOfFaces  << "\n";

        file.close();

        throw std::runtime_error("\nERROR: Restart file does not match current setup");
    }

    file.read( (char*) dataBase->dataHost.data(), LENGTH_CELL_DATA * dataBase->numberOfCells * sizeof( real ) );

    file.close();

    *logging::out << logging::Logger::INFO_HIGH << "done!\n";

    return true;
}

} // namespace GksGpu