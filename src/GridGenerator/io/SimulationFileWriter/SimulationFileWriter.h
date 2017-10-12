#ifndef SimulationFileWriter_H
#define SimulationFileWriter_H

#include "GridGenerator_EXPORT.h"

#include <string>
#include <vector>
#include <fstream>
#include <memory>

#include <GridGenerator/global.h>
#include <GridGenerator/grid/Grid.cuh>

class Transformator;
class UnstructuredGridBuilder;
class GridBuilder;

class SimulationFileWriter
{
public:
	GridGenerator_EXPORT static void writeSimulationFiles(std::string folder, std::shared_ptr<GridBuilder> builder, bool binaer, std::shared_ptr<Transformator> trans);

private:
    SimulationFileWriter() {};
    ~SimulationFileWriter() {};
    static void writeCoordFiles(bool binaer, std::shared_ptr<GridBuilder> builder, std::shared_ptr<Transformator> trans);
    static void writeBoundaryQsFile(std::vector<std::vector<std::vector<doubflo> > > qFiles);

    static void openFiles();
    static void writeLevelAndLevelSize(int sizeCoords, std::vector<std::vector<std::vector<doubflo> > > qFiles);
    static void writeCoordsNeighborsGeo(const int& index, bool binaer, std::shared_ptr<GridBuilder> builder, std::shared_ptr<Transformator> trans);
    static void writeBoundary(std::vector<doubflo> boundary, int rb);
    static void closeFiles();

	static std::ofstream xCoordFile;
	static std::ofstream yCoordFile;
	static std::ofstream zCoordFile;
	static std::ofstream xNeighborFile;
	static std::ofstream yNeighborFile;
	static std::ofstream zNeighborFile;
	static std::ofstream geoVecFile;

    static std::vector<std::shared_ptr<std::ofstream> > qStreams;
    static std::vector<std::shared_ptr<std::ofstream> > valueStreams;

    static std::vector<std::ofstream> qFiles;
    static std::vector<std::ofstream> valueFiles;
    static std::string folder;
};


#endif
