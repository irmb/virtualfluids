#ifndef SimulationFileWriter_H
#define SimulationFileWriter_H

#include <string>
#include <vector>
#include <fstream>
#include <memory>
#include <vector>

#include <GridGenerator/global.h>
#include <core/PointerDefinitions.h>
#include <core/DataTypes.h>
#include <core/NonCreatable.h>

class UnstructuredGridBuilder;
class GridBuilder;
class Grid;

enum class FILEFORMAT
{
    BINARY, ASCII
};

class SimulationFileWriter : private NonCreatable
{
public:
    VF_PUBLIC static void write(std::string folder, SPtr<GridBuilder> builder, FILEFORMAT format);

private:
    static void write(SPtr<GridBuilder> builder, FILEFORMAT format);
    static void openFiles();
    static void writeLevel(uint numberOfLevels);
    static void writeLevelSize(uint numberOfNodes, std::vector<std::vector<std::vector<real> > > qFiles);
    static void writeCoordFiles(SPtr<GridBuilder> builder, uint level, FILEFORMAT format);
    static void writeCoordsNeighborsGeo(SPtr<GridBuilder> builder, int index, uint level, FILEFORMAT format);
    static void writeLevelSizeGridInterface(uint sizeCF, uint sizeFC);
    static void writeGridInterfaceToFile(SPtr<GridBuilder> builder, uint level);
    static void writeGridInterfaceToFile(uint numberOfNodes, std::ofstream& coarseFile, uint* coarse, std::ofstream& fineFile, uint* fine, std::ofstream& offsetFile);
    static void writeBoundaryQsFile(std::vector<std::vector<std::vector<real> > > qFiles);
    static std::vector<std::vector<std::vector<real> > > createBCVectors(SPtr<Grid> grid);
    static void addShortQsToVector(int index, std::vector<std::vector<std::vector<real> > > &qs, SPtr<Grid> grid);
    static void addQsToVector(int index, std::vector<std::vector<std::vector<real> > > &qs, SPtr<Grid> grid);
    static void fillRBForNode(int index, int direction, int directionSign, int rb, std::vector<std::vector<std::vector<real> > > &qs, SPtr<Grid> grid);
    static void writeBoundary(std::vector<real> boundary, int rb);
	static void writeBoundaryShort(std::vector<real> boundary, int rb);
	static void closeFiles();


    static std::ofstream xCoordFile;
    static std::ofstream yCoordFile;
    static std::ofstream zCoordFile;
    static std::ofstream xNeighborFile;
    static std::ofstream yNeighborFile;
    static std::ofstream zNeighborFile;
    static std::ofstream geoVecFile;

    static std::ofstream scaleCF_coarse_File;
    static std::ofstream scaleCF_fine_File;
    static std::ofstream scaleFC_coarse_File;
    static std::ofstream scaleFC_fine_File;

    static std::ofstream offsetVecCF_File;
    static std::ofstream offsetVecFC_File;

    static std::vector<SPtr<std::ofstream> > qStreams;
    static std::vector<SPtr<std::ofstream> > valueStreams;

    static std::vector<std::ofstream> qFiles;
    static std::vector<std::ofstream> valueFiles;
    static std::string folder;
};


#endif
