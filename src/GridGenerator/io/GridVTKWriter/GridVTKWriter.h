#ifndef GridVTKWriter_h
#define GridVTKWriter_h

#include <string>

#include <core/PointerDefinitions.h>
#include <VirtualFluidsDefinitions.h>

enum class WRITING_FORMAT { BINARY, ASCII };

class Grid;

class VF_PUBLIC GridVTKWriter
{
public:
    static void writeSparseGridToVTK(SPtr<Grid> grid, const std::string& name, WRITING_FORMAT format = WRITING_FORMAT::ASCII);
    static void writeGridToVTKXML(SPtr<Grid> grid, const std::string& name, WRITING_FORMAT format = WRITING_FORMAT::ASCII);
    static void writeInterpolationCellsToVTKXML(SPtr<Grid> grid, SPtr<Grid> gridCoarse, const std::string& name, WRITING_FORMAT format = WRITING_FORMAT::ASCII);

private:
    GridVTKWriter() {}
    ~GridVTKWriter() {}

    static FILE *file;
    static WRITING_FORMAT format;

    static void initalVtkWriter(WRITING_FORMAT format, const std::string& name);

    static bool isBinaryWritingFormat();

    static void writeVtkFile(SPtr<Grid> grid);

    static void openFile(const std::string& name, const std::string& mode);
    static void closeFile();

    static void writeHeader();
    static void writePoints(SPtr<Grid> grid);
    static void writeCells(const unsigned int &size);
    static void writeTypeHeader(const unsigned int &size);
    static void writeTypes(SPtr<Grid> grid);

    static void end_line();
    static void force_big_endian(unsigned char *bytes);
    static void write_int(int val);
    static void write_float(float val);
};


#endif
