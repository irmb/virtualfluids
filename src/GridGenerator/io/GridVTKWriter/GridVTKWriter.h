#ifndef GridVTKWriter_h
#define GridVTKWriter_h

#include <string>
#include <memory>

#include <GridGenerator/utilities/Transformator/TransformatorImp.h>


#include <VirtualFluidsDefinitions.h>

class Transformator;
struct Grid;

class VF_PUBLIC GridVTKWriter
{
public:
    static void writeGridToVTK(const Grid &grid, std::string name, std::shared_ptr<const Transformator> trans = std::make_shared<TransformatorImp>(), bool binaer = true);
    static void writeSparseGridToVTK(const Grid &grid, std::string name, std::shared_ptr<const Transformator> trans = std::make_shared<TransformatorImp>(), bool binaer = false);

private:
    GridVTKWriter();
    ~GridVTKWriter();

    static FILE *file;
    static bool binaer;

    static void initalVtkWriter(bool binaer, std::string name);
    static void writeVtkFile(std::shared_ptr<const Transformator> trans, const Grid & grid, const unsigned int &size);

    static void openFile(std::string name, std::string mode);
    static void closeFile();

    static void writeHeader();
    static void writePoints(std::shared_ptr<const Transformator> trans, const struct Grid &grid, const unsigned int &size);
    static void writeCells(const  unsigned  int &size);
    static void writeTypeHeader(const unsigned int &size);
    static void writeTypes(const Grid &grid, const unsigned int &size);

    static void end_line();
    static void force_big_endian(unsigned char *bytes);
    static void write_int(int val);
    static void write_float(float val);
};


#endif
