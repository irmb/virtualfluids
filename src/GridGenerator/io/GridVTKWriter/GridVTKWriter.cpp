#define _CRT_SECURE_NO_DEPRECATE
#include "GridVTKWriter.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>

#include <GridGenerator/utilities/Transformator/Transformator.h>
#include <utilities/logger/Logger.h>
#include <GridGenerator/grid/Grid.cuh>
#include <GridGenerator/geometries/Vertex/Vertex.cuh>
#include "basics/writer/WbWriterVtkXmlBinary.h"
#include "basics/container/CbArray3D.h"


FILE* GridVTKWriter::file = 0;
bool GridVTKWriter::binaer = true;

GridVTKWriter::GridVTKWriter(){}
GridVTKWriter::~GridVTKWriter(){}

void GridVTKWriter::writeSparseGridToVTK(const Grid &grid, std::string name, std::shared_ptr<const Transformator> trans, bool binaer)
{
    initalVtkWriter(binaer, name);
    writeVtkFile(trans, grid);
}

void GridVTKWriter::writeGridToVTKXML(const Grid& grid, const std::string name, bool binaer)
{
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleUInt8> cells;
    std::vector<std::string> nodedatanames;
    std::vector< std::vector<double> > nodedata;

    nodedatanames.push_back("cfc");
    nodedatanames.push_back("cff");
    nodedatanames.push_back("stopper");

    nodedata.resize(nodedatanames.size());

    CbArray3D<int> nodeNumbers(grid.nx, grid.ny, grid.nz, -1);
    int nr = 0;

    for (real x = grid.startX; x <= grid.endX; x += grid.delta)
    {
        for (real y = grid.startY; y <= grid.endY; y += grid.delta)
        {
            for (real z = grid.startZ; z <= grid.endZ; z += grid.delta)
            {
                const int xTranslate = (x - grid.startX) / grid.delta;
                const int yTranslate = (y - grid.startY) / grid.delta;
                const int zTranslate = (z - grid.startZ) / grid.delta;
                nodeNumbers(xTranslate, yTranslate, zTranslate) = nr++;
                nodes.push_back(UbTupleFloat3(float(x), float(y),float(z)));

                nodedata[0].push_back(grid.field[grid.transCoordToIndex(x, y, z)]);
            }
        }
    }

    uint SWB, SEB, NEB, NWB, SWT, SET, NET, NWT;
    for (real x = grid.startX; x <= grid.endX - grid.delta; x += grid.delta)
    {
        for (real y = grid.startY; y <= grid.endY - grid.delta; y += grid.delta)
        {
            for (real z = grid.startZ; z <= grid.endZ - grid.delta; z += grid.delta)
            {
                const int xTranslate = (x - grid.startX) / grid.delta;
                const int yTranslate = (y - grid.startY) / grid.delta;
                const int zTranslate = (z - grid.startZ) / grid.delta;

                if ((SWB = nodeNumbers(xTranslate, yTranslate, zTranslate)) >= 0
                    && (SEB = nodeNumbers(xTranslate + 1, yTranslate, zTranslate)) >= 0
                    && (NEB = nodeNumbers(xTranslate + 1, yTranslate + 1, zTranslate)) >= 0
                    && (NWB = nodeNumbers(xTranslate, yTranslate + 1, zTranslate)) >= 0
                    && (SWT = nodeNumbers(xTranslate, yTranslate, zTranslate + 1)) >= 0
                    && (SET = nodeNumbers(xTranslate + 1, yTranslate, zTranslate + 1)) >= 0
                    && (NET = nodeNumbers(xTranslate + 1, yTranslate + 1, zTranslate + 1)) >= 0
                    && (NWT = nodeNumbers(xTranslate, yTranslate + 1, zTranslate + 1)) >= 0)
                {
                    cells.push_back(makeUbTuple(SWB, SEB, NEB, NWB, SWT, SET, NET, NWT));
                }
            }
        }
    }


    WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(name, nodes, cells, nodedatanames, nodedata);

}


/*#################################################################################*/
/*---------------------------------private methods---------------------------------*/
/*---------------------------------------------------------------------------------*/
void GridVTKWriter::initalVtkWriter(bool binaer, std::string name)
{
    GridVTKWriter::binaer = binaer;
    name += ".vtk";

    *logging::out << logging::Logger::INTERMEDIATE << "Write Grid to vtk output file : " + name + "\n";

    std::string mode = "w";
    if (binaer)
        mode = "wb";
    GridVTKWriter::openFile(name, mode);

    *logging::out << logging::Logger::INTERMEDIATE << "  Output file opened ...\n";
}

void GridVTKWriter::writeVtkFile(std::shared_ptr<const Transformator> trans, const Grid & grid)
{
    GridVTKWriter::writeHeader();
    GridVTKWriter::writePoints(trans, grid);
    GridVTKWriter::writeCells(grid.reducedSize);
    GridVTKWriter::writeTypeHeader(grid.reducedSize);
    GridVTKWriter::writeTypes(grid);
    GridVTKWriter::closeFile();

    *logging::out << logging::Logger::INTERMEDIATE << "Output file closed\n";
}

void GridVTKWriter::openFile(std::string name, std::string mode)
{
    file = fopen(name.c_str(), mode.c_str());
    if(file==NULL)
        *logging::out << logging::Logger::HIGH << "  cannot open file ...\n";
}

void GridVTKWriter::closeFile()
{
    GridVTKWriter::end_line();
	fclose(file);
}

void GridVTKWriter::writeHeader()
{
	fprintf(file, "# vtk DataFile Version 3.0\n");
	fprintf(file, "by MeshGenerator\n");
	if (binaer)
		fprintf(file, "BINARY\n");
	else
		fprintf(file, "ASCII\n");
	fprintf(file, "DATASET UNSTRUCTURED_GRID\n");
}

void GridVTKWriter::writePoints(std::shared_ptr<const Transformator> trans, const struct Grid &grid)
{
    fprintf(file, "POINTS %d float\n", grid.reducedSize);
    real x, y, z;
    for (unsigned int i = 0; i < grid.size; i++) {

        if (grid.matrixIndex[i] == -1)
            continue;

        grid.transIndexToCoords(i, x, y, z);
        Vertex v(x,y,z);
        trans->transformGridToWorld(v);
        if (binaer) {
            write_float((float)v.x);
            write_float((float)v.y);
            write_float((float)v.z);
        }
        else
            fprintf(file, "%f %f %f\n", v.x, v.y, v.z);
    }
}

void GridVTKWriter::writeCells(const unsigned int &size)
{
	fprintf(file, "\nCELLS %d %d\n", size, size * 2);
	for (unsigned int i = 0; i < size; ++i)
	{
		if (binaer){
			write_int(1);
			write_int(i);
		}
		else
			fprintf(file, "1 %d\n", i);
	}

	fprintf(file, "\nCELL_TYPES %d\n", size);
	for (unsigned int i = 0; i < size; ++i)
	{
		if (binaer)
			write_int(1);
		else
			fprintf(file, "1 ");
	}
	if (!binaer)
        GridVTKWriter::end_line();
}

void GridVTKWriter::writeTypeHeader(const unsigned int &size)
{
	fprintf(file, "\nPOINT_DATA %d\n", size);
	fprintf(file, "SCALARS type int\n");
	fprintf(file, "LOOKUP_TABLE default\n");
}

void GridVTKWriter::writeTypes(const Grid &grid)
{
    for (unsigned int i = 0; i < grid.size; i++) 
    {
        if (grid.matrixIndex[i] == -1)
            continue;

        if (binaer)
            write_int(grid.field[i]);
        else
            fprintf(file, "%d ", grid.field[i]);
    }
}

void GridVTKWriter::end_line()
{
	char str2[8] = "\n";
	fprintf(file, str2);
}

void GridVTKWriter::write_int(int val)
{
	force_big_endian((unsigned char *)&val);
	fwrite(&val, sizeof(int), 1, file);
}

void GridVTKWriter::write_float(float val)
{
	force_big_endian((unsigned char *)&val);
	fwrite(&val, sizeof(float), 1, file);
}


void GridVTKWriter::force_big_endian(unsigned char *bytes)
{
	bool shouldSwap = false;
	int tmp1 = 1;
	unsigned char *tmp2 = (unsigned char *)&tmp1;
	if (*tmp2 != 0)
		shouldSwap = true;

	if (shouldSwap)
	{
		unsigned char tmp = bytes[0];
		bytes[0] = bytes[3];
		bytes[3] = tmp;
		tmp = bytes[1];
		bytes[1] = bytes[2];
		bytes[2] = tmp;
	}
}
