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

FILE* GridVTKWriter::file = 0;
bool GridVTKWriter::binaer = true;

GridVTKWriter::GridVTKWriter(){}
GridVTKWriter::~GridVTKWriter(){}

void GridVTKWriter::writeGridToVTK(const Grid &grid, std::string name, std::shared_ptr<const Transformator> trans, bool binaer)
{
    initalVtkWriter(binaer, name);
    writeVtkFile(trans, grid, grid.size);
}

void GridVTKWriter::writeSparseGridToVTK(const Grid &grid, std::string name, std::shared_ptr<const Transformator> trans, bool binaer)
{
    initalVtkWriter(binaer, name);
    writeVtkFile(trans, grid, grid.reducedSize);
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

void GridVTKWriter::writeVtkFile(std::shared_ptr<const Transformator> trans, const Grid & grid, const unsigned int &size)
{
    GridVTKWriter::writeHeader();
    GridVTKWriter::writePoints(trans, grid, size);
    GridVTKWriter::writeCells(size);
    GridVTKWriter::writeTypeHeader(size);
    GridVTKWriter::writeTypes(grid, size);
    GridVTKWriter::closeFile();

    *logging::out << logging::Logger::INTERMEDIATE << "Output file closed\n";
}

void GridVTKWriter::openFile(std::string name, std::string mode)
{
    file = fopen(name.c_str(), mode.c_str());
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

void GridVTKWriter::writePoints(std::shared_ptr<const Transformator> trans, const struct Grid &grid, const unsigned int &size)
{
    fprintf(file, "POINTS %d float\n", size);

    unsigned int x, y, z;
    for (unsigned int i = 0; i < size; i++) {
        grid.transIndexToCoords(grid.matrixIndex[i], x, y, z);
        Vertex v((real)x, (real)y, (real)z);
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

void GridVTKWriter::writeTypes(const Grid &grid, const unsigned int &size)
{
    for (unsigned int i = 0; i < size; i++) {
        if (binaer)
            write_int(grid.field[grid.matrixIndex[i]]);
        else
            fprintf(file, "%d ", grid.field[grid.matrixIndex[i]]);
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
