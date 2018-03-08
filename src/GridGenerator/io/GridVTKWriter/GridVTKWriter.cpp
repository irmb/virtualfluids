#define _CRT_SECURE_NO_DEPRECATE
#include "GridVTKWriter.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>

#include <GridGenerator/utilities/Transformator/Transformator.h>
#include <utilities/logger/Logger.h>
#include <GridGenerator/grid/Grid.h>
#include <GridGenerator/geometries/Vertex/Vertex.cuh>
#include "basics/writer/WbWriterVtkXmlBinary.h"
#include "basics/container/CbArray3D.h"
#include "grid/NodeValues.h"
#include "grid/Cell.h"

FILE* GridVTKWriter::file = 0;
bool GridVTKWriter::binaer = true;

GridVTKWriter::GridVTKWriter(){}
GridVTKWriter::~GridVTKWriter(){}

void GridVTKWriter::writeSparseGridToVTK(SPtr<Grid> grid, std::string name, std::shared_ptr<const Transformator> trans, bool binaer)
{
    initalVtkWriter(binaer, name);
    writeVtkFile(trans, grid);
}

void GridVTKWriter::writeGridToVTKXML(SPtr<Grid> grid, const std::string name, bool binaer)
{
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleUInt8> cells;
    std::vector<std::string> nodedatanames;
    std::vector< std::vector<double> > nodedata;

    nodedatanames.push_back("types");


    nodedata.resize(nodedatanames.size());

    CbArray3D<int> nodeNumbers(grid->getNumberOfNodesX(), grid->getNumberOfNodesY(), grid->getNumberOfNodesZ(), -1);
    int nr = 0;

    for (real x = grid->getStartX(); x <= grid->getEndX(); x += grid->getDelta())
    {
        for (real y = grid->getStartY(); y <= grid->getEndY(); y += grid->getDelta())
        {
            for (real z = grid->getStartZ(); z <= grid->getEndZ(); z += grid->getDelta())
            {
                const auto xTranslate = int((x - grid->getStartX()) / grid->getDelta());
                const auto yTranslate = int((y - grid->getStartY()) / grid->getDelta());
                const auto zTranslate = int((z - grid->getStartZ()) / grid->getDelta());
                nodeNumbers(xTranslate, yTranslate, zTranslate) = nr++;
                nodes.push_back(UbTupleFloat3(float(x), float(y),float(z)));

                const char type = grid->getFieldEntry(grid->transCoordToIndex(x, y, z));
                nodedata[0].push_back(type);
            }
        }
    }

    int SWB, SEB, NEB, NWB, SWT, SET, NET, NWT;
    for (real x = grid->getStartX(); x < grid->getEndX(); x += grid->getDelta())
    {
        for (real y = grid->getStartY(); y < grid->getEndY(); y += grid->getDelta())
        {
            for (real z = grid->getStartZ(); z < grid->getEndZ(); z += grid->getDelta())
            {
                const auto xTranslate = int((x - grid->getStartX()) / grid->getDelta());
                const auto yTranslate = int((y - grid->getStartY()) / grid->getDelta());
                const auto zTranslate = int((z - grid->getStartZ()) / grid->getDelta());

                if ((SWB = nodeNumbers(xTranslate, yTranslate, zTranslate)) >= 0
                    && (SEB = nodeNumbers(xTranslate + 1, yTranslate, zTranslate)) >= 0
                    && (NEB = nodeNumbers(xTranslate + 1, yTranslate + 1, zTranslate)) >= 0
                    && (NWB = nodeNumbers(xTranslate, yTranslate + 1, zTranslate)) >= 0
                    && (SWT = nodeNumbers(xTranslate, yTranslate, zTranslate + 1)) >= 0
                    && (SET = nodeNumbers(xTranslate + 1, yTranslate, zTranslate + 1)) >= 0
                    && (NET = nodeNumbers(xTranslate + 1, yTranslate + 1, zTranslate + 1)) >= 0
                    && (NWT = nodeNumbers(xTranslate, yTranslate + 1, zTranslate + 1)) >= 0)
                {
                    Cell cell(x, y, z, grid->getDelta());
                    if(grid->nodeInCellIs(cell, OUT_OF_GRID) || grid->nodeInCellIs(cell, INVALID_NODE))
                        continue;

                    cells.push_back(makeUbTuple(uint(SWB), uint(SEB), uint(NEB), uint(NWB), uint(SWT), uint(SET), uint(NET), uint(NWT)));
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

void GridVTKWriter::writeVtkFile(std::shared_ptr<const Transformator> trans, SPtr<Grid> grid)
{
    GridVTKWriter::writeHeader();
    GridVTKWriter::writePoints(trans, grid);
    GridVTKWriter::writeCells(grid->getSize());
    GridVTKWriter::writeTypeHeader(grid->getSize());
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

void GridVTKWriter::writePoints(std::shared_ptr<const Transformator> trans, SPtr<Grid> grid)
{
    fprintf(file, "POINTS %d float\n", grid->getSize());
    real x, y, z;
    for (unsigned int i = 0; i < grid->getSize(); i++) {

        /*if (grid->getSparseIndex(i) == -1)
            continue;*/

        grid->transIndexToCoords(i, x, y, z);
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

void GridVTKWriter::writeTypes(SPtr<Grid> grid)
{
    for (unsigned int i = 0; i < grid->getSize(); i++) 
    {
        /*if (grid->getSparseIndex(i) == -1)
            continue;*/

        if (binaer)
            write_int(grid->getFieldEntry(i));
        else
            fprintf(file, "%d ", grid->getFieldEntry(i));
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
