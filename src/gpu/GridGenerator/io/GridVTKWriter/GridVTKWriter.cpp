//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_io io
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#endif
#include "GridVTKWriter.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>

#include "basics/writer/WbWriterVtkXmlBinary.h"
#include "basics/container/CbArray3D.h"

#include "geometries/Vertex/Vertex.h"

#include "grid/Grid.h"
#include "grid/NodeValues.h"
#include "grid/Cell.h"

using namespace vf::gpu;

FILE* GridVTKWriter::file = nullptr;
WRITING_FORMAT GridVTKWriter::format = WRITING_FORMAT::ASCII;


void GridVTKWriter::writeSparseGridToVTK(SPtr<Grid> grid, const std::string& name, WRITING_FORMAT format)
{
    initalVtkWriter(format, name);
    writeVtkFile(grid);
}

void GridVTKWriter::writeGridToVTKXML(SPtr<Grid> grid, const std::string& name)
{
    
    const uint chunkSize = 20000000; 

    const uint chunkSizeZ = chunkSize / ( grid->getNumberOfNodesX() * grid->getNumberOfNodesY() );

    for( uint startZ_Loop = 0, endZ_Loop = chunkSizeZ, part = 0; 
         endZ_Loop < grid->getNumberOfNodesZ() + chunkSizeZ; 
         startZ_Loop += chunkSizeZ, endZ_Loop += chunkSizeZ, part++ )

    {
        int endZ   = endZ_Loop;
        int startZ = startZ_Loop - 1;

        if(endZ >= (int)grid->getNumberOfNodesZ()) 
            endZ = grid->getNumberOfNodesZ();

        if( startZ < 0 )
            startZ = 0;

        std::vector<UbTupleFloat3> nodes;
        std::vector<UbTupleUInt8> cells;
        std::vector<std::string> nodedatanames;
        std::vector< std::vector<double> > nodedata;

        VF_LOG_INFO("Write Grid to XML VTK (*.vtu) output file : {}_Part_{}", name, part);

        nodedatanames.emplace_back("types");
        nodedatanames.emplace_back("sparse_id");
        nodedatanames.emplace_back("matrix_id");

        nodedata.resize(nodedatanames.size());

        CbArray3D<int> nodeNumbers(grid->getNumberOfNodesX(), grid->getNumberOfNodesY(), grid->getNumberOfNodesZ(), -1);
        int nr = 0;

        for (uint xIndex = 0; xIndex < grid->getNumberOfNodesX(); xIndex++)
        {
            for (uint yIndex = 0; yIndex < grid->getNumberOfNodesY(); yIndex++)
            {
                for (int zIndex = startZ; zIndex < endZ; zIndex++)
                {
                    real x, y, z;
                    uint index =
                        grid->getNumberOfNodesX() * grid->getNumberOfNodesY() * zIndex
                        + grid->getNumberOfNodesX() *                             yIndex
                        + xIndex;

                    grid->transIndexToCoords(index, x, y, z);

                    nodeNumbers(xIndex, yIndex, zIndex) = nr++;
                    nodes.emplace_back(UbTupleFloat3(float(x), float(y), float(z)));

                    const char type = grid->getFieldEntry(grid->transCoordToIndex(x, y, z));
                    nodedata[0].push_back(type);
                    nodedata[1].push_back(grid->getSparseIndex(index));
                    nodedata[2].push_back(index);
                }
            }
        }

        int SWB, SEB, NEB, NWB, SWT, SET, NET, NWT;
        for (uint xIndex = 0; xIndex < grid->getNumberOfNodesX() - 1; xIndex++)
        {
            for (uint yIndex = 0; yIndex < grid->getNumberOfNodesY() - 1; yIndex++)
            {
                for (int zIndex = startZ; zIndex < endZ - 1; zIndex++)
                {
                    real x, y, z;
                    uint index = grid->getNumberOfNodesX() * grid->getNumberOfNodesY() * zIndex
                        + grid->getNumberOfNodesX() *                             yIndex
                        + xIndex;

                    grid->transIndexToCoords(index, x, y, z);

                    if ((SWB = nodeNumbers(xIndex, yIndex, zIndex)) >= 0
                        && (SEB = nodeNumbers(xIndex + 1, yIndex, zIndex)) >= 0
                        && (NEB = nodeNumbers(xIndex + 1, yIndex + 1, zIndex)) >= 0
                        && (NWB = nodeNumbers(xIndex, yIndex + 1, zIndex)) >= 0
                        && (SWT = nodeNumbers(xIndex, yIndex, zIndex + 1)) >= 0
                        && (SET = nodeNumbers(xIndex + 1, yIndex, zIndex + 1)) >= 0
                        && (NET = nodeNumbers(xIndex + 1, yIndex + 1, zIndex + 1)) >= 0
                        && (NWT = nodeNumbers(xIndex, yIndex + 1, zIndex + 1)) >= 0)
                    {
                        Cell cell(x, y, z, grid->getDelta());
                        //if (grid->nodeInCellIs(cell, INVALID_OUT_OF_GRID) || grid->nodeInCellIs(cell, INVALID_COARSE_UNDER_FINE))
                        // continue;

                        cells.push_back(makeUbTuple(uint(SWB), uint(SEB), uint(NEB), uint(NWB), uint(SWT), uint(SET), uint(NET), uint(NWT)));
                    }
                }
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(name + "_Part_" + std::to_string(part), nodes, cells, nodedatanames, nodedata);
        VF_LOG_INFO("done.");
    }

}

void GridVTKWriter::writeInterpolationCellsToVTKXML(SPtr<Grid> grid, SPtr<Grid> gridCoarse, const std::string& name)
{
    std::vector<char> nodeInterpolationCellType( grid->getSize() );
    for( auto& type : nodeInterpolationCellType ) type = -1;

    std::vector<char> nodeOffset( grid->getSize() );
    for( auto& offset : nodeOffset ) offset = -1;

    std::vector<uint> matrixIndices( grid->getSparseSize() );

    for( uint matrixIndex = 0; matrixIndex < grid->getSize(); matrixIndex++ ){
        uint sparseIndex = grid->getSparseIndex(matrixIndex);
        if( sparseIndex != INVALID_INDEX )
            matrixIndices[ sparseIndex ] = matrixIndex;
    }

    for( uint index = 0; index < grid->getNumberOfNodesCF(); index++ ){
        nodeInterpolationCellType[ matrixIndices[ grid->getCF_coarse()[index] ] ] = grid->getFieldEntry( matrixIndices[ grid->getCF_coarse()[index] ] );

        nodeOffset               [ matrixIndices[ grid->getCF_coarse()[index] ] ] = grid->getCF_offset()[index];
    }

    //for( int index = 0; index < grid->getNumberOfNodesFC(); index++ ){
    //    nodeInterpolationCellType[ grid->getFC_coarse()[index] ] = grid->getFieldEntry( grid->getFC_coarse()[index] );
    //}

    if( gridCoarse ){

        for( uint index = 0; index < gridCoarse->getNumberOfNodesCF(); index++ ){
            nodeInterpolationCellType[ matrixIndices[ gridCoarse->getCF_fine()[index] ] ] = grid->getFieldEntry( matrixIndices[ gridCoarse->getCF_fine()[index] ] );
        }

        for( uint index = 0; index < gridCoarse->getNumberOfNodesFC(); index++ ){
            nodeInterpolationCellType[ matrixIndices[ gridCoarse->getFC_fine()[index] ] ] = grid->getFieldEntry( matrixIndices[ gridCoarse->getFC_fine()[index] ] );

            nodeOffset               [ matrixIndices[ gridCoarse->getFC_fine()[index] ] ] = gridCoarse->getFC_offset()[index];
        }
    }

    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleInt8> cells;
    std::vector<std::string> celldatanames;
    std::vector< std::vector<double> > celldata;

    celldatanames.emplace_back("InterpolationCells");
    celldatanames.emplace_back("Offset");

    celldata.resize(celldatanames.size());

    CbArray3D<int> nodeNumbers(grid->getNumberOfNodesX(), grid->getNumberOfNodesY(), grid->getNumberOfNodesZ(), -1);
    int nr = 0;

    for (uint xIndex = 0; xIndex < grid->getNumberOfNodesX(); xIndex++)
    {
        for (uint yIndex = 0; yIndex < grid->getNumberOfNodesY(); yIndex++)
        {
            for (uint zIndex = 0; zIndex < grid->getNumberOfNodesZ(); zIndex++)
            {
                real x, y, z;
                uint index = 
                      grid->getNumberOfNodesX() * grid->getNumberOfNodesY() * zIndex
                    + grid->getNumberOfNodesX() *                             yIndex
                    + xIndex;

                grid->transIndexToCoords(index, x, y, z);

                nodeNumbers(xIndex, yIndex, zIndex) = nr++;
                nodes.emplace_back(UbTupleFloat3(float(x), float(y), float(z)));
            }
        }
    }

    int SWB, SEB, NEB, NWB, SWT, SET, NET, NWT;
    for (uint xIndex = 0; xIndex < grid->getNumberOfNodesX() - 1; xIndex++)
    {
        for (uint yIndex = 0; yIndex < grid->getNumberOfNodesY() - 1; yIndex++)
        {
            for (uint zIndex = 0; zIndex < grid->getNumberOfNodesZ() - 1; zIndex++)
            {
                real x, y, z;
                uint index = grid->getNumberOfNodesX() * grid->getNumberOfNodesY() * zIndex
                    + grid->getNumberOfNodesX() *                             yIndex
                    + xIndex;

                grid->transIndexToCoords(index, x, y, z);

                if ((SWB = nodeNumbers(xIndex, yIndex, zIndex)) >= 0
                    && (SEB = nodeNumbers(xIndex + 1, yIndex, zIndex)) >= 0
                    && (NEB = nodeNumbers(xIndex + 1, yIndex + 1, zIndex)) >= 0
                    && (NWB = nodeNumbers(xIndex, yIndex + 1, zIndex)) >= 0
                    && (SWT = nodeNumbers(xIndex, yIndex, zIndex + 1)) >= 0
                    && (SET = nodeNumbers(xIndex + 1, yIndex, zIndex + 1)) >= 0
                    && (NET = nodeNumbers(xIndex + 1, yIndex + 1, zIndex + 1)) >= 0
                    && (NWT = nodeNumbers(xIndex, yIndex + 1, zIndex + 1)) >= 0)
                {
                    Cell cell(x, y, z, grid->getDelta());
                    //if (grid->nodeInCellIs(cell, INVALID_OUT_OF_GRID) || grid->nodeInCellIs(cell, INVALID_COARSE_UNDER_FINE))
                    //    continue;

                    cells.push_back(makeUbTuple(SWB, SEB, NEB, NWB, SWT, SET, NET, NWT));

                    //const char type = grid->getFieldEntry(grid->transCoordToIndex(nodes[SWB].v1, nodes[SWB].v2.v1, nodes[SWB].v2.v2));
                    //const char type = grid->getFieldEntry(grid->transCoordToIndex(val<1>(nodes[SWB]), val<2>(nodes[SWB]), val<3>(nodes[SWB])));
                    const char type = nodeInterpolationCellType[ grid->transCoordToIndex(val<1>(nodes[SWB]), val<2>(nodes[SWB]), val<3>(nodes[SWB])) ];
                    const char offset = nodeOffset               [ grid->transCoordToIndex(val<1>(nodes[SWB]), val<2>(nodes[SWB]), val<3>(nodes[SWB])) ];

                    celldata[0].push_back( type );
                    celldata[1].push_back( offset );
                }
            }
        }
    }
    WbWriterVtkXmlBinary::getInstance()->writeOctsWithCellData(name, nodes, cells, celldatanames, celldata);
}


/*#################################################################################*/
/*---------------------------------private methods---------------------------------*/
/*---------------------------------------------------------------------------------*/
void GridVTKWriter::initalVtkWriter(WRITING_FORMAT format, const std::string& name)
{
    GridVTKWriter::format = format;

    VF_LOG_INFO("Write Grid to vtk output file: {}", name);

    std::string mode = "w";
    if (isBinaryWritingFormat())
        mode = "wb";
    GridVTKWriter::openFile(name, mode);

    VF_LOG_INFO("Output file opened ...");
}

bool GridVTKWriter::isBinaryWritingFormat()
{
    return GridVTKWriter::format == WRITING_FORMAT::BINARY;
}

void GridVTKWriter::writeVtkFile(SPtr<Grid> grid)
{
    GridVTKWriter::writeHeader();
    GridVTKWriter::writePoints(grid);
    GridVTKWriter::writeCells(grid->getSize());
    GridVTKWriter::writeTypeHeader(grid->getSize());
    GridVTKWriter::writeTypes(grid);
    GridVTKWriter::closeFile();

    VF_LOG_INFO("Output file closed");
}

void GridVTKWriter::openFile(const std::string& name, const std::string& mode)
{
    file = fopen(name.c_str(), mode.c_str());
    if(file==NULL)
        VF_LOG_CRITICAL("cannot open file {}", name);
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
    if (isBinaryWritingFormat())
        fprintf(file, "BINARY\n");
    else
        fprintf(file, "ASCII\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n");
}

void GridVTKWriter::writePoints(SPtr<Grid> grid)
{
    fprintf(file, "POINTS %d float\n", grid->getSize());
    real x, y, z;
    for (unsigned int i = 0; i < grid->getSize(); i++) {

        /*if (grid->getSparseIndex(i) == -1)
            continue;*/

        grid->transIndexToCoords(i, x, y, z);

        if (isBinaryWritingFormat()) {
            write_float(float(x));
            write_float(float(y));
            write_float(float(z));
        }
        else
            fprintf(file, "%f %f %f\n", x, y, z);
    }
}

void GridVTKWriter::writeCells(const unsigned int &size)
{
    fprintf(file, "\nCELLS %u %u\n", size, size * 2);
    for (unsigned int i = 0; i < size; ++i)
    {
        if (isBinaryWritingFormat()){
            write_int(1);
            write_int(i);
        }
        else
            fprintf(file, "1 %u\n", i);
    }

    fprintf(file, "\nCELL_TYPES %u\n", size);
    for (unsigned int i = 0; i < size; ++i)
    {
        if (isBinaryWritingFormat())
            write_int(1);
        else
            fprintf(file, "1 ");
    }
    if (!isBinaryWritingFormat())
        GridVTKWriter::end_line();
}

void GridVTKWriter::writeTypeHeader(const unsigned int &size)
{
    fprintf(file, "\nPOINT_DATA %u\n", size);
    fprintf(file, "SCALARS type int\n");
    fprintf(file, "LOOKUP_TABLE default\n");
}

void GridVTKWriter::writeTypes(SPtr<Grid> grid)
{
    for (unsigned int i = 0; i < grid->getSize(); i++) 
    {
        /*if (grid->getSparseIndex(i) == -1)
            continue;*/

        if (isBinaryWritingFormat())
            write_int(grid->getFieldEntry(i));
        else
            fprintf(file, "%d ", grid->getFieldEntry(i));
    }
}

void GridVTKWriter::end_line()
{
    char str2[8] = "\n";
    fprintf(file, "%s", str2);
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

//! \}
