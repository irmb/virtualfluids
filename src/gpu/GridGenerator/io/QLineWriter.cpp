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
#include "QLineWriter.h"

#include <vector>
#include <string>
#include <fstream>

#include "basics/utilities/UbTuple.h"

#include "geometries/Vertex/Vertex.h"

#include "grid/BoundaryConditions/BoundaryCondition.h"
#include "grid/Grid.h"

using namespace std;
void writeLines(std::string filename, std::vector<UbTupleFloat3> nodes, std::vector<UbTupleInt2> lines);

void QLineWriter::writeArrows(std::string fileName, SPtr<GeometryBoundaryCondition> geometryBoundaryCondition, SPtr<Grid> grid)
{
    if (geometryBoundaryCondition == nullptr)
    {
        VF_LOG_WARNING("(QLineWriter::writeArrows) no geometry bc on this grid level.");
        return;
    }
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleInt2> cells;

    int actualNodeNumber = 0;
    for (std::size_t index = 0; index < geometryBoundaryCondition->indices.size(); index++)
    {
        Vertex startNode = getVertex(geometryBoundaryCondition->indices[index], grid);
        for (int qi = 0; qi <= 26; qi++)
        {
            real qval = geometryBoundaryCondition->qs[index][qi];
            if (qval > 0.0f)
            {
                Vertex dir((real)grid->getDirection()[qi * DIMENSION + 0], (real)grid->getDirection()[qi * DIMENSION + 1], (real)grid->getDirection()[qi * DIMENSION + 2]);
                Vertex nodeOnGeometry(startNode + (dir * qval)*grid->getDelta());

                nodes.push_back(makeUbTuple(float(startNode.x), float(startNode.y), float(startNode.z)));
                nodes.push_back(makeUbTuple(float(nodeOnGeometry.x), float(nodeOnGeometry.y), float(nodeOnGeometry.z)));
                actualNodeNumber += 2;
                cells.push_back(makeUbTuple(actualNodeNumber - 2, actualNodeNumber - 1));
            }
        }
    }

    writeLines(fileName, nodes, cells);
}

Vertex QLineWriter::getVertex(int matrixIndex, SPtr<Grid> grid)
{
    real x, y, z;
    grid->transIndexToCoords(matrixIndex, x, y, z);
    return Vertex(x, y, z);
}


void writeLines(std::string filename, std::vector<UbTupleFloat3> nodes, std::vector<UbTupleInt2> lines)
{
    string vtkfilename = filename + ".bin.vtu";

    ofstream out(vtkfilename.c_str(), ios::out | ios::binary);

    int nofNodes = (int)nodes.size();
    int nofCells = (int)lines.size();

    int bytesPerByteVal = 4; //==sizeof(int)
    int bytesPoints = 3 /*x1/x2/x3        */ * nofNodes * sizeof(float);
    int bytesCellConnectivty = 2 /*nodes per line */ * nofCells * sizeof(int);
    int bytesCellOffsets = 1 /*offset per line */ * nofCells * sizeof(int);
    int bytesCellTypes = 1 /*type of line */ * nofCells * sizeof(unsigned char);

    int offset = 0;
    //VTK FILE
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >" << "\n";
    out << "   <UnstructuredGrid>" << "\n";
    out << "      <Piece NumberOfPoints=\"" << nofNodes << "\" NumberOfCells=\"" << nofCells << "\">\n";

    //POINTS SECTION
    out << "         <Points>\n";
    out << "            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\"  />\n";
    out << "         </Points>\n";
    offset += (bytesPerByteVal + bytesPoints);

    //CELLS SECTION
    out << "         <Cells>\n";
    out << "            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += (bytesPerByteVal + bytesCellConnectivty);
    out << "            <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += (bytesPerByteVal + bytesCellOffsets);
    out << "            <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"" << offset << "\" />\n ";
    offset += (bytesPerByteVal + bytesCellTypes);
    out << "         </Cells>\n";

    out << "      </Piece>\n";
    out << "   </UnstructuredGrid>\n";

    // AppendedData SECTION
    out << "   <AppendedData encoding=\"raw\">\n";
    out << "_";

    //POINTS SECTION
    out.write((char*)&bytesPoints, bytesPerByteVal);
    for (int n = 0; n < nofNodes; n++)
    {
        out.write((char*)&val<1>(nodes[n]), sizeof(float));
        out.write((char*)&val<2>(nodes[n]), sizeof(float));
        out.write((char*)&val<3>(nodes[n]), sizeof(float));
    }

    //CELLS SECTION
    //cellConnectivity
    out.write((char*)&bytesCellConnectivty, bytesPerByteVal);
    for (int c = 0; c < nofCells; c++)
    {
        out.write((char*)&val<1>(lines[c]), sizeof(int));
        out.write((char*)&val<2>(lines[c]), sizeof(int));

    }

    //cellOffsets
    out.write((char*)&bytesCellOffsets, bytesPerByteVal);
    int itmp;
    for (int c = 1; c <= nofCells; c++)
    {
        itmp = 2 * c;
        out.write((char*)&itmp, sizeof(int));
    }

    //cellTypes
    out.write((char*)&bytesCellTypes, bytesPerByteVal);
    unsigned char vtkCellType = 3;
    for (int c = 0; c < nofCells; c++)
    {
        out.write((char*)&vtkCellType, sizeof(unsigned char));
    }
    out << "\n</AppendedData>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
}

//! \}
