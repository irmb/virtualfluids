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
//! \addtogroup writer
//! \ingroup basics
//! \{
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#include <basics/utilities/UbLogger.h>
#include <basics/writer/WbWriterVtkASCII.h>
#include <basics/writer/WbWriterVtkBinary.h>
#include <cstring>

using namespace std;

std::string WbWriterVtkBinary::writeLines(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                          std::vector<UbTupleInt2> &lines)
{
    return WbWriterVtkASCII::getInstance()->writeLines(filename, nodes, lines);
}
/*===============================================================================*/
std::string WbWriterVtkBinary::writeTriangles(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                              std::vector<UbTupleInt3> &cells)
{
    return WbWriterVtkASCII::getInstance()->writeTriangles(filename, nodes, cells);
}
/*===============================================================================*/
std::string WbWriterVtkBinary::writeTrianglesWithNodeData(const std::string &filename,
                                                          std::vector<UbTupleFloat3> &nodes,
                                                          std::vector<UbTupleInt3> &cells,
                                                          std::vector<std::string> &datanames,
                                                          std::vector<std::vector<double>> &nodedata)
{
    return WbWriterVtkASCII::getInstance()->writeTrianglesWithNodeData(filename, nodes, cells, datanames, nodedata);
}
/*===============================================================================*/
std::string WbWriterVtkBinary::writeQuads(const string &filename, vector<UbTupleFloat3> &nodes,
                                          vector<UbTupleInt4> &cells)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkBinary::writeQuads to " << vtkfilename << " - start");

    ofstream out(vtkfilename.c_str(), ofstream::out | ofstream::binary);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str(), ios::out | ios::binary);
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    // HEADER-SECTION
    // WRITE BIGENDIAN VtkBinary FILE
    bool swapByte = UbSystem::isLittleEndian();
    int nofNodes  = (int)nodes.size();
    int nofCells  = (int)cells.size();

    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "D3Q19MasterNodeGrid"
        << "\n";
    out << ""
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++) {
        float x1 = (float)val<1>(nodes[n]);
        float x2 = (float)val<2>(nodes[n]);
        float x3 = (float)val<3>(nodes[n]);

        if (swapByte) {
            UbSystem::swapByteOrder((unsigned char *)&x1, sizeof(float));
            UbSystem::swapByteOrder((unsigned char *)&x2, sizeof(float));
            UbSystem::swapByteOrder((unsigned char *)&x3, sizeof(float));
        }

        out.write((char *)&x1, sizeof(float));
        out.write((char *)&x2, sizeof(float));
        out.write((char *)&x3, sizeof(float));
    }
    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofCells << " " << nofCells * 5 << "\n";

    int nodesPerCellDummy = 4; // nofNodesPerCell
    if (swapByte)
        UbSystem::swapByteOrder((unsigned char *)&nodesPerCellDummy, sizeof(int));
    for (int c = 0; c < (int)cells.size(); c++) {
        int SW = val<1>(cells[c]);
        int SE = val<2>(cells[c]);
        int NE = val<3>(cells[c]);
        int NW = val<4>(cells[c]);
        if (swapByte) {
            UbSystem::swapByteOrder((unsigned char *)&SW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&SE, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&NW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&NE, sizeof(int));
        }

        out.write((char *)&nodesPerCellDummy, sizeof(int));
        out.write((char *)&SW, sizeof(int));
        out.write((char *)&SE, sizeof(int));
        out.write((char *)&NW, sizeof(int));
        out.write((char *)&NE, sizeof(int));
    }
    out << "\n";

    out << "CELL_TYPES " << (int)cells.size() << "\n";
    int celltype = 8;
    if (swapByte)
        UbSystem::swapByteOrder((unsigned char *)&celltype, sizeof(int));
    for (int c = 0; c < nofCells; c++)
        out.write((char *)&celltype, sizeof(int));

    out << endl;
    out.close();

    UBLOG(logDEBUG1, "WbWriterVtkBinary::writeQuads to " << vtkfilename << " - end");
    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkBinary::writeQuadsWithNodeData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                      vector<UbTupleInt4> &cells, vector<string> &datanames,
                                                      vector<vector<double>> &nodedata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkBinary::writeQuadsWithNodeData to " << vtkfilename << " - start");

    ofstream out(vtkfilename.c_str(), ofstream::out | ofstream::binary);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str(), ios::out | ios::binary);
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    // WRITE BIGENDIAN VtkBinary FILE
    bool swapByte = UbSystem::isLittleEndian();
    int nofNodes  = (int)nodes.size();
    int nofCells  = (int)cells.size();

    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "D3Q19MasterNodeGrid"
        << "\n";
    out << ""
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++) {
        float x1 = (float)val<1>(nodes[n]);
        float x2 = (float)val<2>(nodes[n]);
        float x3 = (float)val<3>(nodes[n]);

        if (swapByte) {
            UbSystem::swapByteOrder((unsigned char *)&x1, sizeof(float));
            UbSystem::swapByteOrder((unsigned char *)&x2, sizeof(float));
            UbSystem::swapByteOrder((unsigned char *)&x3, sizeof(float));
        }

        out.write((char *)&x1, sizeof(float));
        out.write((char *)&x2, sizeof(float));
        out.write((char *)&x3, sizeof(float));
    }
    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofCells << " " << nofCells * 5 << "\n";

    int nodesPerCellDummy = 4; // nofNodesPerCell
    if (swapByte)
        UbSystem::swapByteOrder((unsigned char *)&nodesPerCellDummy, sizeof(int));
    for (int c = 0; c < (int)cells.size(); c++) {
        int SW = val<1>(cells[c]);
        int SE = val<2>(cells[c]);
        int NE = val<3>(cells[c]);
        int NW = val<4>(cells[c]);
        if (swapByte) {
            UbSystem::swapByteOrder((unsigned char *)&SW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&SE, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&NW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&NE, sizeof(int));
        }

        out.write((char *)&nodesPerCellDummy, sizeof(int));
        out.write((char *)&SW, sizeof(int));
        out.write((char *)&SE, sizeof(int));
        out.write((char *)&NW, sizeof(int));
        out.write((char *)&NE, sizeof(int));
    }
    out << "\n";

    out << "CELL_TYPES " << (int)cells.size() << "\n";
    int celltype = 8;
    if (swapByte)
        UbSystem::swapByteOrder((unsigned char *)&celltype, sizeof(int));
    for (int c = 0; c < nofCells; c++)
        out.write((char *)&celltype, sizeof(int));

    out << endl;

    // DATA SECTION
    // write data section
    out << "POINT_DATA " << nofNodes << "\n";
    for (int s = 0; s < (int)datanames.size(); ++s) {
        if ((int)nodedata[s].size() != nofNodes)
            throw UbException(UB_EXARGS, "datasetsize must be equal to nofNodes");
        out << "SCALARS " << datanames[s] << " float 1 \n LOOKUP_TABLE default \n";
        for (int d = 0; d < (int)nodedata[s].size(); d++) {
            float dummy = (float)nodedata[s][d];
            if (swapByte)
                UbSystem::swapByteOrder((unsigned char *)&dummy, sizeof(float));
            out.write((const char *)&dummy, sizeof(float));
        }
        out << endl;
    }
    out.close();

    UBLOG(logDEBUG1, "WbWriterVtkBinary::writeQuadsWithNodeData to " << vtkfilename << " - end");
    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkBinary::writeQuadsWithCellData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                      vector<UbTupleInt4> &cells, vector<string> &datanames,
                                                      vector<vector<double>> &celldata)
{
    // HEADER-SECTION
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkBinary::writeQuadsWithCellData to " << vtkfilename << " - start");

    ofstream out(vtkfilename.c_str(), ofstream::out | ofstream::binary);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str(), ios::out | ios::binary);
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    // WRITE BIGENDIAN VtkBinary FILE
    bool swapByte = UbSystem::isLittleEndian();
    int nofNodes  = (int)nodes.size();
    int nofCells  = (int)cells.size();

    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "D3Q19MasterNodeGrid"
        << "\n";
    out << ""
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++) {
        float x1 = (float)val<1>(nodes[n]);
        float x2 = (float)val<2>(nodes[n]);
        float x3 = (float)val<3>(nodes[n]);

        if (swapByte) {
            UbSystem::swapByteOrder((unsigned char *)&x1, sizeof(float));
            UbSystem::swapByteOrder((unsigned char *)&x2, sizeof(float));
            UbSystem::swapByteOrder((unsigned char *)&x3, sizeof(float));
        }

        out.write((char *)&x1, sizeof(float));
        out.write((char *)&x2, sizeof(float));
        out.write((char *)&x3, sizeof(float));
    }
    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofCells << " " << nofCells * 5 << "\n";

    int nodesPerCellDummy = 4; // nofNodesPerCell
    if (swapByte)
        UbSystem::swapByteOrder((unsigned char *)&nodesPerCellDummy, sizeof(int));
    for (int c = 0; c < (int)cells.size(); c++) {
        int SW = val<1>(cells[c]);
        int SE = val<2>(cells[c]);
        int NE = val<3>(cells[c]);
        int NW = val<4>(cells[c]);
        if (swapByte) {
            UbSystem::swapByteOrder((unsigned char *)&SW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&SE, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&NW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&NE, sizeof(int));
        }

        out.write((char *)&nodesPerCellDummy, sizeof(int));
        out.write((char *)&SW, sizeof(int));
        out.write((char *)&SE, sizeof(int));
        out.write((char *)&NW, sizeof(int));
        out.write((char *)&NE, sizeof(int));
    }
    out << "\n";

    out << "CELL_TYPES " << (int)cells.size() << "\n";
    int celltype = 8;
    if (swapByte)
        UbSystem::swapByteOrder((unsigned char *)&celltype, sizeof(int));
    for (int c = 0; c < nofCells; c++)
        out.write((char *)&celltype, sizeof(int));

    out << endl;

    // DATA SECTION
    // write data section
    out << "CELL_DATA " << nofCells << "\n";
    for (int s = 0; s < (int)datanames.size(); ++s) {
        if ((int)celldata[s].size() != nofCells)
            throw UbException(UB_EXARGS, "datasetsize must be equal to nofNodes");
        out << "SCALARS " << datanames[s] << " float 1 \n LOOKUP_TABLE default \n";
        for (int d = 0; d < (int)celldata[s].size(); d++) {
            float dummy = (float)celldata[s][d];
            if (swapByte)
                UbSystem::swapByteOrder((unsigned char *)&dummy, sizeof(float));
            out.write((const char *)&dummy, sizeof(float));
        }
        out << endl;
    }
    out.close();

    UBLOG(logDEBUG1, "WbWriterVtkBinary::writeQuadsWithCellData to " << vtkfilename << " - end");
    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkBinary::writeQuadsWithNodeAndCellData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                             vector<UbTupleInt4> &cells, vector<string> &nodedatanames,
                                                             vector<vector<double>> &nodedata,
                                                             vector<string> &celldatanames,
                                                             vector<vector<double>> &celldata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkBinary::writeQuadsWithNodeAndCellData to " << vtkfilename << " - start");

    ofstream out(vtkfilename.c_str(), ofstream::out | ofstream::binary);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str(), ios::out | ios::binary);
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    // HEADER-SECTION
    // WRITE BIGENDIAN VtkBinary FILE
    bool swapByte = UbSystem::isLittleEndian();
    int nofNodes  = (int)nodes.size();
    int nofCells  = (int)cells.size();

    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "D3Q19MasterNodeGrid"
        << "\n";
    out << ""
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++) {
        float x1 = (float)val<1>(nodes[n]);
        float x2 = (float)val<2>(nodes[n]);
        float x3 = (float)val<3>(nodes[n]);

        if (swapByte) {
            UbSystem::swapByteOrder((unsigned char *)&x1, sizeof(float));
            UbSystem::swapByteOrder((unsigned char *)&x2, sizeof(float));
            UbSystem::swapByteOrder((unsigned char *)&x3, sizeof(float));
        }

        out.write((char *)&x1, sizeof(float));
        out.write((char *)&x2, sizeof(float));
        out.write((char *)&x3, sizeof(float));
    }
    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofCells << " " << nofCells * 5 << "\n";

    int nodesPerCellDummy = 4; // nofNodesPerCell
    if (swapByte)
        UbSystem::swapByteOrder((unsigned char *)&nodesPerCellDummy, sizeof(int));
    for (int c = 0; c < (int)cells.size(); c++) {
        int SW = val<1>(cells[c]);
        int SE = val<2>(cells[c]);
        int NE = val<3>(cells[c]);
        int NW = val<4>(cells[c]);
        if (swapByte) {
            UbSystem::swapByteOrder((unsigned char *)&SW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&SE, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&NW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&NE, sizeof(int));
        }

        out.write((char *)&nodesPerCellDummy, sizeof(int));
        out.write((char *)&SW, sizeof(int));
        out.write((char *)&SE, sizeof(int));
        out.write((char *)&NW, sizeof(int));
        out.write((char *)&NE, sizeof(int));
    }
    out << "\n";

    out << "CELL_TYPES " << (int)cells.size() << "\n";
    int celltype = 8;
    if (swapByte)
        UbSystem::swapByteOrder((unsigned char *)&celltype, sizeof(int));
    for (int c = 0; c < nofCells; c++)
        out.write((char *)&celltype, sizeof(int));

    out << endl;

    // NODE DATA SECTION
    // write data section
    out << "POINT_DATA " << nofNodes << "\n";
    for (int s = 0; s < (int)nodedatanames.size(); ++s) {
        if ((int)nodedata[s].size() != nofNodes)
            throw UbException(UB_EXARGS, "datasetsize must be equal to nofNodes");
        out << "SCALARS " << nodedatanames[s] << " float 1 \n LOOKUP_TABLE default \n";
        for (int d = 0; d < (int)nodedata[s].size(); d++) {
            float dummy = (float)nodedata[s][d];
            if (swapByte)
                UbSystem::swapByteOrder((unsigned char *)&dummy, sizeof(float));
            out.write((const char *)&dummy, sizeof(float));
        }
        out << endl;
    }

    // CELL DATA SECTION
    // write data section
    out << "CELL_DATA " << nofCells << "\n";
    for (int s = 0; s < (int)celldatanames.size(); ++s) {
        if ((int)celldata[s].size() != nofCells)
            throw UbException(UB_EXARGS, "datasetsize must be equal to nofNodes");
        out << "SCALARS " << celldatanames[s] << " float 1 \n LOOKUP_TABLE default \n";
        for (int d = 0; d < (int)celldata[s].size(); d++) {
            float dummy = (float)celldata[s][d];
            if (swapByte)
                UbSystem::swapByteOrder((unsigned char *)&dummy, sizeof(float));
            out.write((const char *)&dummy, sizeof(float));
        }
        out << endl;
    }

    out.close();

    UBLOG(logDEBUG1, "WbWriterVtkBinary::writeQuadsWithNodeAndCellData to " << vtkfilename << " - end");
    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkBinary::writeOctsWithCellData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                     vector<UbTupleInt8> &cells, vector<string> &datanames,
                                                     vector<vector<double>> &celldata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkBinary::writeOctsWithCellData to " << vtkfilename << " - start");

    ofstream out(vtkfilename.c_str(), ofstream::out | ofstream::binary);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str(), ios::out | ios::binary);
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    // HEADER-SECTION
    // WRITE BIGENDIAN VtkBinary FILE
    bool swapByte = UbSystem::isLittleEndian();
    int nofNodes  = (int)nodes.size();
    int nofCells  = (int)cells.size();

    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "D3Q19MasterNodeGrid"
        << "\n";
    out << ""
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++) {
        float x1 = (float)val<1>(nodes[n]);
        float x2 = (float)val<2>(nodes[n]);
        float x3 = (float)val<3>(nodes[n]);

        if (swapByte) {
            UbSystem::swapByteOrder((unsigned char *)&x1, sizeof(float));
            UbSystem::swapByteOrder((unsigned char *)&x2, sizeof(float));
            UbSystem::swapByteOrder((unsigned char *)&x3, sizeof(float));
        }

        out.write((char *)&x1, sizeof(float));
        out.write((char *)&x2, sizeof(float));
        out.write((char *)&x3, sizeof(float));
    }
    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofCells << " " << nofCells * 9 << "\n";

    int nodesPerCellDummy = 8; // nofNodesPerCell
    if (swapByte)
        UbSystem::swapByteOrder((unsigned char *)&nodesPerCellDummy, sizeof(int));
    for (int c = 0; c < (int)cells.size(); c++) {
        int BSW = val<1>(cells[c]);
        int TSW = val<5>(cells[c]);
        int BSE = val<2>(cells[c]);
        int TSE = val<6>(cells[c]);
        int BNW = val<3>(cells[c]);
        int TNW = val<7>(cells[c]);
        int BNE = val<4>(cells[c]);
        int TNE = val<8>(cells[c]);
        if (swapByte) {
            UbSystem::swapByteOrder((unsigned char *)&BSW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&BSE, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&BNW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&BNE, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&TSW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&TSE, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&TNW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&TNE, sizeof(int));
        }

        out.write((char *)&nodesPerCellDummy, sizeof(int));
        out.write((char *)&BSW, sizeof(int));
        out.write((char *)&BSE, sizeof(int));
        out.write((char *)&BNE, sizeof(int));
        out.write((char *)&BNW, sizeof(int));
        out.write((char *)&TSW, sizeof(int));
        out.write((char *)&TSE, sizeof(int));
        out.write((char *)&TNE, sizeof(int));
        out.write((char *)&TNW, sizeof(int));
    }
    out << "\n";

    out << "CELL_TYPES " << (int)cells.size() << "\n";
    int celltype = 11;
    if (swapByte)
        UbSystem::swapByteOrder((unsigned char *)&celltype, sizeof(int));
    for (int c = 0; c < nofCells; c++)
        out.write((char *)&celltype, sizeof(int));

    out << endl;

    // CELL DATA SECTION
    // write data section
    out << "CELL_DATA " << nofCells << "\n";
    for (int s = 0; s < (int)datanames.size(); ++s) {
        if ((int)celldata[s].size() != nofCells)
            throw UbException(UB_EXARGS, "datasetsize must be equal to nofNodes");
        out << "SCALARS " << datanames[s] << " float 1 \n LOOKUP_TABLE default \n";
        for (int d = 0; d < (int)celldata[s].size(); d++) {
            float dummy = (float)celldata[s][d];
            if (swapByte)
                UbSystem::swapByteOrder((unsigned char *)&dummy, sizeof(float));
            out.write((const char *)&dummy, sizeof(float));
        }
        out << endl;
    }
    out.close();

    UBLOG(logDEBUG1, "WbWriterVtkBinary::writeOctsWithCellData to " << vtkfilename << " - end");
    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkBinary::writeOctsWithNodeData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                     vector<UbTupleUInt8> &cells, vector<string> &datanames,
                                                     vector<vector<double>> &nodedata)
{
    // HEADER-SECTION
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkBinary::writeOctsWithNodeData to " << vtkfilename << " - start");

    ofstream out(vtkfilename.c_str(), ofstream::out | ofstream::binary);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str(), ios::out | ios::binary);
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    // WRITE BIGENDIAN VtkBinary FILE
    bool swapByte = UbSystem::isLittleEndian();
    int nofNodes  = (int)nodes.size();
    int nofCells  = (int)cells.size();

    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "D3Q19MasterNodeGrid"
        << "\n";
    out << ""
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++) {
        float x1 = val<1>(nodes[n]);
        float x2 = val<2>(nodes[n]);
        float x3 = val<3>(nodes[n]);

        if (swapByte) {
            UbSystem::swapByteOrder((unsigned char *)&x1, sizeof(float));
            UbSystem::swapByteOrder((unsigned char *)&x2, sizeof(float));
            UbSystem::swapByteOrder((unsigned char *)&x3, sizeof(float));
        }

        out.write((char *)&x1, sizeof(float));
        out.write((char *)&x2, sizeof(float));
        out.write((char *)&x3, sizeof(float));
    }
    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofCells << " " << nofCells * 9 << "\n";

    int nodesPerCellDummy = 8; // nofNodesPerCell
    if (swapByte)
        UbSystem::swapByteOrder((unsigned char *)&nodesPerCellDummy, sizeof(int));
    for (int c = 0; c < (int)cells.size(); c++) {
        int BSW = val<1>(cells[c]);
        int TSW = val<5>(cells[c]);
        int BSE = val<2>(cells[c]);
        int TSE = val<6>(cells[c]);
        int BNW = val<3>(cells[c]);
        int TNW = val<7>(cells[c]);
        int BNE = val<4>(cells[c]);
        int TNE = val<8>(cells[c]);
        if (swapByte) {
            UbSystem::swapByteOrder((unsigned char *)&BSW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&BSE, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&BNW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&BNE, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&TSW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&TSE, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&TNW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&TNE, sizeof(int));
        }

        out.write((char *)&nodesPerCellDummy, sizeof(int));
        out.write((char *)&BSW, sizeof(int));
        out.write((char *)&BSE, sizeof(int));
        out.write((char *)&BNE, sizeof(int));
        out.write((char *)&BNW, sizeof(int));
        out.write((char *)&TSW, sizeof(int));
        out.write((char *)&TSE, sizeof(int));
        out.write((char *)&TNE, sizeof(int));
        out.write((char *)&TNW, sizeof(int));
    }
    out << "\n";

    out << "CELL_TYPES " << (int)cells.size() << "\n";
    int celltype = 11;
    if (swapByte)
        UbSystem::swapByteOrder((unsigned char *)&celltype, sizeof(int));
    for (int c = 0; c < nofCells; c++)
        out.write((char *)&celltype, sizeof(int));

    out << endl;

    // NODE DATA SECTION
    // write data section
    out << "POINT_DATA " << nofNodes << "\n";
    for (int s = 0; s < (int)datanames.size(); ++s) {
        if ((int)nodedata[s].size() != nofNodes)
            throw UbException(UB_EXARGS, "datasetsize must be equal to nofNodes");
        out << "SCALARS " << datanames[s] << " float 1 \n LOOKUP_TABLE default \n";
        for (int d = 0; d < (int)nodedata[s].size(); d++) {
            float dummy = (float)nodedata[s][d];
            if (swapByte)
                UbSystem::swapByteOrder((unsigned char *)&dummy, sizeof(float));
            out.write((const char *)&dummy, sizeof(float));
        }
        out << endl;
    }

    out.close();

    UBLOG(logDEBUG1, "WbWriterVtkBinary::writeOctsWithNodeData to " << vtkfilename << " - end");
    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkBinary::writeOcts(const string &filename, vector<UbTupleFloat3> &nodes,
                                         vector<UbTupleInt8> &cells)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkBinary::writeOcts to " << vtkfilename << " - start");

    ofstream out(vtkfilename.c_str(), ofstream::out | ofstream::binary);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str(), ios::out | ios::binary);
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    // HEADER-SECTION
    // WRITE BIGENDIAN VtkBinary FILE
    bool swapByte = UbSystem::isLittleEndian();
    int nofNodes  = (int)nodes.size();
    int nofCells  = (int)cells.size();

    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "D3Q19MasterNodeGrid"
        << "\n";
    out << ""
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++) {
        float x1 = val<1>(nodes[n]);
        float x2 = val<2>(nodes[n]);
        float x3 = val<3>(nodes[n]);

        if (swapByte) {
            UbSystem::swapByteOrder((unsigned char *)&x1, sizeof(float));
            UbSystem::swapByteOrder((unsigned char *)&x2, sizeof(float));
            UbSystem::swapByteOrder((unsigned char *)&x3, sizeof(float));
        }

        out.write((char *)&x1, sizeof(float));
        out.write((char *)&x2, sizeof(float));
        out.write((char *)&x3, sizeof(float));
    }
    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofCells << " " << nofCells * 9 << "\n";

    int nodesPerCellDummy = 8; // nofNodesPerCell
    if (swapByte)
        UbSystem::swapByteOrder((unsigned char *)&nodesPerCellDummy, sizeof(int));
    for (int c = 0; c < (int)cells.size(); c++) {
        int BSW = val<1>(cells[c]);
        int TSW = val<5>(cells[c]);
        int BSE = val<2>(cells[c]);
        int TSE = val<6>(cells[c]);
        int BNW = val<3>(cells[c]);
        int TNW = val<7>(cells[c]);
        int BNE = val<4>(cells[c]);
        int TNE = val<8>(cells[c]);
        if (swapByte) {
            UbSystem::swapByteOrder((unsigned char *)&BSW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&BSE, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&BNW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&BNE, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&TSW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&TSE, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&TNW, sizeof(int));
            UbSystem::swapByteOrder((unsigned char *)&TNE, sizeof(int));
        }

        out.write((char *)&nodesPerCellDummy, sizeof(int));
        out.write((char *)&BSW, sizeof(int));
        out.write((char *)&BSE, sizeof(int));
        out.write((char *)&BNE, sizeof(int));
        out.write((char *)&BNW, sizeof(int));
        out.write((char *)&TSW, sizeof(int));
        out.write((char *)&TSE, sizeof(int));
        out.write((char *)&TNE, sizeof(int));
        out.write((char *)&TNW, sizeof(int));
    }
    out << "\n";

    out << "CELL_TYPES " << (int)cells.size() << "\n";
    int celltype = 11;
    if (swapByte)
        UbSystem::swapByteOrder((unsigned char *)&celltype, sizeof(int));
    for (int c = 0; c < nofCells; c++)
        out.write((char *)&celltype, sizeof(int));

    out << endl;
    out.close();

    UBLOG(logDEBUG1, "WbWriterVtkBinary::writeOcts to " << vtkfilename << " - end");
    return vtkfilename;
}

//! \}
