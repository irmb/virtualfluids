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
#include <cstring>

using namespace std;

std::string WbWriterVtkASCII::writeQuads(const string &filename, vector<UbTupleFloat3> &nodes,
                                         vector<UbTupleInt4> &cells)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeQuads to " << vtkfilename << " - start");

    std::ofstream out(vtkfilename.c_str());
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str());
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    // VtkASCII FILE
    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "GeoFile"
        << "\n";
    out << "ASCII"
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << " \n";

    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofCells << " " << 5 * nofCells << "\n";
    for (int c = 0; c < (int)cells.size(); c++)
        out << "4 " << val<1>(cells[c]) << " " << val<2>(cells[c]) << " " << val<4>(cells[c]) << " " << val<3>(cells[c])
            << " \n";
    out << "\n";

    out << "CELL_TYPES " << nofCells << "\n";
    for (int i = 0; i < nofCells; i++)
        out << "8" << endl;
    out << endl;

    out.close();

    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeQuads to " << vtkfilename << " - end");

    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkASCII::writeQuadsWithNodeData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                     vector<UbTupleInt4> &cells, vector<string> &datanames,
                                                     vector<vector<double>> &nodedata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeQuadsWithNodeData to " << vtkfilename << " - start");

    std::ofstream out(vtkfilename.c_str());
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str());
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    // write geo
    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    // VtkASCII FILE
    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "GeoFile"
        << "\n";
    out << "ASCII"
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << " \n";

    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofCells << " " << 5 * nofCells << "\n";
    for (int c = 0; c < (int)cells.size(); c++)
        out << "4 " << val<1>(cells[c]) << " " << val<2>(cells[c]) << " " << val<4>(cells[c]) << " " << val<3>(cells[c])
            << " \n";
    out << "\n";

    out << "CELL_TYPES " << nofCells << "\n";
    for (int i = 0; i < nofCells; i++)
        out << "8" << endl;
    out << endl;

    // write data section
    out << "POINT_DATA " << nofNodes << "\n";
    for (int s = 0; s < (int)datanames.size(); ++s) {
        out << "SCALARS " << datanames[s] << " float 1 \n LOOKUP_TABLE default \n";
        for (int d = 0; d < (int)nodedata[s].size(); d++)
            out << nodedata[s][d] << "\n";

        out << endl;
    }

    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeQuadsWithNodeData to " << vtkfilename << " - end");
    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkASCII::writeQuadsWithCellData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                     vector<UbTupleInt4> &cells, vector<string> &datanames,
                                                     vector<vector<double>> &celldata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeQuadsWithCellData to " << vtkfilename << " - start");

    std::ofstream out(vtkfilename.c_str());
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str());
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    // write geo
    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    // VtkASCII FILE
    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "GeoFile"
        << "\n";
    out << "ASCII"
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << " \n";

    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofCells << " " << 5 * nofCells << "\n";
    for (int c = 0; c < (int)cells.size(); c++)
        out << "4 " << val<1>(cells[c]) << " " << val<2>(cells[c]) << " " << val<4>(cells[c]) << " " << val<3>(cells[c])
            << " \n";
    out << "\n";

    out << "CELL_TYPES " << nofCells << "\n";
    for (int i = 0; i < nofCells; i++)
        out << "8" << endl;
    out << endl;

    // write data section
    out << "CELL_DATA " << nofCells << "\n";
    for (int s = 0; s < (int)datanames.size(); ++s) {
        out << "SCALARS " << datanames[s] << " float 1 \n LOOKUP_TABLE default \n";
        for (int d = 0; d < (int)celldata[s].size(); d++)
            out << celldata[s][d] << "\n";

        out << endl;
    }

    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeQuadsWithCellData to " << vtkfilename << " - end");
    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkASCII::writeQuadsWithNodeAndCellData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                            vector<UbTupleInt4> &cells, vector<string> &nodedatanames,
                                                            vector<vector<double>> &nodedata,
                                                            vector<string> &celldatanames,
                                                            vector<vector<double>> &celldata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeQuadsWithNodeAndCellData to " << vtkfilename << " - start");

    std::ofstream out(vtkfilename.c_str());
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str());
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    // write geo
    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    // VtkASCII FILE
    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "GeoFile"
        << "\n";
    out << "ASCII"
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << " \n";

    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofCells << " " << 5 * nofCells << "\n";
    for (int c = 0; c < (int)cells.size(); c++)
        out << "4 " << val<1>(cells[c]) << " " << val<2>(cells[c]) << " " << val<4>(cells[c]) << " " << val<3>(cells[c])
            << " \n";
    out << "\n";

    out << "CELL_TYPES " << nofCells << "\n";
    for (int i = 0; i < nofCells; i++)
        out << "8" << endl;
    out << endl;

    // write node data section
    out << "POINT_DATA " << nofNodes << "\n";
    for (int s = 0; s < (int)nodedatanames.size(); ++s) {
        out << "SCALARS " << nodedatanames[s] << " float 1 \n LOOKUP_TABLE default \n";
        for (int d = 0; d < (int)nodedata[s].size(); d++)
            out << nodedata[s][d] << "\n";

        out << endl;
    }

    // write cell data section
    out << "CELL_DATA " << nofCells << "\n";
    for (int s = 0; s < (int)celldatanames.size(); ++s) {
        out << "SCALARS " << celldatanames[s] << " float 1 \n LOOKUP_TABLE default \n";
        for (int d = 0; d < (int)celldata[s].size(); d++)
            out << celldata[s][d] << "\n";

        out << endl;
    }

    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeQuadsWithNodeAndCellData to " << vtkfilename << " - end");
    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkASCII::writeLines(const string &filename, vector<UbTupleFloat3> &nodes,
                                         vector<UbTupleInt2> &lines)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeLines to " << vtkfilename << " - start");

    std::ofstream out(vtkfilename.c_str());
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str());
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    int nofNodes = (int)nodes.size();
    int nofLines = (int)lines.size();

    // VtkASCII FILE
    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "GeoFile"
        << "\n";
    out << "ASCII"
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++) {
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << " \n";
    }
    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofLines << " " << 3 * nofLines << "\n";
    int nr = 0;
    for (int l = 0; l < nofLines; l++) {
        int el = nr + 1;
        out << "2 " << val<1>(lines[l]) << " " << val<2>(lines[l]) << " " << endl;
        nr = el + 1;
    }
    out << "\n";

    out << "CELL_TYPES " << nofLines << "\n";
    for (int l = 0; l < nofLines; l++)
        out << "3" << endl;
    out << endl;

    out.close();

    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeLines to " << vtkfilename << " - end");
    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkASCII::writeTriangles(const string &filename, vector<UbTupleFloat3> &nodes,
                                             vector<UbTupleInt3> &triangles)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeTriangles to " << vtkfilename << " - start");

    std::ofstream out(vtkfilename.c_str());
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str());
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    int nofNodes     = (int)nodes.size();
    int nofTriangles = (int)triangles.size();

    // VtkASCII FILE
    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "GeoFile"
        << "\n";
    out << "ASCII"
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++) {
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << " \n";
    }
    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofTriangles << " " << 4 * nofTriangles << "\n";
    int nr = 0;
    for (int t = 0; t < nofTriangles; t++) {
        int el = nr + 1;
        out << "3 " << val<1>(triangles[t]) << " " << val<2>(triangles[t]) << " " << val<3>(triangles[t]) << " "
            << endl;
        nr = el + 1;
    }
    out << "\n";

    out << "CELL_TYPES " << nofTriangles << "\n";
    for (int l = 0; l < nofTriangles; l++)
        out << "5" << endl;
    out << endl;

    out.close();

    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeTriangles to " << vtkfilename << " - end");
    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkASCII::writeTrianglesWithNodeData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                         vector<UbTupleInt3> &cells, vector<string> &datanames,
                                                         vector<vector<double>> &nodedata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeTrianglesWithNodeData to " << vtkfilename << " - start");

    std::ofstream out(vtkfilename.c_str());
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str());
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    // write geo
    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    // VtkASCII FILE
    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "GeoFile"
        << "\n";
    out << "ASCII"
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << " \n";

    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofCells << " " << 4 * nofCells << "\n";
    for (int c = 0; c < (int)cells.size(); c++)
        out << "3 " << val<1>(cells[c]) << " " << val<2>(cells[c]) << " " << val<3>(cells[c]) << " \n";
    out << "\n";

    out << "CELL_TYPES " << nofCells << "\n";
    for (int i = 0; i < nofCells; i++)
        out << "5" << endl;
    out << endl;

    // write data section
    out << "POINT_DATA " << nofNodes << "\n";
    for (int s = 0; s < (int)datanames.size(); ++s) {
        out << "SCALARS " << datanames[s] << " float 1 \n LOOKUP_TABLE default \n";
        for (int d = 0; d < (int)nodedata[s].size(); d++)
            out << nodedata[s][d] << "\n";

        out << endl;
    }

    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeTrianglesWithNodeData to " << vtkfilename << " - end");
    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkASCII::writeOctsWithCellData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                    vector<UbTupleInt8> &cells, vector<string> &datanames,
                                                    vector<vector<double>> &celldata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeOctsWithCellData to " << vtkfilename << " - start");

    std::ofstream out(vtkfilename.c_str());
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str());
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    // write geo
    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    // VtkASCII FILE
    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "GeoFile"
        << "\n";
    out << "ASCII"
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << " \n";

    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofCells << " " << 9 * nofCells << "\n";
    for (int c = 0; c < (int)cells.size(); c++) {
        out << "8 " << val<1>(cells[c]) << " " << val<2>(cells[c]) << " " << val<4>(cells[c]) << " " << val<3>(cells[c])
            << " " << val<5>(cells[c]) << " " << val<6>(cells[c]) << " " << val<8>(cells[c]) << " " << val<7>(cells[c])
            << " \n";
    }

    out << "\n";

    out << "CELL_TYPES " << nofCells << "\n";
    for (int i = 0; i < nofCells; i++)
        out << "11 " << endl;
    out << endl;

    // write data section
    out << "CELL_DATA " << nofCells << "\n";
    for (int s = 0; s < (int)datanames.size(); ++s) {
        out << "SCALARS " << datanames[s] << " float 1 \n LOOKUP_TABLE default \n";
        for (int d = 0; d < (int)celldata[s].size(); d++)
            out << celldata[s][d] << "\n";

        out << endl;
    }

    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeOctsWithCellData to " << vtkfilename << " - end");

    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkASCII::writeOctsWithNodeData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                    vector<UbTupleUInt8> &cells, vector<string> &datanames,
                                                    vector<vector<double>> &nodedata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeOctsWithNodeData to " << vtkfilename << " - start");

    std::ofstream out(vtkfilename.c_str());
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str());
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    // write geo
    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    // VtkASCII FILE
    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "GeoFile"
        << "\n";
    out << "ASCII"
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << " \n";

    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofCells << " " << 9 * nofCells << "\n";
    for (int c = 0; c < (int)cells.size(); c++) {
        out << "8 " << val<1>(cells[c]) << " " << val<2>(cells[c]) << " " << val<4>(cells[c]) << " " << val<3>(cells[c])
            << " " << val<5>(cells[c]) << " " << val<6>(cells[c]) << " " << val<8>(cells[c]) << " " << val<7>(cells[c])
            << " \n";
    }
    out << "\n";

    out << "CELL_TYPES " << nofCells << "\n";
    for (int i = 0; i < nofCells; i++)
        out << "11" << endl;
    out << endl;

    // write data section
    out << "POINT_DATA " << nofNodes << "\n";
    for (int s = 0; s < (int)datanames.size(); ++s) {
        out << "SCALARS " << datanames[s] << " float 1 \n LOOKUP_TABLE default \n";
        for (int d = 0; d < (int)nodedata[s].size(); d++)
            out << nodedata[s][d] << "\n";

        out << endl;
    }

    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeOctsWithNodeData to " << vtkfilename << " - end");
    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkASCII::writeOcts(const string &filename, vector<UbTupleFloat3> &nodes,
                                        vector<UbTupleInt8> &cells)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeOcts to " << vtkfilename << " - start");

    std::ofstream out(vtkfilename.c_str());
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(vtkfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(vtkfilename.c_str());
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + vtkfilename);
    }

    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    // VtkASCII FILE
    out << "# vtk DataFile Version 4.0"
        << "\n";
    out << "GeoFile"
        << "\n";
    out << "ASCII"
        << "\n";

    // POINTS SECTION
    out << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    out << "POINTS " << nofNodes << " float"
        << "\n";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << " \n";

    out << "\n";

    // CELLS SECTION
    out << "CELLS " << nofCells << " " << 9 * nofCells << "\n";
    for (int c = 0; c < (int)cells.size(); c++)
        out << "8 " << val<1>(cells[c]) << " " << val<2>(cells[c]) << " " << val<4>(cells[c]) << " " << val<3>(cells[c])
            << " " << val<5>(cells[c]) << " " << val<6>(cells[c]) << " " << val<8>(cells[c]) << " " << val<7>(cells[c])
            << " \n";
    out << "\n";

    out << "CELL_TYPES " << nofCells << "\n";
    for (int i = 0; i < nofCells; i++)
        out << "11" << endl;
    out << endl;

    out.close();

    UBLOG(logDEBUG1, "WbWriterVtkASCII::writeOcts to " << vtkfilename << " - end");
    return vtkfilename;
}

//! \}
