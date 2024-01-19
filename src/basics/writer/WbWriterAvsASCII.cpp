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
#include <basics/utilities/UbSystem.h>
#include <basics/writer/WbWriterAvsASCII.h>
#include <cstring>

using namespace std;

std::string WbWriterAvsASCII::writeQuads(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                         std::vector<UbTupleInt4> &cells)
{
    string avsfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeQuads to " << avsfilename << " - start");

    ofstream out(avsfilename.c_str(), ios::out | ios::binary);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
        string path = UbSystem::getPathFromString(avsfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(avsfilename.c_str(), ios::out | ios::binary);
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + avsfilename);
    }
    char magic = (char)7;
    int idummy;
    float fdummy;

    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    int nofNodeData     = 0;
    int nofCellData     = 0;
    int nofModelData    = 0;
    int cellType        = 3; //=quad
    int nofNodesPerCell = 4;

    out.write((char *)&magic, sizeof(char));
    out.write((char *)&nofNodes, sizeof(int));
    out.write((char *)&nofCells, sizeof(int));
    out.write((char *)&nofNodeData, sizeof(int));
    out.write((char *)&nofCellData, sizeof(int));
    out.write((char *)&nofModelData, sizeof(int));

    idummy = (int)nofCells * nofNodesPerCell;
    out.write((char *)&idummy, sizeof(int)); //(nof nodes) * (nodes per cell)
    for (int c = 0; c < nofCells; c++) {
        idummy = c + 1;
        out.write((char *)&idummy, sizeof(int)); // cell id
        idummy = 1;
        out.write((char *)&idummy, sizeof(int)); // mat
        idummy = nofNodesPerCell;
        out.write((char *)&idummy, sizeof(int)); // nodes per cell
        idummy = cellType;
        out.write((char *)&idummy, sizeof(int)); // cell type
    }
    // knotennummern der einzelnen zellen
    for (int c = 0; c < nofCells; c++) {
        idummy = val<1>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<2>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<3>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<4>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
    }

    // coords
    // x1-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<1>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // x2-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<2>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // x3-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<3>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // out<<"\n";

    out.close();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeQuads to " << avsfilename << " - end");

    return avsfilename;
}
/*===============================================================================*/
std::string WbWriterAvsASCII::writeOcts(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                        std::vector<UbTupleInt8> &cells)
{
    string avsfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeOcts to " << avsfilename << " - start");

    ofstream out(avsfilename.c_str(), ios::out | ios::binary);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
        string path = UbSystem::getPathFromString(avsfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(avsfilename.c_str(), ios::out | ios::binary);
        }
        if (!out)
            throw UbException(UB_EXARGS, "file konnte nicht geschrieben werden " + avsfilename);
    }

    char magic = (char)7;
    int idummy;
    float fdummy;

    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    int nofNodeData     = 0;
    int nofCellData     = 0;
    int nofModelData    = 0;
    int cellType        = 7; //=hex
    int nofNodesPerCell = 8;

    out.write((char *)&magic, sizeof(char));
    out.write((char *)&nofNodes, sizeof(int));
    out.write((char *)&nofCells, sizeof(int));
    out.write((char *)&nofNodeData, sizeof(int));
    out.write((char *)&nofCellData, sizeof(int));
    out.write((char *)&nofModelData, sizeof(int));

    idummy = (int)nofCells * nofNodesPerCell;
    out.write((char *)&idummy, sizeof(int)); //(nof nodes) * (nodes per cell)
    for (int c = 0; c < nofCells; c++) {
        idummy = c + 1;
        out.write((char *)&idummy, sizeof(int)); // cell id
        idummy = 1;
        out.write((char *)&idummy, sizeof(int)); // mat
        idummy = nofNodesPerCell;
        out.write((char *)&idummy, sizeof(int)); // nodes per cell
        idummy = cellType;
        out.write((char *)&idummy, sizeof(int)); // cell type
    }
    // knotennummern der einzelnen zellen
    for (int c = 0; c < nofCells; c++) {
        idummy = val<1>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<2>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<3>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<4>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<5>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<6>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<7>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<8>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
    }

    // coords
    // x1-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<1>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // x2-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<2>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // x3-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<3>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // out<<"\n";

    out.close();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeOcts to " << avsfilename << " - end");

    return avsfilename;
}
/*===============================================================================*/
std::string WbWriterAvsASCII::writeQuadsWithNodeData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                     vector<UbTupleInt4> &cells, vector<string> &datanames,
                                                     vector<vector<double>> &nodedata)
{
    string avsfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeQuadsWithNodeData to " << avsfilename << " - start");

    ofstream out(avsfilename.c_str(), ios::out | ios::binary);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
        string path = UbSystem::getPathFromString(avsfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(avsfilename.c_str(), ios::out | ios::binary);
        }
        if (!out)
            throw UbException(UB_EXARGS, "write UCD File " + avsfilename + " konnte nicht geschrieben werden");
    }
    char magic = (char)7;
    int idummy;
    float fdummy;

    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    int nofNodeData     = (int)datanames.size();
    int nofCellData     = 0;
    int nofModelData    = 0;
    int cellType        = 3; //=quad
    int nofNodesPerCell = 4;

    out.write((char *)&magic, sizeof(char));
    out.write((char *)&nofNodes, sizeof(int));
    out.write((char *)&nofCells, sizeof(int));
    out.write((char *)&nofNodeData, sizeof(int));
    out.write((char *)&nofCellData, sizeof(int));
    out.write((char *)&nofModelData, sizeof(int));

    idummy = (int)nofCells * nofNodesPerCell;
    out.write((char *)&idummy, sizeof(int)); //(nof nodes) * (nodes per cell)
    for (int c = 0; c < nofCells; c++) {
        idummy = c + 1;
        out.write((char *)&idummy, sizeof(int)); // cell id
        idummy = 1;
        out.write((char *)&idummy, sizeof(int)); // mat
        idummy = nofNodesPerCell;
        out.write((char *)&idummy, sizeof(int)); // nodes per cell
        idummy = cellType;
        out.write((char *)&idummy, sizeof(int)); // cell type
    }
    // knotennummern der einzelnen zellen
    for (int c = 0; c < nofCells; c++) {
        idummy = val<1>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<2>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<3>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<4>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
    }

    // coords
    // x1-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<1>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // x2-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<2>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // x3-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<3>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // out<<"\n";

    // NODE DATA
    char labels[1024];
    char units[1024];
    strcpy(labels, "");
    strcpy(units, "");

    for (int d = 0; d < nofNodeData - 1; ++d) {
        strcat(labels, datanames[d].c_str());
        strcat(labels, ".");
    }
    strcat(labels, datanames[nofNodeData - 1].c_str());

    for (int i = 0; i < (nofNodeData - 1); i++)
        strcat(units, "no_unit.");
    strcat(units, "no_unit");

    out.write((char *)&labels, sizeof(labels));
    out.write((char *)&units, sizeof(units));

    // nof and type of data
    idummy = nofNodeData;
    out.write((char *)&idummy, sizeof(int)); // Datentypen pro knoten (hier = nof_node_data, da NUR skalare)

    idummy = 1;
    for (int i = 0; i < nofNodeData; i++)
        out.write((char *)&idummy, sizeof(int)); // jeder Datentyp ist ein skalarer Wert

    // min and max of data
    fdummy = 0.0;
    for (int i = 0; i < nofNodeData; i++)
        out.write((char *)&fdummy, sizeof(float)); // min Wert pro Datentyp
    fdummy = 1.0;
    for (int i = 0; i < nofNodeData; i++)
        out.write((char *)&fdummy, sizeof(float)); // max Wert pro Datentyp

    // daten ins file schreiben
    for (int d = 0; d < nofNodeData; ++d)
        for (int n = 0; n < (int)nodedata[d].size(); n++) {
            fdummy = (float)nodedata[d][n];
            out.write((char *)&fdummy, sizeof(float));
        }

    fdummy = 1.;
    for (int i = 0; i < nofNodeData; i++)
        out.write((char *)&fdummy, sizeof(float)); // max Wert pro Datentyp

    out.close();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeQuadsWithNodeData to " << avsfilename << " - end");

    return avsfilename;
}
/*===============================================================================*/
std::string WbWriterAvsASCII::writeQuadsWithCellData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                     vector<UbTupleInt4> &cells, vector<string> &datanames,
                                                     vector<vector<double>> &celldata)
{
    string avsfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeQuadsWithCellData to " << avsfilename << " - start");

    ofstream out(avsfilename.c_str(), ios::out | ios::binary);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
        string path = UbSystem::getPathFromString(avsfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(avsfilename.c_str(), ios::out | ios::binary);
        }
        if (!out)
            throw UbException(UB_EXARGS,
                              "write_OutputFile-UCD File  " + avsfilename + " konnte nicht geschrieben werden");
    }

    char magic = (char)7;
    int idummy;
    float fdummy;

    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    int nofNodeData     = 0;
    int nofCellData     = (int)datanames.size();
    int nofModelData    = 0;
    int cellType        = 3; //=quad
    int nofNodesPerCell = 4;

    out.write((char *)&magic, sizeof(char));
    out.write((char *)&nofNodes, sizeof(int));
    out.write((char *)&nofCells, sizeof(int));
    out.write((char *)&nofNodeData, sizeof(int));
    out.write((char *)&nofCellData, sizeof(int));
    out.write((char *)&nofModelData, sizeof(int));

    idummy = (int)nofCells * nofNodesPerCell;
    out.write((char *)&idummy, sizeof(int)); //(nof nodes) * (nodes per cell)
    for (int c = 0; c < nofCells; c++) {
        idummy = c + 1;
        out.write((char *)&idummy, sizeof(int)); // cell id
        idummy = 1;
        out.write((char *)&idummy, sizeof(int)); // mat
        idummy = nofNodesPerCell;
        out.write((char *)&idummy, sizeof(int)); // nodes per cell
        idummy = cellType;
        out.write((char *)&idummy, sizeof(int)); // cell type
    }
    // knotennummern der einzelnen zellen
    for (int c = 0; c < nofCells; c++) {
        idummy = val<1>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<2>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<3>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<4>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
    }

    // coords
    // x1-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<1>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // x2-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<2>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // x3-coords
    fdummy = 0.0;
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<3>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // out<<"\n";

    // CELL DATA
    char labels[1024];
    char units[1024];
    strcpy(labels, "");
    strcpy(units, "");

    for (int d = 0; d < nofCellData - 1; ++d) {
        strcat(labels, datanames[d].c_str());
        strcat(labels, ".");
    }
    strcat(labels, datanames[nofCellData - 1].c_str());

    for (int d = 0; d < nofCellData - 1; ++d)
        strcat(units, "no_unit.");
    strcat(units, "no_unit");

    out.write((char *)&labels, sizeof(labels));
    out.write((char *)&units, sizeof(units));

    // nof and type of data
    idummy = nofCellData;
    out.write((char *)&idummy, sizeof(int)); // Datentypen pro knoten (hier = nof_node_data, da NUR skalare)

    idummy = 1;
    for (int i = 0; i < nofCellData; i++)
        out.write((char *)&idummy, sizeof(int)); // jeder Datentyp ist ein skalarer Wert

    // min and max of data
    fdummy = 0.0;
    for (int i = 0; i < nofCellData; i++)
        out.write((char *)&fdummy, sizeof(float)); // min Wert pro Datentyp
    fdummy = 1.0;
    for (int i = 0; i < nofCellData; i++)
        out.write((char *)&fdummy, sizeof(float)); // max Wert pro Datentyp

    // daten ins file schreiben
    for (int d = 0; d < nofCellData; ++d)
        for (int n = 0; n < (int)celldata[d].size(); n++) {
            fdummy = (float)celldata[d][n];
            out.write((char *)&fdummy, sizeof(float));
        }

    fdummy = 1.;
    for (int i = 0; i < nofCellData; i++)
        out.write((char *)&fdummy, sizeof(float)); // max Wert pro Datentyp

    out.close();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeQuadsWithCellData to " << avsfilename << " - end");

    return avsfilename;
}
/*===============================================================================*/
std::string WbWriterAvsASCII::writeQuadsWithNodeAndCellData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                            vector<UbTupleInt4> &cells, vector<string> &nodedatanames,
                                                            vector<vector<double>> &nodedata,
                                                            vector<string> &celldatanames,
                                                            vector<vector<double>> &celldata)
{
    string avsfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeQuadsWithNodeAndCellData to " << avsfilename << " - start");

    ofstream out(avsfilename.c_str(), ios::out | ios::binary);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
        string path = UbSystem::getPathFromString(avsfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(avsfilename.c_str(), ios::out | ios::binary);
        }
        if (!out)
            throw UbException(UB_EXARGS,
                              "write_OutputFile-UCD File  " + avsfilename + " konnte nicht geschrieben werden");
    }

    char magic = (char)7;
    int idummy;
    float fdummy;

    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    int nofNodeData     = (int)nodedatanames.size();
    int nofCellData     = (int)celldatanames.size();
    int nofModelData    = 0;
    int cellType        = 3; //=quad
    int nofNodesPerCell = 4;

    out.write((char *)&magic, sizeof(char));
    out.write((char *)&nofNodes, sizeof(int));
    out.write((char *)&nofCells, sizeof(int));
    out.write((char *)&nofNodeData, sizeof(int));
    out.write((char *)&nofCellData, sizeof(int));
    out.write((char *)&nofModelData, sizeof(int));

    idummy = (int)nofCells * nofNodesPerCell;
    out.write((char *)&idummy, sizeof(int)); //(nof nodes) * (nodes per cell)
    for (int c = 0; c < nofCells; c++) {
        idummy = c + 1;
        out.write((char *)&idummy, sizeof(int)); // cell id
        idummy = 1;
        out.write((char *)&idummy, sizeof(int)); // mat
        idummy = nofNodesPerCell;
        out.write((char *)&idummy, sizeof(int)); // nodes per cell
        idummy = cellType;
        out.write((char *)&idummy, sizeof(int)); // cell type
    }
    // knotennummern der einzelnen zellen
    for (int c = 0; c < nofCells; c++) {
        idummy = val<1>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<2>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<3>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<4>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
    }

    // coords
    // x1-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<1>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // x2-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<2>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // x3-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<3>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // out<<"\n";

    // NODE DATA
    char nodelabels[1024];
    char nodeunits[1024];
    strcpy(nodelabels, "");
    strcpy(nodeunits, "");

    for (int d = 0; d < nofNodeData - 1; ++d) {
        strcat(nodelabels, nodedatanames[d].c_str());
        strcat(nodelabels, ".");
    }
    strcat(nodelabels, nodedatanames[nofNodeData - 1].c_str());

    for (int i = 0; i < (nofNodeData - 1); i++)
        strcat(nodeunits, "no_unit.");
    strcat(nodeunits, "no_unit");

    out.write((char *)&nodelabels, sizeof(nodelabels));
    out.write((char *)&nodeunits, sizeof(nodeunits));

    // nof and type of data
    idummy = nofNodeData;
    out.write((char *)&idummy, sizeof(int)); // Datentypen pro knoten (hier = nof_node_data, da NUR skalare)

    idummy = 1;
    for (int i = 0; i < nofNodeData; i++)
        out.write((char *)&idummy, sizeof(int)); // jeder Datentyp ist ein skalarer Wert

    // min and max of data
    fdummy = 0.0;
    for (int i = 0; i < nofNodeData; i++)
        out.write((char *)&fdummy, sizeof(float)); // min Wert pro Datentyp
    fdummy = 1.0;
    for (int i = 0; i < nofNodeData; i++)
        out.write((char *)&fdummy, sizeof(float)); // max Wert pro Datentyp

    // daten ins file schreiben
    for (int d = 0; d < nofNodeData; ++d)
        for (int n = 0; n < (int)nodedata[d].size(); n++) {
            fdummy = (float)nodedata[d][n];
            out.write((char *)&fdummy, sizeof(float));
        }

    fdummy = 1.;
    for (int i = 0; i < nofNodeData; i++)
        out.write((char *)&fdummy, sizeof(float)); // max Wert pro Datentyp

    // CELL DATA
    char celllabels[1024];
    char cellunits[1024];
    strcpy(celllabels, "");
    strcpy(cellunits, "");

    for (int d = 0; d < nofCellData - 1; ++d) {
        strcat(celllabels, celldatanames[d].c_str());
        strcat(celllabels, ".");
    }
    strcat(celllabels, celldatanames[nofCellData - 1].c_str());

    for (int d = 0; d < nofCellData - 1; ++d)
        strcat(cellunits, "no_unit.");
    strcat(cellunits, "no_unit");

    out.write((char *)&celllabels, sizeof(celllabels));
    out.write((char *)&cellunits, sizeof(cellunits));

    // nof and type of data
    idummy = nofCellData;
    out.write((char *)&idummy, sizeof(int)); // Datentypen pro knoten (hier = nof_node_data, da NUR skalare)

    idummy = 1;
    for (int i = 0; i < nofCellData; i++)
        out.write((char *)&idummy, sizeof(int)); // jeder Datentyp ist ein skalarer Wert

    // min and max of data
    fdummy = 0.0;
    for (int i = 0; i < nofCellData; i++)
        out.write((char *)&fdummy, sizeof(float)); // min Wert pro Datentyp
    fdummy = 1.0;
    for (int i = 0; i < nofCellData; i++)
        out.write((char *)&fdummy, sizeof(float)); // max Wert pro Datentyp

    // daten ins file schreiben
    for (int d = 0; d < nofCellData; ++d)
        for (int n = 0; n < (int)celldata[d].size(); n++) {
            fdummy = (float)celldata[d][n];
            out.write((char *)&fdummy, sizeof(float));
        }

    fdummy = 1.;
    for (int i = 0; i < nofCellData; i++)
        out.write((char *)&fdummy, sizeof(float)); // max Wert pro Datentyp

    out.close();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeQuadsWithNodeAndCellData to " << avsfilename << " - end");

    return avsfilename;
}
/*===============================================================================*/
std::string WbWriterAvsASCII::writeLines(const string &filename, vector<UbTupleFloat3> &nodes,
                                         vector<UbTupleInt2> &lines)
{
    string avsfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeLines to " << avsfilename << " - start");

    ofstream out(avsfilename.c_str(), ios::out);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
        string path = UbSystem::getPathFromString(avsfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(avsfilename.c_str(), ios::out);
        }
        if (!out)
            throw UbException(UB_EXARGS, avsfilename + " konnte nicht geschrieben werden");
    }

    int nofNodes = (int)nodes.size();
    int nofLines = (int)lines.size();

    out << "# UCD-File created by WbWriterAvsASCII\n";
    out << nofNodes << " " << nofLines << " 0 0 0 " << endl;

    for (int n = 0; n < nofNodes; n++)
        out << n + 1 << " " << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << " \n";

    for (int l = 0; l < nofLines; l++)
        out << l + 1 << " 2 line " << val<1>(lines[l]) + 1 << " " << val<2>(lines[l]) + 1 << " " << endl;

    out.close();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeLines to " << avsfilename << " - end");

    return avsfilename;
}
/*===============================================================================*/
std::string WbWriterAvsASCII::writeTriangles(const string &filename, vector<UbTupleFloat3> &nodes,
                                             vector<UbTupleInt3> &triangles)
{
    string avsfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeTriangles to " << avsfilename << " - start");

    ofstream out(avsfilename.c_str(), ios::out);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
        string path = UbSystem::getPathFromString(avsfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(avsfilename.c_str(), ios::out);
        }
        if (!out)
            throw UbException(UB_EXARGS, "file konnte nicht geschrieben werden " + avsfilename);
    }

    int nofNodes = (int)nodes.size();
    int nofTrian = (int)triangles.size();

    out << "# UCD-File created by WbWriterAvsASCII\n";
    out << nofNodes << " " << nofTrian << " 0 0 0 " << endl;

    for (int n = 0; n < nofNodes; n++)
        out << n + 1 << " " << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << " \n";

    for (int l = 0; l < nofTrian; l++)
        out << l + 1 << " 2 tri " << val<1>(triangles[l]) + 1 << " " << val<2>(triangles[l]) + 1 << " "
            << val<3>(triangles[l]) + 1 << " " << endl;

    out.close();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeTriangles to " << avsfilename << " - end");

    return avsfilename;
}
/*===============================================================================*/
std::string WbWriterAvsASCII::writeTrianglesWithNodeData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                                         std::vector<UbTupleInt3> &cells,
                                                         std::vector<std::string> &datanames,
                                                         std::vector<std::vector<double>> &nodedata)
{
    string avsfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeTrianglesWithNodeData to " << avsfilename << " - end");

    ofstream out(avsfilename.c_str(), ios::out | ios::binary);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
        string path = UbSystem::getPathFromString(avsfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(avsfilename.c_str(), ios::out | ios::binary);
        }
        if (!out)
            throw UbException(UB_EXARGS,
                              "write_OutputFile-UCD File " + avsfilename + " konnte nicht geschrieben werden");
    }
    char magic = (char)7;
    int idummy;
    float fdummy;

    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    int nofNodeData     = (int)datanames.size();
    int nofCellData     = 0;
    int nofModelData    = 0;
    int cellType        = 2; // triangle
    int nofNodesPerCell = 3;

    out.write((char *)&magic, sizeof(char));
    out.write((char *)&nofNodes, sizeof(int));
    out.write((char *)&nofCells, sizeof(int));
    out.write((char *)&nofNodeData, sizeof(int));
    out.write((char *)&nofCellData, sizeof(int));
    out.write((char *)&nofModelData, sizeof(int));

    idummy = (int)nofCells * nofNodesPerCell;
    out.write((char *)&idummy, sizeof(int)); //(nof nodes) * (nodes per cell)
    for (int c = 0; c < nofCells; c++) {
        idummy = c + 1;
        out.write((char *)&idummy, sizeof(int)); // cell id
        idummy = 1;
        out.write((char *)&idummy, sizeof(int)); // mat
        idummy = nofNodesPerCell;
        out.write((char *)&idummy, sizeof(int)); // nodes per cell
        idummy = cellType;
        out.write((char *)&idummy, sizeof(int)); // cell type
    }
    // knotennummern der einzelnen zellen
    for (int c = 0; c < nofCells; c++) {
        idummy = val<1>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<2>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<3>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
    }

    // coords
    // x1-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<1>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // x2-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<2>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // x3-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<3>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // out<<"\n";

    // NODE DATA
    char labels[1024];
    char units[1024];
    strcpy(labels, "");
    strcpy(units, "");

    for (int d = 0; d < nofNodeData - 1; ++d) {
        strcat(labels, datanames[d].c_str());
        strcat(labels, ".");
    }
    strcat(labels, datanames[nofNodeData - 1].c_str());

    for (int i = 0; i < (nofNodeData - 1); i++)
        strcat(units, "no_unit.");
    strcat(units, "no_unit");

    out.write((char *)&labels, sizeof(labels));
    out.write((char *)&units, sizeof(units));

    // nof and type of data
    idummy = nofNodeData;
    out.write((char *)&idummy, sizeof(int)); // Datentypen pro knoten (hier = nof_node_data, da NUR skalare)

    idummy = 1;
    for (int i = 0; i < nofNodeData; i++)
        out.write((char *)&idummy, sizeof(int)); // jeder Datentyp ist ein skalarer Wert

    // min and max of data
    fdummy = 0.0;
    for (int i = 0; i < nofNodeData; i++)
        out.write((char *)&fdummy, sizeof(float)); // min Wert pro Datentyp
    fdummy = 1.0;
    for (int i = 0; i < nofNodeData; i++)
        out.write((char *)&fdummy, sizeof(float)); // max Wert pro Datentyp

    // daten ins file schreiben
    for (int d = 0; d < nofNodeData; ++d) {
        for (int n = 0; n < (int)nodedata[d].size(); n++) {
            fdummy = (float)nodedata[d][n];
            out.write((char *)&fdummy, sizeof(float));
        }
    }

    fdummy = 1.;
    for (int i = 0; i < nofNodeData; i++)
        out.write((char *)&fdummy, sizeof(float)); // max Wert pro Datentyp

    out.close();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeTrianglesWithNodeData to " << avsfilename << " - end");

    return avsfilename;
}
/*===============================================================================*/
std::string WbWriterAvsASCII::writeOctsWithCellData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                    vector<UbTupleInt8> &cells, vector<string> &datanames,
                                                    vector<vector<double>> &celldata)
{
    string avsfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeOctsWithCellData to " << avsfilename << " - start");

    ofstream out(avsfilename.c_str(), ios::out | ios::binary);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
        string path = UbSystem::getPathFromString(avsfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(avsfilename.c_str(), ios::out | ios::binary);
        }
        if (!out)
            throw UbException(UB_EXARGS, "file konnte nicht geschrieben werden " + avsfilename);
    }

    char magic = (char)7;
    int idummy;
    float fdummy;

    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    int nofNodeData     = 0;
    int nofCellData     = (int)datanames.size();
    int nofModelData    = 0;
    int cellType        = 7; //=hex
    int nofNodesPerCell = 8;

    out.write((char *)&magic, sizeof(char));
    out.write((char *)&nofNodes, sizeof(int));
    out.write((char *)&nofCells, sizeof(int));
    out.write((char *)&nofNodeData, sizeof(int));
    out.write((char *)&nofCellData, sizeof(int));
    out.write((char *)&nofModelData, sizeof(int));

    idummy = (int)nofCells * nofNodesPerCell;
    out.write((char *)&idummy, sizeof(int)); //(nof nodes) * (nodes per cell)
    for (int c = 0; c < nofCells; c++) {
        idummy = c + 1;
        out.write((char *)&idummy, sizeof(int)); // cell id
        idummy = 1;
        out.write((char *)&idummy, sizeof(int)); // mat
        idummy = nofNodesPerCell;
        out.write((char *)&idummy, sizeof(int)); // nodes per cell
        idummy = cellType;
        out.write((char *)&idummy, sizeof(int)); // cell type
    }
    // knotennummern der einzelnen zellen
    for (int c = 0; c < nofCells; c++) {
        idummy = val<1>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<2>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<3>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<4>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<5>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<6>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<7>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<8>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
    }

    // coords
    // x1-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<1>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // x2-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<2>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // x3-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<3>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // out<<"\n";

    // CELL DATA
    char labels[1024];
    char units[1024];
    strcpy(labels, "");
    strcpy(units, "");

    for (int d = 0; d < nofCellData - 1; ++d) {
        strcat(labels, datanames[d].c_str());
        strcat(labels, ".");
    }
    strcat(labels, datanames[nofCellData - 1].c_str());

    for (int d = 0; d < nofCellData - 1; ++d)
        strcat(units, "no_unit.");
    strcat(units, "no_unit");

    out.write((char *)&labels, sizeof(labels));
    out.write((char *)&units, sizeof(units));

    // nof and type of data
    idummy = nofCellData;
    out.write((char *)&idummy, sizeof(int)); // Datentypen pro knoten (hier = nof_node_data, da NUR skalare)

    idummy = 1;
    for (int i = 0; i < nofCellData; i++)
        out.write((char *)&idummy, sizeof(int)); // jeder Datentyp ist ein skalarer Wert

    // min and max of data
    fdummy = 0.0;
    for (int i = 0; i < nofCellData; i++)
        out.write((char *)&fdummy, sizeof(float)); // min Wert pro Datentyp
    fdummy = 1.0;
    for (int i = 0; i < nofCellData; i++)
        out.write((char *)&fdummy, sizeof(float)); // max Wert pro Datentyp

    // daten ins file schreiben
    for (int d = 0; d < nofCellData; ++d)
        for (int n = 0; n < (int)celldata[d].size(); n++) {
            fdummy = (float)celldata[d][n];
            out.write((char *)&fdummy, sizeof(float));
        }

    fdummy = 1.;
    for (int i = 0; i < nofCellData; i++)
        out.write((char *)&fdummy, sizeof(float)); // max Wert pro Datentyp

    out.close();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeOctsWithCellData to " << avsfilename << " - end");

    return avsfilename;
}
/*===============================================================================*/
std::string WbWriterAvsASCII::writeOctsWithNodeData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                    vector<UbTupleUInt8> &cells, vector<string> &datanames,
                                                    vector<vector<double>> &nodedata)
{
    string avsfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeOctsWithNodeData to " << avsfilename << " - start");

    ofstream out(avsfilename.c_str(), ios::out | ios::binary);
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
        string path = UbSystem::getPathFromString(avsfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(avsfilename.c_str(), ios::out | ios::binary);
        }
        if (!out)
            throw UbException(UB_EXARGS, "file konnte nicht geschrieben werden " + avsfilename);
    }

    if ((int)nodedata.size() == 0)
        throw UbException(UB_EXARGS, "no nodedata!!!");
    if (nodes.size() != nodedata[0].size())
        throw UbException(UB_EXARGS, "nodedata != nofNodes!!!");

    char magic = (char)7;
    int idummy;
    float fdummy;

    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    int nofNodeData     = (int)datanames.size();
    int nofCellData     = 0;
    int nofModelData    = 0;
    int cellType        = 7; //=hex
    int nofNodesPerCell = 8;

    out.write((char *)&magic, sizeof(char));
    out.write((char *)&nofNodes, sizeof(int));
    out.write((char *)&nofCells, sizeof(int));
    out.write((char *)&nofNodeData, sizeof(int));
    out.write((char *)&nofCellData, sizeof(int));
    out.write((char *)&nofModelData, sizeof(int));

    idummy = (int)nofCells * nofNodesPerCell;
    out.write((char *)&idummy, sizeof(int)); //(nof nodes) * (nodes per cell)
    for (int c = 0; c < nofCells; c++) {
        idummy = c + 1;
        out.write((char *)&idummy, sizeof(int)); // cell id
        idummy = 1;
        out.write((char *)&idummy, sizeof(int)); // mat
        idummy = nofNodesPerCell;
        out.write((char *)&idummy, sizeof(int)); // nodes per cell
        idummy = cellType;
        out.write((char *)&idummy, sizeof(int)); // cell type
    }
    // knotennummern der einzelnen zellen
    for (int c = 0; c < nofCells; c++) {
        idummy = val<1>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<2>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<3>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<4>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<5>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<6>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<7>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
        idummy = val<8>(cells[c]) + 1;
        out.write((char *)&idummy, sizeof(int));
    }

    // coords
    // x1-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<1>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // x2-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<2>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }
    // x3-coords
    for (int n = 0; n < nofNodes; n++) {
        fdummy = (float)(val<3>(nodes[n]));
        out.write((char *)&fdummy, sizeof(float));
    }

    // NODE DATA
    char labels[1024];
    char units[1024];
    strcpy(labels, "");
    strcpy(units, "");

    for (int d = 0; d < nofNodeData - 1; ++d) {
        strcat(labels, datanames[d].c_str());
        strcat(labels, ".");
    }
    strcat(labels, datanames[nofNodeData - 1].c_str());

    for (int i = 0; i < (nofNodeData - 1); i++)
        strcat(units, "no_unit.");
    strcat(units, "no_unit");

    out.write((char *)&labels, sizeof(labels));
    out.write((char *)&units, sizeof(units));

    // nof and type of data
    idummy = nofNodeData;
    out.write((char *)&idummy, sizeof(int)); // Datentypen pro knoten (hier = nof_node_data, da NUR skalare)

    idummy = 1;
    for (int i = 0; i < nofNodeData; i++)
        out.write((char *)&idummy, sizeof(int)); // jeder Datentyp ist ein skalarer Wert

    // min and max of data
    fdummy = 0.0;
    for (int i = 0; i < nofNodeData; i++)
        out.write((char *)&fdummy, sizeof(float)); // min Wert pro Datentyp
    fdummy = 1.0;
    for (int i = 0; i < nofNodeData; i++)
        out.write((char *)&fdummy, sizeof(float)); // max Wert pro Datentyp

    // daten ins file schreiben
    for (int d = 0; d < nofNodeData; ++d)
        for (int n = 0; n < (int)nodedata[d].size(); n++) {
            fdummy = (float)nodedata[d][n];
            out.write((char *)&fdummy, sizeof(float));
        }

    fdummy = 1.;
    for (int i = 0; i < nofNodeData; i++)
        out.write((char *)&fdummy, sizeof(float)); // max Wert pro Datentyp

    out.close();
    UBLOG(logDEBUG1, "WbWriterAvsASCII::writeOctsWithNodeData to " << avsfilename << " - end");

    return avsfilename;
}

//! \}
