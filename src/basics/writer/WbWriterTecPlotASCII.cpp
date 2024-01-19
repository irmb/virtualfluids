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
#include <basics/writer/WbWriterTecPlotASCII.h>

using namespace std;

/*===============================================================================*/
string WbWriterTecPlotASCII::writeOctsWithNodeData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                   vector<UbTupleUInt8> &cells, vector<string> &datanames,
                                                   vector<vector<double>> &nodedata)
{
    string tecplotfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterTecPlotASCII::writeOctsWithNodeData to " << tecplotfilename << " - start");

    ofstream out(tecplotfilename.c_str());
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(tecplotfilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(tecplotfilename.c_str());
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + tecplotfilename);
    }

    int nofNodes = (int)nodes.size();
    int nofCells = (int)cells.size();

    out << "TITLE = VirtualFluids OctGrid from " << UbSystem::getTimeStamp() << endl;

    out << "VARIABLES = \"X\", \"Y\", \"Z\"";
    for (size_t d = 0; d < datanames.size(); d++)
        out << ", \"" << datanames[d] << "\"";
    out << endl;

    out << "ZONE NODES=" << nofNodes << ", ELEMENTS=" << nofCells << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" << endl;
    for (size_t n = 0; n < nodes.size(); n++) {
        UbTupleFloat3 &coords = nodes[n];
        out << val<1>(coords) << " " << val<2>(coords) << " " << val<3>(coords);
        for (size_t d = 0; d < datanames.size(); d++)
            out << " " << nodedata[d][n];
        out << endl;
    }

    for (size_t c = 0; c < cells.size(); c++) {
        UbTupleUInt8 &cell = cells[c];
        out << val<1>(cell) << " " << val<2>(cell) << " " << val<3>(cell) << " " << val<4>(cell) << " " << val<5>(cell)
            << " " << val<6>(cell) << " " << val<7>(cell) << " " << val<8>(cell) << endl;
    }

    out.close();
    UBLOG(logDEBUG1, "WbWriterTecPlotASCII::writeOctsWithNodeData to " << tecplotfilename << " - end");

    return tecplotfilename;
}
/*===============================================================================*/
string WbWriterTecPlotASCII::writeOctsU(const string &filename, vector<UbTupleFloat3> &nodes,
                                        vector<UbTupleUInt8> &cells)
{
    vector<string> datanames;
    vector<vector<double>> nodedata;
    return writeOctsWithNodeData(filename, nodes, cells, datanames, nodedata);
}
/*===============================================================================*/

//! \}
