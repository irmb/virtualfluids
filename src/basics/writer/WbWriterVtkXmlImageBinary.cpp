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
//! \author Soeren Freudiger, Sebastian Geller, Henry Korb, Henrik Asmuth
//=======================================================================================
#include <basics/utilities/UbLogger.h>
#include <basics/writer/WbWriterVtkXmlImageBinary.h>
#include <cstring>

using namespace std;

/*===============================================================================*/
const std::string WbWriterVtkXmlImageBinary::pvdEndTag = "   </Collection>\n</VTKFile>";
/*===============================================================================*/
string WbWriterVtkXmlImageBinary::writeCollection(const string &filename, const vector<string> &filenames,
                                                  const double &timeStep, const bool &sepGroups)
{
    string vtkfilename = filename + ".pvd";
    ofstream out(vtkfilename.c_str());
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

    string endian;
    if (UbSystem::isLittleEndian())
        endian = "LittleEndian";
    else
        endian = "BigEndian";
    out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"" << endian << "\" >" << endl;
    out << "   <Collection>" << endl;

    int group = 0, part = 0;
    for (size_t i = 0; i < filenames.size(); i++) {
        out << "       <DataSet timestep=\"" << timeStep << "\" group=\"" << group << "\" part=\"" << part
            << "\" file=\"" << filenames[i] << "\"/>" << endl;
        if (sepGroups)
            group++;
        else
            part++;
    }
    out << pvdEndTag;
    out.close();

    return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlImageBinary::addFilesToCollection(const string &filename, const vector<string> &filenames,
                                                       const double &timeStep, const bool &sepGroups)
{
    string vtkfilename = filename;
    fstream test(vtkfilename.c_str(), ios::in);
    if (!test) {
        test.clear();
        vtkfilename += ".pvd";
        test.open(vtkfilename.c_str(), ios::in);
        if (!test)
            return this->writeCollection(filename, filenames, timeStep, sepGroups);
    }

    fstream out(vtkfilename.c_str(), ios::in | ios::out);
    out.seekp(-(int)pvdEndTag.size() - 1, ios_base::end);

    int group = 0;
    for (size_t i = 0; i < filenames.size(); i++) {
        out << "       <DataSet timestep=\"" << timeStep << "\" group=\"" << group << "\" part=\"" << i << "\" file=\""
            << filenames[i] << "\"/>" << endl;
        if (sepGroups)
            group++;
    }
    out << pvdEndTag;

    return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlImageBinary::writeParallelFile(const string &filename, const UbTupleInt6 &wholeExtent,
                                                    const UbTupleFloat3 &origin, const UbTupleFloat3 &spacing,
                                                    vector<string> &pieceSources, vector<UbTupleInt6> &pieceExtents,
                                                    vector<string> &pointDataNames, vector<string> &cellDataNames)
{
    string vtkfilename = filename + ".pvti";
    UBLOG(logDEBUG1, "WbWriterVtkXmlImageBinary::writeParallelFile to " << vtkfilename << " - start");

    ofstream out(vtkfilename.c_str());
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

    // VTK FILE
    out << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">"
        << "\n";
    out << "  <PImageData "
            << "WholeExtent=\"" << val<1>(wholeExtent) << " "
                                << val<2>(wholeExtent) << " " 
                                << val<3>(wholeExtent) << " " 
                                << val<4>(wholeExtent) << " " 
                                << val<5>(wholeExtent) << " "
                                << val<6>(wholeExtent) << "\" "
            << "GhostLevel=\"0\" "
            << "Origin=\""  << val<1>(origin) << " "
                            << val<2>(origin) << " "
                            << val<3>(origin) << "\" "
            << "Spacing=\"" << val<1>(spacing) << " "
                            << val<2>(spacing) << " "
                            << val<3>(spacing) << "\" "
        << "> \n";
    out << "    <PPointData>\n";
    for (size_t s = 0; s < pointDataNames.size(); s++)
        out << "      <PDataArray type=\"Float32\" Name=\"" << pointDataNames[s] << "\"/>\n";
    out << "    </PPointData>\n";
    if (cellDataNames.size() > 0) {
        out << "    <PCellData>\n";
        for (size_t s = 0; s < cellDataNames.size(); s++)
            out << "      <PDataArray type=\"Float32\" Name=\"" << cellDataNames[s] << "\"/>\n";
        out << "    </PCellData>\n";
    }
    for (size_t s = 0; s < pieceSources.size(); s++)
        out << "    <Piece Extent=\""   << val<1>(pieceExtents[s]) << " " 
                                        << val<2>(pieceExtents[s]) << " " 
                                        << val<3>(pieceExtents[s]) << " " 
                                        << val<4>(pieceExtents[s]) << " " 
                                        << val<5>(pieceExtents[s]) << " "
                                        << val<6>(pieceExtents[s]) << "\" Source=\"" << pieceSources[s] << "\"/>\n";
    out << "  </PImageData>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkXmlImageBinary::writeParallelFile to " << vtkfilename << " - end");

    return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlImageBinary::writeOctsWithCellData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                        vector<UbTupleInt8> & /*cells*/, vector<string> &datanames,
                                                        vector<vector<double>> &celldata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkXmlImageBinary::writeOctsWithCellData to " << vtkfilename << " - start");

    vector<string> nodeDataNames;
    vector<vector<double>> nodedata;

    UbTupleFloat3 origin, spacing;
    UbTupleInt6 extent;

    getMetaDataOfImage(nodes, origin, spacing, extent);

    this->writeData(vtkfilename, nodeDataNames, datanames, nodedata, celldata, extent, origin, spacing, extent);
    UBLOG(logDEBUG1, "WbWriterVtkXmlImageBinary::writeOctsWithCellData to " << vtkfilename << " - end");

    return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlImageBinary::writeOctsWithNodeData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                        vector<UbTupleUInt8> & /*cells*/, vector<string> &datanames,
                                                        vector<vector<double>> &nodedata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkXmlImageBinary::writeOctsWithNodeData to " << vtkfilename << " - start");

    vector<string> cellDataNames;
    vector<vector<double>> cellData;

    UbTupleFloat3 origin, spacing;
    UbTupleInt6 extent;

    getMetaDataOfImage(nodes, origin, spacing, extent);

    this->writeData(vtkfilename, datanames, cellDataNames, nodedata, cellData, extent, origin, spacing, extent);

    UBLOG(logDEBUG1, "WbWriterVtkXmlImageBinary::writeOctsWithNodeData to " << vtkfilename << " - end");

    return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlImageBinary::writeNodesWithNodeData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                         vector<string> &datanames, vector<vector<double>> &nodedata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkXmlImageBinary::writeNodesWithNodeData to " << vtkfilename << " - start");

    vector<string> cellDataNames;
    vector<vector<double>> cellData;

    UbTupleFloat3 origin, spacing;
    UbTupleInt6 extent;

    getMetaDataOfImage(nodes, origin, spacing, extent);
    this->writeData(vtkfilename, datanames, cellDataNames, nodedata, cellData, extent, origin, spacing, extent);

    return vtkfilename;
}

void WbWriterVtkXmlImageBinary::getMetaDataOfImage(vector<UbTupleFloat3> &nodes, UbTupleFloat3 &origin,
                                                   UbTupleFloat3 &spacing, UbTupleInt6 &extent)
{
    int nofNodes = (int)nodes.size();
    val<1>(origin) = val<1>(nodes[0]);
    val<2>(origin) = val<2>(nodes[0]);
    val<3>(origin) = val<3>(nodes[0]);

    float l_x = val<1>(nodes[nofNodes-1])-val<1>(origin);
    float l_y = val<2>(nodes[nofNodes-1])-val<2>(origin);

    val<1>(spacing) = val<1>(nodes[1])-val<1>(nodes[0]);
    int nx = (l_x) / val<1>(spacing);
    val<2>(spacing) = val<2>(nodes[nx])-val<2>(nodes[0]);    
    int ny = (l_y) / val<2>(spacing);
    val<3>(spacing) = val<3>(nodes[nx*ny])-val<3>(nodes[0]);

    val<1>(extent) = val<1>(origin) / val<1>(spacing); val<2>(extent) = val<1>(nodes[nofNodes - 1]) / val<1>(spacing);    
    val<3>(extent) = val<2>(origin) / val<2>(spacing); val<4>(extent) = val<2>(nodes[nofNodes - 1]) / val<2>(spacing);    
    val<5>(extent) = val<3>(origin) / val<3>(spacing); val<6>(extent) = val<3>(nodes[nofNodes - 1]) / val<3>(spacing);    

}

void WbWriterVtkXmlImageBinary::writeData(const string &vtkfilename, vector<string> &pointDataNames,
                                          vector<string> &cellDataNames, vector<vector<double>> &nodedata,
                                          vector<vector<double>> &celldata, UbTupleInt6 &wholeExtent,
                                          UbTupleFloat3 &origin, UbTupleFloat3 &spacing, UbTupleInt6 &extent,
                                          unsigned int precision)
{
    ofstream out(vtkfilename.c_str(), ios::out | ios::binary);
    out.precision(precision);

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

    size_t nPoints = pointDataNames.size() > 0 ? nodedata[0].size() : celldata[0].size();

    int bytesPerByteVal = 4; //==sizeof(int)

    int bytesScalarData = 1 /*scalar         */ * (int)nPoints * sizeof(double);

    int offset = 0;

    // VTK FILE
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" >"
        << "\n";
    out << "   <ImageData "
            << "WholeExtent=\"" << val<1>(wholeExtent) << " " 
                                << val<2>(wholeExtent) << " " 
                                << val<3>(wholeExtent) << " " 
                                << val<4>(wholeExtent) << " " 
                                << val<5>(wholeExtent) << " "
                                << val<6>(wholeExtent) << "\" "
            << "Origin=\""  << val<1>(origin) << " " 
                            << val<2>(origin) << " "
                            << val<3>(origin) << "\" "
            << "Spacing=\"" << val<1>(spacing) << " " 
                            << val<2>(spacing) << " " 
                            << val<3>(spacing) << "\""
        << "> \n";
    out << "      <Piece Extent=\"" << val<1>(extent) << " " 
                                    << val<2>(extent) << " " 
                                    << val<3>(extent) << " " 
                                    << val<4>(extent) << " " 
                                    << val<5>(extent) << " "
                                    << val<6>(extent) << "\">\n";

    // DATA SECTION
    if (pointDataNames.size() > 0) {
        out << "         <PointData>\n";
        for (size_t s = 0; s < pointDataNames.size(); ++s) {
            out << "            <DataArray type=\"Float64\" Name=\"" << pointDataNames[s]
                << "\" format=\"appended\" offset=\"" << offset << "\" /> \n";
            offset += (bytesPerByteVal + bytesScalarData);
        }
        out << "         </PointData>\n";
    }

    if (cellDataNames.size() > 0) {
        out << "         <CellData>\n";
        for (size_t s = 0; s < cellDataNames.size(); ++s) {
            out << "            <DataArray type=\"Float64\" Name=\"" << cellDataNames[s]
                << "\" format=\"appended\" offset=\"" << offset << "\" /> \n";
            offset += (bytesPerByteVal + bytesScalarData);
        }
        out << "         </CellData>\n";
    }

    out << "      </Piece>\n";
    out << "   </ImageData>\n";

    // AppendedData SECTION
    out << "   <AppendedData encoding=\"raw\">\n";
    out << "_";

    // DATA SECTION
    // pointData
    for (size_t s = 0; s < pointDataNames.size(); ++s) {
        out.write((char *)&bytesScalarData, bytesPerByteVal);
        for (size_t d = 0; d < nodedata[s].size(); ++d) {
            double tmp = nodedata[s][d];
            out.write((char *)&tmp, sizeof(double));
        }
    }

    // cellData
    for (size_t s = 0; s < cellDataNames.size(); ++s) {
        out.write((char *)&bytesScalarData, bytesPerByteVal);
        for (size_t d = 0; d < celldata[s].size(); ++d) {
            double tmp = celldata[s][d];
            out.write((char *)&tmp, sizeof(double));
        }
    }
    out << "\n   </AppendedData>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
}

//! \}
