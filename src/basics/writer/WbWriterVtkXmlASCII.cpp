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
#include <basics/writer/WbWriterVtkXmlASCII.h>
#include <cstring>
#include <limits>

using namespace std;

/*===============================================================================*/
std::string WbWriterVtkXmlASCII::pvdEndTag = "   </Collection>\n</VTKFile>";
/*===============================================================================*/
std::string WbWriterVtkXmlASCII::writeCollection(const std::string &filename, const std::vector<std::string> &filenames,
                                                 const double &timeStep, const bool &sepGroups)
{
    std::string vtkfilename = filename + ".pvd";
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

    std::string endian;
    if (UbSystem::isLittleEndian())
        endian = "LittleEndian";
    else
        endian = "BigEndian";
    out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"" << endian << "\" >\n";
    out << "   <Collection>" << endl;

    int group = 0, part = 0;
    for (std::size_t i = 0; i < filenames.size(); i++) {
        out << "       <DataSet timestep=\"" << timeStep << "\" group=\"" << group << "\" part=\"" << part
            << "\" file=\"" << filenames[i] << "\"/>\n";
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
std::string WbWriterVtkXmlASCII::addFilesToCollection(const std::string &filename,
                                                      const std::vector<std::string> &filenames, const double &timeStep,
                                                      const bool &sepGroups)
{
    std::string vtkfilename = filename;
    std::fstream test(vtkfilename.c_str(), ios::in);
    if (!test) {
        test.clear();
        vtkfilename += ".pvd";
        test.open(vtkfilename.c_str(), ios::in);
        if (!test)
            return this->writeCollection(filename, filenames, timeStep, sepGroups);
    }

    std::fstream out(vtkfilename.c_str(), ios::in | ios::out);
    out.seekp(-(int)pvdEndTag.size() - 1, ios_base::end);

    int group = 0;
    for (std::size_t i = 0; i < filenames.size(); i++) {
        out << "       <DataSet timestep=\"" << timeStep << "\" group=\"" << group << "\" part=\"" << i << "\" file=\""
            << filenames[i] << "\"/>\n";
        if (sepGroups)
            group++;
    }
    out << pvdEndTag;

    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkXmlASCII::writeParallelFile(const string &filename, vector<string> &pieceSources,
                                                   vector<string> &pointDataNames, vector<string> &cellDataNames)
{
    string vtkfilename = filename + ".pvtu";
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeParallelFile to " << vtkfilename << " - start");

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

    // VTK FILE
    out << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
    out << "    <PPoints>\n";
    out << "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
    out << "    </PPoints>\n";
    out << "    <PPointData>\n";
    for (size_t s = 0; s < pointDataNames.size(); s++)
        out << "      <PDataArray type=\"Float64\" Name=\"" << pointDataNames[s] << "\"/>\n";
    out << "    </PPointData>\n";
    if (cellDataNames.size() > 0) {
        out << "    <PCellData>\n";
        for (size_t s = 0; s < cellDataNames.size(); s++)
            out << "      <PDataArray type=\"Float32\" Name=\"" << cellDataNames[s] << "\"/>\n";
        out << "    </PCellData>\n";
    }

    for (size_t s = 0; s < pieceSources.size(); s++)
        out << "    <Piece Source=\"" << pieceSources[s] << "\"/>\n";
    out << "  </PUnstructuredGrid>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeParallelFile to " << vtkfilename << " - end");

    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkXmlASCII::writeQuads(const string &filename, vector<UbTupleFloat3> &nodes,
                                            vector<UbTupleInt4> &cells)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeQuads to " << vtkfilename << " - start");

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

    // VTK FILE
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"
        << "\n";
    out << "   <UnstructuredGrid>"
        << "\n";
    out << "      <Piece NumberOfPoints=\"" << nofNodes << "\"    NumberOfCells=\"" << nofCells << "\">   \n";

    // POINTS SECTION
    out << "      <Points>\n";
    out << "         <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << "   ";

    out << "\n";
    out << "         </DataArray>\n";
    out << "      </Points>\n";

    // CELLS SECTION
    out << "      <Cells>\n";
    out << "         <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

    for (int c = 0; c < nofCells; c++)
        out << val<1>(cells[c]) << " " << val<2>(cells[c]) << " " << val<4>(cells[c]) << " " << val<3>(cells[c])
            << "   ";
    out << "\n";
    out << "      </DataArray>\n";
    out << "         <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int c = 1; c < nofCells + 1; c++)
        out << c * 4 << " ";

    out << "\n";
    out << "         </DataArray>\n";

    out << "      <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";

    for (int c = 0; c < nofCells; c++)
        out << "8 ";
    out << "\n";
    out << "      </DataArray>\n";
    out << "      </Cells>\n";
    out << "      </Piece>\n";
    out << "   </UnstructuredGrid>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeQuads to " << vtkfilename << " - end");

    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkXmlASCII::writeQuadsWithNodeData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                        vector<UbTupleInt4> &cells, vector<string> &datanames,
                                                        vector<vector<double>> &nodedata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeQuadsWithNodeData to " << vtkfilename << " - start");

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

    // VTK FILE
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"
        << "\n";
    out << "   <UnstructuredGrid>"
        << "\n";
    out << "      <Piece NumberOfPoints=\"" << nofNodes << "\"    NumberOfCells=\"" << nofCells << "\">   \n";

    // POINTS SECTION
    out << "         <Points>\n";
    out << "            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n               ";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << "   ";

    out << "\n";
    out << "            </DataArray>\n";
    out << "         </Points>\n";

    // CELLS SECTION
    out << "         <Cells>\n";
    out << "            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n              ";

    for (int c = 0; c < nofCells; c++)
        out << val<1>(cells[c]) << " " << val<2>(cells[c]) << " " << val<4>(cells[c]) << " " << val<3>(cells[c])
            << "   ";
    out << "\n";
    out << "            </DataArray>\n";
    out << "            <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n              ";
    for (int c = 1; c < nofCells + 1; c++)
        out << c * 4 << " ";

    out << "\n";
    out << "            </DataArray>\n";

    out << "            <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n              ";

    for (int c = 0; c < nofCells; c++)
        out << "8 ";
    out << "\n";
    out << "            </DataArray>\n";

    out << "         </Cells>\n";

    // write data section
    out << "         <PointData Scalars=\"Scalars\"> \n";
    for (int s = 0; s < (int)datanames.size(); ++s) {
        out << "           <DataArray type=\"Float32\" Name=\"" << datanames[s] << "\" format=\"ascii\"> \n";

        for (int d = 0; d < (int)nodedata[s].size(); d++)
            out << nodedata[s][d] << " ";

        out << "\n          </DataArray>\n";
    }
    out << "         </PointData>\n";
    out << "      </Piece>\n";
    out << "   </UnstructuredGrid>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeQuadsWithNodeData to " << vtkfilename << " - end");
    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkXmlASCII::writeQuadsWithCellData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                        vector<UbTupleInt4> &cells, vector<string> &datanames,
                                                        vector<vector<double>> &celldata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeQuadsWithCellData to " << vtkfilename << " - start");

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

    // VTK FILE
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"
        << "\n";
    out << "   <UnstructuredGrid>"
        << "\n";
    out << "      <Piece NumberOfPoints=\"" << nofNodes << "\"    NumberOfCells=\"" << nofCells << "\">   \n";

    // POINTS SECTION
    out << "         <Points>\n";
    out << "            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n               ";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << "   ";

    out << "\n";
    out << "            </DataArray>\n";
    out << "         </Points>\n";

    // CELLS SECTION
    out << "         <Cells>\n";
    out << "            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n              ";

    for (int c = 0; c < nofCells; c++)
        out << val<1>(cells[c]) << " " << val<2>(cells[c]) << " " << val<4>(cells[c]) << " " << val<3>(cells[c])
            << "   ";
    out << "\n";
    out << "            </DataArray>\n";
    out << "            <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n              ";
    for (int c = 1; c < nofCells + 1; c++)
        out << c * 4 << " ";

    out << "\n";
    out << "            </DataArray>\n";

    out << "            <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n              ";

    for (int c = 0; c < nofCells; c++)
        out << "8 ";
    out << "\n";
    out << "            </DataArray>\n";

    out << "         </Cells>\n";

    // write data section
    out << "         <CellData Scalars=\"Scalars\"> \n";
    for (int s = 0; s < (int)datanames.size(); ++s) {
        out << "           <DataArray type=\"Float32\" Name=\"" << datanames[s] << "\" format=\"ascii\"> \n";

        for (int d = 0; d < (int)celldata[s].size(); d++)
            out << celldata[s][d] << " ";

        out << "\n          </DataArray>\n";
    }
    out << "         </CellData>\n";
    out << "      </Piece>\n";
    out << "   </UnstructuredGrid>\n";
    out << "</VTKFile>";
    out << endl;

    out.close();

    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeQuadsWithCellData to " << vtkfilename << " - end");

    return vtkfilename;
}
/*===============================================================================*/
string WbWriterVtkXmlASCII::writeQuadsWithNodeAndCellData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                          vector<UbTupleInt4> &cells, vector<string> &nodedatanames,
                                                          vector<vector<double>> &nodedata,
                                                          vector<string> &celldatanames,
                                                          vector<vector<double>> &celldata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeQuadsWithNodeAndCellData to " << vtkfilename << " - start");

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

    // VTK FILE
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"
        << "\n";
    out << "   <UnstructuredGrid>"
        << "\n";
    out << "      <Piece NumberOfPoints=\"" << nofNodes << "\"    NumberOfCells=\"" << nofCells << "\">   \n";

    // POINTS SECTION
    out << "         <Points>\n";
    out << "            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n               ";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << "   ";

    out << "\n";
    out << "            </DataArray>\n";
    out << "         </Points>\n";

    // CELLS SECTION
    out << "         <Cells>\n";
    out << "            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n              ";

    for (int c = 0; c < nofCells; c++)
        out << val<1>(cells[c]) << " " << val<2>(cells[c]) << " " << val<4>(cells[c]) << " " << val<3>(cells[c])
            << "   ";
    out << "\n";
    out << "            </DataArray>\n";
    out << "            <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n              ";
    for (int c = 1; c < nofCells + 1; c++)
        out << c * 4 << " ";

    out << "\n";
    out << "            </DataArray>\n";

    out << "            <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n              ";

    for (int c = 0; c < nofCells; c++)
        out << "8 ";
    out << "\n";
    out << "            </DataArray>\n";

    out << "         </Cells>\n";

    // write PointData section
    out << "         <PointData Scalars=\"PScalars\"> \n";
    for (int s = 0; s < (int)nodedatanames.size(); ++s) {
        out << "           <DataArray type=\"Float32\" Name=\"" << nodedatanames[s] << "\" format=\"ascii\"> \n";

        for (int d = 0; d < (int)nodedata[s].size(); d++)
            out << nodedata[s][d] << " ";

        out << "\n          </DataArray>\n";
    }
    out << "         </PointData>\n";

    // write celldata section
    out << "         <CellData Scalars=\"CScalars\"> \n";
    for (int s = 0; s < (int)celldatanames.size(); ++s) {
        out << "           <DataArray type=\"Float32\" Name=\"" << celldatanames[s] << "\" format=\"ascii\"> \n";

        for (int d = 0; d < (int)celldata[s].size(); d++)
            out << celldata[s][d] << " ";

        out << "\n          </DataArray>\n";
    }
    out << "         </CellData>\n";
    out << "      </Piece>\n";
    out << "   </UnstructuredGrid>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeQuadsWithNodeAndCellData to " << vtkfilename << " - end");

    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkXmlASCII::writeLines(const string &filename, vector<UbTupleFloat3> &nodes,
                                            vector<UbTupleInt2> &lines)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeLines to " << vtkfilename << " - start");

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

    // VTK FILE
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"
        << "\n";
    out << "   <UnstructuredGrid>"
        << "\n";
    out << "      <Piece NumberOfPoints=\"" << nofNodes << "\"    NumberOfCells=\"" << nofLines << "\">   \n";

    // POINTS SECTION
    out << "      <Points>\n";
    out << "         <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << "   ";

    out << "\n";
    out << "         </DataArray>\n";
    out << "      </Points>\n";

    // CELLS SECTION
    out << "      <Cells>\n";
    out << "         <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

    for (int c = 0; c < nofLines; c++)
        out << val<1>(lines[c]) << " " << val<2>(lines[c]) << "  ";
    out << "\n";
    out << "      </DataArray>\n";
    out << "         <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int c = 1; c <= nofLines; c++)
        out << c * 2 << " ";

    out << "\n";
    out << "         </DataArray>\n";

    out << "      <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";

    for (int c = 0; c < nofLines; c++)
        out << "3 ";
    out << "\n";
    out << "      </DataArray>\n";
    out << "      </Cells>\n";
    out << "      </Piece>\n";
    out << "   </UnstructuredGrid>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeLines to " << vtkfilename << " - end");

    return vtkfilename;
}
/*===============================================================================*/
// std::string WbWriterVtkXmlASCII::writeLinesWithNodeData(const string& filename,vector<UbTupleFloat3 >& nodes,
// vector<UbTupleInt2 >& lines, std::vector< std::string >& datanames, std::vector< std::vector< double > >& nodedata)
//{
//   string vtkfilename=filename+getFileExtension();
//   UBLOG(logDEBUG1,"WbWriterVtkXmlASCII::writeLinesWithNodeData to "<<vtkfilename<<" - start");
//
//   std::ofstream out(vtkfilename.c_str());
//   if(!out)
//   {
//      out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
//      string path = UbSystem::getPathFromString(vtkfilename);
//      if(path.size()>0){UbSystem::makeDirectory(path);out.open(vtkfilename.c_str());}
//      if(!out) throw UbException(UB_EXARGS,"couldn't open file "+vtkfilename);
//   }
//
//   int nofNodes = (int)nodes.size();
//   int nofLines = (int)lines.size();
//
//   //VTK FILE
//   out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"<<"\n";
//   out<<"   <UnstructuredGrid>"<<"\n";
//   out<<"      <Piece NumberOfPoints=\""<<nofNodes<<"\"    NumberOfCells=\""<<nofLines<<"\">   \n";
//
//   //POINTS SECTION
//   out<<"      <Points>\n";
//   out<<"         <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
//   for(int n=0; n<nofNodes; n++)
//      out<< val<1>(nodes[n]) <<" "<< val<2>(nodes[n]) <<" "<< val<3>(nodes[n]) <<"   ";
//
//   out<<"\n";
//   out<<"         </DataArray>\n";
//   out<<"      </Points>\n";
//
//   //CELLS SECTION
//   out<<"      <Cells>\n";
//   out<<"         <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
//
//   for(int c=0; c<nofLines; c++)
//      out<< val<1>(lines[c]) <<" "<< val<2>(lines[c])<<"  ";
//   out<<"\n";
//   out<<"      </DataArray>\n";
//   out<<"         <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
//   for(int c=1; c<=nofLines; c++)
//      out<<c*2<<" " ;
//
//   out<<"\n";
//   out<<"         </DataArray>\n";
//
//   out<<"      <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
//
//   for(int c=0; c<nofLines; c++)
//      out<<"3 ";
//   out<<"\n";
//   out<<"      </DataArray>\n";
//   out<<"      </Cells>\n";
//
//   //write data section
//   out<<"         <PointData Scalars=\"Scalars\"> \n";
//   for(int s=0; s<(int)datanames.size(); ++s)
//   {
//      out<< "           <DataArray type=\"Float32\" Name=\""<< datanames[s] <<"\" format=\"ascii\"> \n";
//
//      for(int d=0; d<(int)nodedata[s].size(); d++)
//         out<<nodedata[s][d]<<" ";
//
//      out<<"\n          </DataArray>\n";
//   }
//   out<<"         </PointData>\n";
//   out<<"      </Piece>\n";
//   out<<"   </UnstructuredGrid>\n";
//   out<<"</VTKFile>";
//   out<<endl;
//   out.close();
//   UBLOG(logDEBUG1,"WbWriterVtkXmlASCII::writeLinesWithNodeData to "<<vtkfilename<<" - end");
//
//   return vtkfilename;
//}
/*===============================================================================*/
std::string WbWriterVtkXmlASCII::writeTriangles(const string &filename, vector<UbTupleFloat3> &nodes,
                                                vector<UbTupleInt3> &triangles)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeTriangles to " << vtkfilename << " - start");

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

    // VTK FILE
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"
        << "\n";
    out << "   <UnstructuredGrid>"
        << "\n";
    out << "      <Piece NumberOfPoints=\"" << nofNodes << "\"    NumberOfCells=\"" << nofTriangles << "\">   \n";

    // POINTS SECTION
    out << "      <Points>\n";
    out << "         <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << "   ";

    out << "\n";
    out << "         </DataArray>\n";
    out << "      </Points>\n";

    // CELLS SECTION
    out << "      <Cells>\n";
    out << "         <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

    for (int c = 0; c < nofTriangles; c++)
        out << val<1>(triangles[c]) << " " << val<2>(triangles[c]) << " " << val<3>(triangles[c]) << "  ";
    out << "\n";
    out << "      </DataArray>\n";
    out << "         <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int c = 1; c < nofTriangles + 1; c++)
        out << c * 3 << " ";

    out << "\n";
    out << "         </DataArray>\n";

    out << "      <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";

    for (int c = 0; c < nofTriangles; c++)
        out << "5 ";
    out << "\n";
    out << "      </DataArray>\n";
    out << "      </Cells>\n";
    out << "      </Piece>\n";
    out << "   </UnstructuredGrid>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeTriangles to " << vtkfilename << " - end");

    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkXmlASCII::writeTrianglesWithNodeData(const std::string &filename,
                                                            std::vector<UbTupleFloat3> &nodes,
                                                            std::vector<UbTupleInt3> &cells,
                                                            std::vector<std::string> &datanames,
                                                            std::vector<std::vector<double>> &nodedata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeTrianglesWithNodeData to " << vtkfilename << " - start");

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

    // VTK FILE
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"
        << "\n";
    out << "   <UnstructuredGrid>"
        << "\n";
    out << "      <Piece NumberOfPoints=\"" << nofNodes << "\"    NumberOfCells=\"" << nofCells << "\">   \n";

    // POINTS SECTION
    out << "         <Points>\n";
    out << "            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n               ";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << "   ";

    out << "\n";
    out << "            </DataArray>\n";
    out << "         </Points>\n";

    // CELLS SECTION
    out << "         <Cells>\n";
    out << "            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n              ";

    for (int c = 0; c < nofCells; c++)
        out << val<1>(cells[c]) << " " << val<2>(cells[c]) << " " << val<3>(cells[c]) << "   ";
    out << "\n";
    out << "            </DataArray>\n";
    out << "            <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n              ";
    for (int c = 1; c < nofCells + 1; c++)
        out << c * 3 << " ";

    out << "\n";
    out << "            </DataArray>\n";

    out << "            <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n              ";

    for (int c = 0; c < nofCells; c++)
        out << "5 ";
    out << "\n";
    out << "            </DataArray>\n";

    out << "         </Cells>\n";

    // write data section
    out << "         <PointData Scalars=\"Scalars\"> \n";
    for (int s = 0; s < (int)datanames.size(); ++s) {
        out << "           <DataArray type=\"Float32\" Name=\"" << datanames[s] << "\" format=\"ascii\"> \n";

        for (int d = 0; d < (int)nodedata[s].size(); d++)
            out << nodedata[s][d] << " ";

        out << "\n          </DataArray>\n";
    }
    out << "         </PointData>\n";
    out << "      </Piece>\n";
    out << "   </UnstructuredGrid>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeTrianglesWithNodeData to " << vtkfilename << " - end");

    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkXmlASCII::writeOctsWithCellData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                       vector<UbTupleInt8> &cells, vector<string> &datanames,
                                                       vector<vector<double>> &celldata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeOctsWithCellData to " << vtkfilename << " - start");

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

    // VTK FILE
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"
        << "\n";
    out << "   <UnstructuredGrid>"
        << "\n";
    out << "      <Piece NumberOfPoints=\"" << nofNodes << "\"    NumberOfCells=\"" << nofCells << "\">   \n";

    // POINTS SECTION
    out << "         <Points>\n";
    out << "            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n               ";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << "   ";

    out << "\n";
    out << "            </DataArray>\n";
    out << "         </Points>\n";

    // CELLS SECTION
    out << "         <Cells>\n";
    out << "            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n              ";

    for (int c = 0; c < nofCells; c++)
        out << val<1>(cells[c]) << " " << val<2>(cells[c]) << " " << val<4>(cells[c]) << " " << val<3>(cells[c]) << " "
            << val<5>(cells[c]) << " " << val<6>(cells[c]) << " " << val<8>(cells[c]) << " " << val<7>(cells[c])
            << "  ";
    out << "\n";
    out << "            </DataArray>\n";
    out << "            <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n              ";
    for (int c = 1; c < nofCells + 1; c++)
        out << c * 8 << " ";

    out << "\n";
    out << "            </DataArray>\n";

    out << "            <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n              ";

    for (int c = 0; c < nofCells; c++)
        out << "11 ";
    out << "\n";
    out << "            </DataArray>\n";

    out << "         </Cells>\n";

    // write data section
    out << "         <CellData Scalars=\"Scalars\"> \n";
    for (int s = 0; s < (int)datanames.size(); ++s) {
        out << "           <DataArray type=\"Float32\" Name=\"" << datanames[s] << "\" format=\"ascii\"> \n";

        for (int d = 0; d < (int)celldata[s].size(); d++)
            out << celldata[s][d] << " ";

        out << "\n          </DataArray>\n";
    }
    out << "         </CellData>\n";

    out << "      </Piece>\n";
    out << "   </UnstructuredGrid>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeOctsWithCellData to " << vtkfilename << " - end");

    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkXmlASCII::writeOctsWithNodeData(const string &filename, vector<UbTupleFloat3> &nodes,
                                                       vector<UbTupleUInt8> &cells, vector<string> &datanames,
                                                       vector<vector<double>> &nodedata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeOctsWithNodeData to " << vtkfilename << " - start");

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

    // VTK FILE
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"
        << "\n";
    out << "   <UnstructuredGrid>"
        << "\n";
    out << "      <Piece NumberOfPoints=\"" << nofNodes << "\"    NumberOfCells=\"" << nofCells << "\">   \n";

    // POINTS SECTION
    out << "         <Points>\n";
    out << "            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n               ";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << " ";

    out << "\n";
    out << "            </DataArray>\n";
    out << "         </Points>\n";

    // CELLS SECTION
    out << "         <Cells>\n";
    out << "            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n              ";

    for (int c = 0; c < nofCells; c++)
        out << val<1>(cells[c]) << " " << val<2>(cells[c]) << " " << val<4>(cells[c]) << " " << val<3>(cells[c]) << " "
            << val<5>(cells[c]) << " " << val<6>(cells[c]) << " " << val<8>(cells[c]) << " " << val<7>(cells[c])
            << "  ";

    out << "\n";
    out << "            </DataArray>\n";
    out << "            <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n              ";
    for (int c = 1; c < nofCells + 1; c++)
        out << c * 8 << " ";

    out << "\n";
    out << "            </DataArray>\n";

    out << "            <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n              ";

    for (int c = 0; c < nofCells; c++)
        out << "11 ";
    out << "\n";
    out << "            </DataArray>\n";

    out << "         </Cells>\n";

    // write PointData section
    out << "         <PointData Scalars=\"PScalars\"> \n";
    for (int s = 0; s < (int)datanames.size(); ++s) {
        out << "           <DataArray type=\"Float64\" Name=\"" << datanames[s] << "\" format=\"ascii\">";

        for (int d = 0; d < (int)nodedata[s].size(); d++) {
            // out<<base64_encode((unsigned char*)(&nodedata[s][d]),sizeof(float));
            // out.write((char*)&nodedata[s][d],sizeof(float));
            out << nodedata[s][d] << " ";
        }
        out << "</DataArray>\n";
    }
    out << "         </PointData>\n";
    out << "      </Piece>\n";
    out << "   </UnstructuredGrid>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeOctsWithNodeData to " << vtkfilename << " - end");

    return vtkfilename;
}
/*===============================================================================*/
std::string WbWriterVtkXmlASCII::writeOcts(const string &filename, vector<UbTupleFloat3> &nodes,
                                           vector<UbTupleInt8> &cells)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeOcts to " << vtkfilename << " - start");

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

    // VTK FILE
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"
        << "\n";
    out << "   <UnstructuredGrid>"
        << "\n";
    out << "      <Piece NumberOfPoints=\"" << nofNodes << "\"    NumberOfCells=\"" << nofCells << "\">   \n";

    // POINTS SECTION
    out << "         <Points>\n";
    out << "            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n               ";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << "   ";

    out << "\n";
    out << "            </DataArray>\n";
    out << "         </Points>\n";

    // CELLS SECTION
    out << "         <Cells>\n";
    out << "            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n              ";

    for (int c = 0; c < nofCells; c++)
        out << val<1>(cells[c]) << " " << val<2>(cells[c]) << " " << val<4>(cells[c]) << " " << val<3>(cells[c]) << " "
            << val<5>(cells[c]) << " " << val<6>(cells[c]) << " " << val<8>(cells[c]) << " " << val<7>(cells[c])
            << "   ";
    out << "\n";
    out << "            </DataArray>\n";
    out << "            <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n              ";
    for (int c = 1; c < nofCells + 1; c++)
        out << c * 8 << " ";

    out << "\n";
    out << "            </DataArray>\n";

    out << "            <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n              ";

    for (int c = 0; c < nofCells; c++)
        out << "11 ";
    out << "\n";
    out << "            </DataArray>\n";
    out << "         </Cells>\n";
    out << "      </Piece>\n";
    out << "   </UnstructuredGrid>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeOcts to " << vtkfilename << " - end");

    return vtkfilename;
}
std::string WbWriterVtkXmlASCII::writeNodes(const std::string &filename, std::vector<UbTupleFloat3> &nodes)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeLines to " << vtkfilename << " - start");

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

    // VTK FILE
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"
        << "\n";
    out << "   <UnstructuredGrid>"
        << "\n";
    out << "      <Piece NumberOfPoints=\"" << nofNodes << "\"    NumberOfCells=\"" << nofNodes << "\">   \n";

    // POINTS SECTION
    out << "      <Points>\n";
    out << "         <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << "   ";

    out << "\n";
    out << "         </DataArray>\n";
    out << "      </Points>\n";

    // CELLS SECTION
    out << "      <Cells>\n";
    out << "         <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int n = 0; n < nofNodes; n++)
        out << n << "  ";
    out << "\n";

    out << "      </DataArray>\n";
    out << "         <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int n = 1; n <= nofNodes; n++)
        out << n << " ";

    out << "\n";
    out << "         </DataArray>\n";

    out << "      <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";

    for (int n = 0; n < nofNodes; n++)
        out << "1 ";
    out << "\n";
    out << "      </DataArray>\n";
    out << "      </Cells>\n";
    out << "      </Piece>\n";
    out << "   </UnstructuredGrid>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeLines to " << vtkfilename << " - end");

    return vtkfilename;
}
std::string WbWriterVtkXmlASCII::writeNodesWithNodeData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                                        std::vector<std::string> &datanames,
                                                        std::vector<std::vector<double>> &nodedata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeNodesWithNodeData to " << vtkfilename << " - start");

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

    // VTK FILE
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"
        << "\n";
    out << "   <UnstructuredGrid>"
        << "\n";
    out << "      <Piece NumberOfPoints=\"" << nofNodes << "\"    NumberOfCells=\"" << nofNodes << "\">   \n";

    // POINTS SECTION
    out << "         <Points>\n";
    out << "            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n               ";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << "   ";

    out << "\n";
    out << "            </DataArray>\n";
    out << "         </Points>\n";

    // CELLS SECTION
    out << "         <Cells>\n";
    out << "            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n              ";

    for (int c = 0; c < nofNodes; c++)
        out << c << "   ";
    out << "\n";

    out << "            </DataArray>\n";
    out << "            <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n              ";
    for (int c = 1; c < nofNodes + 1; c++)
        out << c << " ";

    out << "\n";
    out << "            </DataArray>\n";

    out << "            <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n              ";
    for (int c = 0; c < nofNodes; c++)
        out << "1 ";

    out << "\n";
    out << "            </DataArray>\n";

    out << "         </Cells>\n";

    // write data section
    out << "         <PointData Scalars=\"Scalars\"> \n";
    for (int s = 0; s < (int)datanames.size(); ++s) {
        out << "           <DataArray type=\"Float32\" Name=\"" << datanames[s] << "\" format=\"ascii\"> \n";

        for (int d = 0; d < (int)nodedata[s].size(); d++)
            out << nodedata[s][d] << " ";

        out << "\n          </DataArray>\n";
    }
    out << "         </PointData>\n";
    out << "      </Piece>\n";
    out << "   </UnstructuredGrid>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeNodesWithNodeData to " << vtkfilename << " - end");

    return vtkfilename;
}

//////////////////////////////////////////////////////////////////////////
std::string WbWriterVtkXmlASCII::writeNodesWithNodeDataDouble(const std::string &filename,
                                                              std::vector<UbTupleDouble3> &nodes,
                                                              std::vector<std::string> &datanames,
                                                              std::vector<std::vector<double>> &nodedata)
{
    string vtkfilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeNodesWithNodeData to " << vtkfilename << " - start");

    std::ofstream out(vtkfilename.c_str());
    out.precision(std::numeric_limits<double>::digits10 + 1);
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

    // VTK FILE
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >"
        << "\n";
    out << "   <UnstructuredGrid>"
        << "\n";
    out << "      <Piece NumberOfPoints=\"" << nofNodes << "\"    NumberOfCells=\"" << nofNodes << "\">   \n";

    // POINTS SECTION
    out << "         <Points>\n";
    out << "            <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n               ";
    for (int n = 0; n < nofNodes; n++)
        out << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << "   ";

    out << "\n";
    out << "            </DataArray>\n";
    out << "         </Points>\n";

    // CELLS SECTION
    out << "         <Cells>\n";
    out << "            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n              ";

    for (int c = 0; c < nofNodes; c++)
        out << c << "   ";
    out << "\n";

    out << "            </DataArray>\n";
    out << "            <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n              ";
    for (int c = 1; c < nofNodes + 1; c++)
        out << c << " ";

    out << "\n";
    out << "            </DataArray>\n";

    out << "            <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n              ";
    for (int c = 0; c < nofNodes; c++)
        out << "1 ";

    out << "\n";
    out << "            </DataArray>\n";

    out << "         </Cells>\n";

    // write data section
    out << "         <PointData Scalars=\"Scalars\"> \n";
    for (int s = 0; s < (int)datanames.size(); ++s) {
        out << "           <DataArray type=\"Float64\" Name=\"" << datanames[s] << "\" format=\"ascii\"> \n";

        for (int d = 0; d < (int)nodedata[s].size(); d++)
            out << nodedata[s][d] << " ";

        out << "\n          </DataArray>\n";
    }
    out << "         </PointData>\n";
    out << "      </Piece>\n";
    out << "   </UnstructuredGrid>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
    out.close();
    UBLOG(logDEBUG1, "WbWriterVtkXmlASCII::writeNodesWithNodeData to " << vtkfilename << " - end");

    return vtkfilename;
}

//! \}
