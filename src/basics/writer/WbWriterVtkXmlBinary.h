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
#ifndef WBWRITERVTKXMLBINARY_H
#define WBWRITERVTKXMLBINARY_H

#include <string>

#include <basics/writer/WbWriter.h>
#include <basics/DataTypes.h>



class WbWriterVtkXmlBinary : public WbWriter
{
public:
    static WbWriterVtkXmlBinary *getInstance()
    {
        static WbWriterVtkXmlBinary instance;
        return &instance;
    }

    WbWriterVtkXmlBinary(const WbWriterVtkXmlBinary &) = delete;
    const WbWriterVtkXmlBinary &operator=(const WbWriterVtkXmlBinary &) = delete;

private:
    WbWriterVtkXmlBinary() : WbWriter()
    {
        if constexpr (sizeof(unsigned char) != 1)
            throw UbException(UB_EXARGS, "machine error char  type mismatch");
        if constexpr (sizeof(int) != 4)
            throw UbException(UB_EXARGS, "machine error int   type mismatch");
        if constexpr (sizeof(float) != 4)
            throw UbException(UB_EXARGS, "machine error float type mismatch");
    }

    static const std::string pvdEndTag;

public:
    std::string getFileExtension() override { return ".bin.vtu"; }

    // write a metafile
    std::string writeCollection(const std::string &filename, const std::vector<std::string> &filenames,
                                const double &timeStep, const bool &separateGroups) const;
    std::string writeCollectionForTimeSeries(const std::string &filename,
                                             const std::map<uint, std::vector<std::string>> &filesNamesForTimeSteps, bool separateGroups) const;
    std::string addFilesToCollection(const std::string &filename, const std::vector<std::string> &filenames,
                                     const double &timeStep, const bool &separateGroups) const;
    std::string writeParallelFile(const std::string &filename, std::vector<std::string> &pieceSources,
                                  std::vector<std::string> &pointDataNames, std::vector<std::string> &cellDataNames) const;

    //////////////////////////////////////////////////////////////////////////
    // nodes
    std::string writeNodes(const std::string &filename, std::vector<UbTupleFloat3> &nodes) override;
    std::string writeNodesWithNodeData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                       std::vector<std::string> &datanames,
                                       std::vector<std::vector<double>> &nodedata) override;

    //////////////////////////////////////////////////////////////////////////
    // lines
    //     0 ---- 1
    // nodenumbering must start with 0!
    std::string writeLines(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                           std::vector<UbTupleInt2> &lines) override;
    std::string writePolyLines(const std::string &filename, real* coordinatesX,
                                                    real* coordinatesY, real*coordinatesZ, uint numberOfCoordinates) override;
    std::string writePolyLines(const std::string & filename, std::vector<real>& coordinatesX,
                                                    std::vector<real>& coordinatesY,  std::vector<real>& coordinatesZ) override;
    // std::string writeLinesWithNodeData(const std::string& filename,std::vector<UbTupleFloat3 >& nodes,
    // std::vector<UbTupleInt2 >& lines, std::vector< std::string >& datanames, std::vector< std::vector< double > >&
    // nodedata);
    // FIXME: hides function in base class

    std::string writeLinesWithLineData(const std::string &filename, std::vector<UbTupleFloat3> &nodes, std::vector<UbTupleInt2> &lines,
                                       std::vector<std::string> &datanames, std::vector<std::vector<float>> &celldata) override;

    //////////////////////////////////////////////////////////////////////////
    // triangles
    //                    2
    //
    //                  0---1
    // nodenumbering must start with 0!
    std::string writeTriangles(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                               std::vector<UbTupleInt3> &triangles) override;
    std::string writeTrianglesWithNodeData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                           std::vector<UbTupleInt3> &cells, std::vector<std::string> &datanames,
                                           std::vector<std::vector<double>> &nodedata) override;

    //////////////////////////////////////////////////////////////////////////
    // 2D
    // cell numbering:
    //                  3---2
    //                  |   |
    //                  0---1
    // nodenumbering must start with 0!

    std::string writeQuads(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                           std::vector<UbTupleInt4> &cells) override;
    std::string writeQuadsWithNodeData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                       std::vector<UbTupleInt4> &cells, std::vector<std::string> &datanames,
                                       std::vector<std::vector<double>> &nodedata) override;
    std::string writeQuadsWithCellData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                       std::vector<UbTupleInt4> &cells, std::vector<std::string> &datanames,
                                       std::vector<std::vector<double>> &celldata) override;
    std::string writeQuadsWithNodeAndCellData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                              std::vector<UbTupleInt4> &cells, std::vector<std::string> &nodedatanames,
                                              std::vector<std::vector<double>> &nodedata,
                                              std::vector<std::string> &celldatanames,
                                              std::vector<std::vector<double>> &celldata) override;

    //////////////////////////////////////////////////////////////////////////
    // octs
    //     7 ---- 6
    //    /|     /|
    //   4 +--- 5 |
    //   | |    | |
    //   | 3 ---+ 2
    //   |/     |/
    //   0 ---- 1
    std::string writeOcts(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                          std::vector<UbTupleInt8> &cells) override;
    std::string writeOctsWithCellData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                      std::vector<UbTupleInt8> &cells, std::vector<std::string> &datanames,
                                      std::vector<std::vector<double>> &celldata) override;
    std::string writeOctsWithNodeData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                      std::vector<UbTupleUInt8> &cells, std::vector<std::string> &datanames,
                                      std::vector<std::vector<double>> &nodedata) override;

private:
};

#endif // WBWRITERVTKXMLBINARY_H

//! \}
