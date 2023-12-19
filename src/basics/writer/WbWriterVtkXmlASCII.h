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
#ifndef WBWRITERVTKXMLASCII_H
#define WBWRITERVTKXMLASCII_H

#include <string>

#include <basics/writer/WbWriter.h>

class WbWriterVtkXmlASCII : public WbWriter
{
public:
    static WbWriterVtkXmlASCII *getInstance()
    {
        static WbWriterVtkXmlASCII instance;
        return &instance;
    }

    WbWriterVtkXmlASCII(const WbWriterVtkXmlASCII &) = delete;
    const WbWriterVtkXmlASCII &operator=(const WbWriterVtkXmlASCII &) = delete;

private:
    WbWriterVtkXmlASCII() : WbWriter()
    {
        if constexpr (sizeof(unsigned char) != 1)
            throw UbException(UB_EXARGS, "error char  type mismatch");
        if constexpr (sizeof(int) != 4)
            throw UbException(UB_EXARGS, "error int   type mismatch");
        if constexpr (sizeof(float) != 4)
            throw UbException(UB_EXARGS, "error float type mismatch");
    }

    static std::string pvdEndTag;

public:
    std::string getFileExtension() override { return ".ascii.vtu"; }

    // write a metafile
    std::string writeCollection(const std::string &filename, const std::vector<std::string> &filenames,
                                const double &timesteps,
                                const bool &sepGroups); // std::vector<double>& groups, std::vector<double>& parts);
    std::string addFilesToCollection(const std::string &filename, const std::vector<std::string> &filenames,
                                     const double &timestep, const bool &sepGroups);
    std::string writeParallelFile(const std::string &filename, std::vector<std::string> &pieceSources,
                                  std::vector<std::string> &pointDataNames, std::vector<std::string> &cellDataNames);

    //////////////////////////////////////////////////////////////////////////
    // nodes
    std::string writeNodes(const std::string &filename, std::vector<UbTupleFloat3> &nodes) override;
    std::string writeNodesWithNodeData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                       std::vector<std::string> &datanames,
                                       std::vector<std::vector<double>> &nodedata) override;
    std::string writeNodesWithNodeDataDouble(const std::string &filename, std::vector<UbTupleDouble3> &nodes,
                                             std::vector<std::string> &datanames,
                                             std::vector<std::vector<double>> &nodedata) override;

    //////////////////////////////////////////////////////////////////////////
    // lines
    //     0 ---- 1
    // nodenumbering must start with 0!
    std::string writeLines(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                           std::vector<UbTupleInt2> &lines) override;
    // std::string writeLinesWithNodeData(const std::string& filename,std::vector<UbTupleFloat3 >& nodes,
    // std::vector<UbTupleInt2 >& lines, std::vector< std::string >& datanames, std::vector< std::vector< double > >&
    // nodedata);
    // FIXME: hides function in base class

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

#endif // WBWRITERVTKXMLASCII_H

//! \}
