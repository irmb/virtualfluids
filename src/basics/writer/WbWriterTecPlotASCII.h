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
#ifndef WBWRITERTECPLOTASCII_H
#define WBWRITERTECPLOTASCII_H

#include <string>

#include <basics/writer/WbWriter.h>

class WbWriterTecPlotASCII : public WbWriter
{
public:
    static WbWriterTecPlotASCII *getInstance()
    {
        static WbWriterTecPlotASCII instance;
        return &instance;
    }

    WbWriterTecPlotASCII(const WbWriterTecPlotASCII &) = delete;
    const WbWriterTecPlotASCII &operator=(const WbWriterTecPlotASCII &) = delete;

private:
    WbWriterTecPlotASCII() : WbWriter()
    {
        if constexpr (sizeof(unsigned char) != 1)
            throw UbException(UB_EXARGS, "machine error char  type mismatch");
        if constexpr (sizeof(int) != 4)
            throw UbException(UB_EXARGS, "machine error int   type mismatch");
        if constexpr (sizeof(float) != 4)
            throw UbException(UB_EXARGS, "machine error float type mismatch");
    }

    static std::string pvdEndTag;

public:
    std::string getFileExtension() override { return ".ascii.dat"; }

    // write a metafile
    //    std::string writeCollection(const std::string& filename, const std::vector<std::string>& filenames, const
    //    double& timestep, const bool& sepGroups); std::string addFilesToCollection(const std::string& filename, const
    //    std::vector<std::string>& filenames, const double& timestep, const bool& sepGroups); std::string
    //    writeParallelFile(const std::string& filename,std::vector<std::string>& pieceSources,
    //    std::vector<std::string>& pointDataNames, std::vector<std::string>& cellDataNames);

    //////////////////////////////////////////////////////////////////////////
    // nodes
    //    std::string writeNodes(const std::string& filename,std::vector< UbTupleFloat3 >& nodes);
    //    std::string writeNodesWithNodeData(const std::string& filename,std::vector< UbTupleFloat3 >& nodes,
    //    std::vector<std::string >& datanames, std::vector<std::vector<double > >& nodedata);

    //////////////////////////////////////////////////////////////////////////
    // lines
    //     0 ---- 1
    // nodenumbering must start with 0!
    //    std::string writeLines(const std::string& filename,std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt2
    //    >& lines); std::string writeLinesWithNodeData(const std::string& filename,std::vector<UbTupleFloat3 >& nodes,
    //    std::vector<UbTupleInt2 >& lines, std::vector< std::string >& datanames, std::vector< std::vector< double > >&
    //    nodedata);
    //
    //////////////////////////////////////////////////////////////////////////
    // triangles
    //                    2
    //
    //                  0---1
    // nodenumbering must start with 0!
    //    std::string writeTriangles(const std::string& filename,std::vector<UbTupleFloat3 >& nodes,
    //    std::vector<UbTupleInt3 >& triangles); std::string writeTrianglesWithNodeData(const std::string&
    //    filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt3 >& cells, std::vector<std::string >&
    //    datanames, std::vector<std::vector<double > >& nodedata);

    //////////////////////////////////////////////////////////////////////////
    // 2D
    // cell numbering:
    //                  3---2
    //                  |   |
    //                  0---1
    // nodenumbering must start with 0!

    //    std::string writeQuads(const std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector<
    //    UbTupleInt4 >& cells); std::string writeQuadsWithNodeData(const std::string& filename,std::vector<
    //    UbTupleFloat3 >& nodes, std::vector< UbTupleInt4 >& cells, std::vector< std::string >& datanames, std::vector<
    //    std::vector< double > >& nodedata); std::string writeQuadsWithCellData(const std::string&
    //    filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt4 >& cells, std::vector< std::string >&
    //    datanames, std::vector< std::vector< double > >& celldata); std::string writeQuadsWithNodeAndCellData(const
    //    std::string& filename,std::vector< UbTupleFloat3 >& nodes, std::vector< UbTupleInt4 >& cells,
    //                                              std::vector< std::string >& nodedatanames, std::vector< std::vector<
    //                                              double > >& nodedata, std::vector< std::string >& celldatanames,
    //                                              std::vector< std::vector< double > >& celldata );

    //////////////////////////////////////////////////////////////////////////
    // octs
    //     7 ---- 6
    //    /|     /|
    //   4 +--- 5 |
    //   | |    | |
    //   | 3 ---+ 2
    //   |/     |/
    //   0 ---- 1
    std::string writeOctsU(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                           std::vector<UbTupleUInt8> &cells);
    // std::string writeOctsWithCellData(const std::string& filename,std::vector<UbTupleFloat3 >& nodes,
    // std::vector<UbTupleInt8 >& cells, std::vector<std::string >& datanames, std::vector<std::vector<double > >&
    // celldata);
    std::string writeOctsWithNodeData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                      std::vector<UbTupleUInt8> &cells, std::vector<std::string> &datanames,
                                      std::vector<std::vector<double>> &nodedata) override;
};

#endif // WBWRITERTECPLOTASCII_H

//! \}
