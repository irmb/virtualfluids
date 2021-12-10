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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file WbWriterAvsBinary.h
//! \ingroup writer
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef WBWRITERAVSBINARY_H
#define WBWRITERAVSBINARY_H

#include <basics/writer/WbWriter.h>

class WbWriterAvsBinary : public WbWriter
{
public:
    static WbWriterAvsBinary *getInstance()
    {
        static WbWriterAvsBinary instance;
        return &instance;
    }

    WbWriterAvsBinary(const WbWriterAvsBinary &) = delete;
    const WbWriterAvsBinary &operator=(const WbWriterAvsBinary &) = delete;

private:
    WbWriterAvsBinary() = default;

public:
    std::string getFileExtension() override { return ".bin.inp"; }

    //////////////////////////////////////////////////////////////////////////
    // lines
    //     0 ---- 1
    // nodenumbering must start with 0!
    std::string writeLines(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                           std::vector<UbTupleInt2> &lines) override;

    //////////////////////////////////////////////////////////////////////////
    // triangles
    // cell numbering:
    //                    2
    //
    //                  0---1
    // nodenumbering must start with 0!
    std::string writeTriangles(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                               std::vector<UbTuple<int, int, int>> &triangles) override;
    std::string writeTrianglesWithNodeData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                           std::vector<UbTupleInt3> &cells, std::vector<std::string> &datanames,
                                           std::vector<std::vector<double>> &nodedata) override;

    //////////////////////////////////////////////////////////////////////////
    // quads
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
};

#endif // WBWRITERAVSBINARY_H
