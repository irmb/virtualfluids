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
#ifndef WBWRITER_H
#define WBWRITER_H

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <basics/DataTypes.h>
#include <basics/utilities/UbException.h>
#include <basics/utilities/UbSystem.h>
#include <basics/utilities/UbTuple.h>

class WbWriter
{
public:
    //////////////////////////////////////////////////////////////////////////
    virtual ~WbWriter() = default;

    //////////////////////////////////////////////////////////////////////////
    // rein virtuelle Methoden
    virtual std::string getFileExtension() = 0;

    //////////////////////////////////////////////////////////////////////////
    // nodes
    virtual std::string writeNodes(const std::string & /*filename*/, std::vector<UbTupleFloat3> & /*nodes*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual std::string writeNodesWithNodeData(const std::string & /*filename*/, std::vector<UbTupleFloat3> & /*nodes*/,
                                               std::vector<std::string> & /*datanames*/,
                                               std::vector<std::vector<double>> & /*nodedata*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual std::string writeNodesWithNodeDataDouble(const std::string & /*filename*/,
                                                     std::vector<UbTupleDouble3> & /*nodes*/,
                                                     std::vector<std::string> & /*datanames*/,
                                                     std::vector<std::vector<double>> & /*nodedata*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }

    //////////////////////////////////////////////////////////////////////////
    // lines
    //     0 ---- 1
    // nodenumbering must start with 0!
    virtual std::string writeLines(const std::string & /*filename*/, std::vector<UbTupleFloat3> & /*nodes*/,
                                   std::vector<UbTupleInt2> & /*lines*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual std::string writePolyLines(const std::string & /*filename*/, real*  /*coordinatesX*/,
                                                            real*  /*coordinatesY*/, real* /*coordinatesZ*/, uint  /*numberOfCoordinates*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual std::string writePolyLines(const std::string & /*filename*/, std::vector<real>&  /*coordinatesX*/,
                                                            std::vector<real>&  /*coordinatesY*/,  std::vector<real>& /*coordinatesZ*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual std::string writeLinesWithNodeData(const std::string & /*filename*/, std::vector<UbTupleFloat3> & /*nodes*/,
                                               std::vector<UbTupleInt2> & /*lines*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual std::string writeLinesWithLineData(const std::string & /*filename*/, std::vector<UbTupleFloat3> & /*nodes*/,
                                               std::vector<UbTupleInt2> & /*lines*/, std::vector<std::string> & /*datanames*/,
                                               std::vector<std::vector<float>> & /*celldata*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    //////////////////////////////////////////////////////////////////////////
    // triangles
    // cell numbering:
    //                     2
    //
    //                  0 === 1
    // nodenumbering must start with 0!
    virtual std::string writeTriangles(const std::string & /*filename*/, std::vector<UbTupleFloat3> & /*nodes*/,
                                       std::vector<UbTupleInt3> & /*cells*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual std::string writeTrianglesWithNodeData(const std::string & /*filename*/,
                                                   std::vector<UbTupleFloat3> & /*nodes*/,
                                                   std::vector<UbTupleInt3> & /*cells*/,
                                                   std::vector<std::string> & /*datanames*/,
                                                   std::vector<std::vector<double>> & /*nodedata*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }

    //////////////////////////////////////////////////////////////////////////
    // quads
    // cell numbering:
    //                  3---2
    //                  |   |
    //                  0---1
    // nodenumbering must start with 0!
    virtual std::string writeQuads(const std::string & /*filename*/, std::vector<UbTupleFloat3> & /*nodes*/,
                                   std::vector<UbTupleInt4> & /*cells*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual std::string writeQuadsWithNodeData(const std::string & /*filename*/, std::vector<UbTupleFloat3> & /*nodes*/,
                                               std::vector<UbTupleInt4> & /*cells*/,
                                               std::vector<std::string> & /*datanames*/,
                                               std::vector<std::vector<double>> & /*nodedata*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual std::string writeQuadsWithCellData(const std::string & /*filename*/, std::vector<UbTupleFloat3> & /*nodes*/,
                                               std::vector<UbTupleInt4> & /*cells*/,
                                               std::vector<std::string> & /*datanames*/,
                                               std::vector<std::vector<double>> & /*celldata*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual std::string writeQuadsWithNodeAndCellData(
        const std::string & /*filename*/, std::vector<UbTupleFloat3> & /*nodes*/, std::vector<UbTupleInt4> & /*cells*/,
        std::vector<std::string> & /*nodedatanames*/, std::vector<std::vector<double>> & /*nodedata*/,
        std::vector<std::string> & /*celldatanames*/, std::vector<std::vector<double>> & /*celldata*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }

    //////////////////////////////////////////////////////////////////////////
    // octs
    //     7 ---- 6
    //    /|     /|
    //   4 +--- 5 |
    //   | |    | |
    //   | 3 ---+ 2
    //   |/     |/
    //   0 ---- 1
    virtual std::string writeOcts(const std::string & /*filename*/, std::vector<UbTupleFloat3> & /*nodes*/,
                                  std::vector<UbTupleInt8> & /*cells*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual std::string writeOctsWithCellData(const std::string & /*filename*/, std::vector<UbTupleFloat3> & /*nodes*/,
                                              std::vector<UbTupleInt8> & /*cells*/,
                                              std::vector<std::string> & /*datanames*/,
                                              std::vector<std::vector<double>> & /*celldata*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }
    virtual std::string writeOctsWithNodeData(const std::string & /*filename*/, std::vector<UbTupleFloat3> & /*nodes*/,
                                              std::vector<UbTupleUInt8> & /*cells*/,
                                              std::vector<std::string> & /*datanames*/,
                                              std::vector<std::vector<double>> & /*nodedata*/)
    {
        throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
    }

private:
};

#endif // WBWRITER_H

//! \}
