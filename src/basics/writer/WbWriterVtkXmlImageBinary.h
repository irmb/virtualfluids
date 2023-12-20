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
#ifndef WBWRITERVTKXMLIMAGEBINARY_H
#define WBWRITERVTKXMLIMAGEBINARY_H

#include <string>

#include <basics/writer/WbWriter.h>



class WbWriterVtkXmlImageBinary : public WbWriter
{
public:
    static WbWriterVtkXmlImageBinary *getInstance()
    {
        static WbWriterVtkXmlImageBinary instance;
        return &instance;
    }

    WbWriterVtkXmlImageBinary(const WbWriterVtkXmlImageBinary &) = delete;
    const WbWriterVtkXmlImageBinary &operator=(const WbWriterVtkXmlImageBinary &) = delete;

private:
    WbWriterVtkXmlImageBinary() : WbWriter()
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
    std::string getFileExtension() override { return ".bin.vti"; }

    // write a metafile
    std::string writeCollection(const std::string &filename, const std::vector<std::string> &filenames,
                                const double &timestep, const bool &sepGroups);
    std::string addFilesToCollection(const std::string &filename, const std::vector<std::string> &filenames,
                                     const double &timestep, const bool &sepGroups);
    std::string writeParallelFile(const std::string &filename, const UbTupleInt6 &wholeExtent, const UbTupleFloat3 &origin, const UbTupleFloat3 &spacing, 
                                std::vector<std::string> &pieceSources, std::vector<UbTupleInt6> &pieceExtents,
                                std::vector<std::string> &pointDataNames, std::vector<std::string> &cellDataNames);

    //////////////////////////////////////////////////////////////////////////
    // nodes
    std::string writeNodesWithNodeData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                       std::vector<std::string> &datanames,
                                       std::vector<std::vector<double>> &nodedata) override;

    //////////////////////////////////////////////////////////////////////////
    // octs
    //     7 ---- 6
    //    /|     /|
    //   4 +--- 5 |
    //   | |    | |
    //   | 3 ---+ 2
    //   |/     |/
    //   0 ---- 1
    std::string writeOctsWithCellData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                      std::vector<UbTupleInt8> &cells, std::vector<std::string> &datanames,
                                      std::vector<std::vector<double>> &celldata) override;
    std::string writeOctsWithNodeData(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                      std::vector<UbTupleUInt8> &cells, std::vector<std::string> &datanames,
                                      std::vector<std::vector<double>> &nodedata) override;
    void writeData(const std::string &vtkfilename,
                                            std::vector<std::string> &pointDataNames, std::vector<std::string> &cellDataNames,
                                            std::vector<std::vector<double>> &nodedata, std::vector<std::vector<double>> &celldata, 
                                            UbTupleInt6 &wholeExtent,
                                            UbTupleFloat3 &origin, UbTupleFloat3 &spacing, UbTupleInt6 &extent, unsigned int precision=6);

private:
    void getMetaDataOfImage(std::vector<UbTupleFloat3> &nodes, UbTupleFloat3& origin, UbTupleFloat3& spacing, UbTupleInt6& extent);
};

#endif // WBWRITERVTKXMLIMAGEBINARY_H

//! \}
