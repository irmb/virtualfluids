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
//! \addtogroup utilities
//! \ingroup basics
//! \{
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef UBFILEOUTPUTBINARY_H
#define UBFILEOUTPUTBINARY_H

#include <fstream>
#include <iomanip>
#include <iostream>

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbFileOutput.h>

/*=========================================================================*/
/*  UbFileOutputBinary                                                             */
/*                                                                         */
/**
...
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 23.11.04
*/

/*
usage: ...
*/

class UbFileOutputBinary : public UbFileOutput
{
public:
    UbFileOutputBinary() : UbFileOutput() {}
    UbFileOutputBinary(const std::string &filename, const bool &createPath = true);
    UbFileOutputBinary(const std::string &filename, UbFileOutput::CREATEOPTION opt, const bool &createPath);

    bool open(const std::string &filename, UbFileOutput::CREATEOPTION opt = OUTFILE) override;

    void writeInteger(const int &value, const int &width = 0) override;
    void writeDouble(const double &value, const int &width = 0) override;
    void writeFloat(const float &value, const int &width = 0) override;
    void writeBool(const bool &value, const int &width = 0) override;
    void writeChar(const char &value, const int &width = 0) override;
    void writeSize_t(const std::size_t &value, const int &width = 0) override;
    void writeString(const std::string &value, const int &width = 0) override;
    void writeStringOnly(const std::string &value) override;
    void writeLine(const std::string &value, const int &width = 0) override;
    void writeLine() override;
    void writeCommentLine(const std::string &line) override;
    void writeCommentLine(char indicator, const std::string &line) override;
    void writeCopyOfFile(const std::string &filename) override;

    void setPrecision(const int &precision) override;
    int getPrecision() override;

    FILETYPE getFileType() override { return BINARY; }

    template <typename T>
    friend inline UbFileOutputBinary &operator<<(UbFileOutputBinary &file, const T &data)
    {
        file.outfile.write((char *)&data, sizeof(T));
        return file;
    }

    template <typename T>
    void writeVector(std::vector<T> &v)
    {
        size_t size = v.size();
        if (size > 0) {
            outfile.write((char *)&v[0], sizeof(T) * size);
        }
    }
};

#endif

//! \}
