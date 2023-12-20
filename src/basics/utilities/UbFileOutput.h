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
#ifndef UBFILEOUTPUT_H
#define UBFILEOUTPUT_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <basics/utilities/UbException.h>

/*=========================================================================*/
/*  UbFileOutput                                                             */
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

class UbFileOutput
{
public:
    enum CREATEOPTION { OUTFILE = 0, INANDOUTFILE = 1, APPENDFILE = 2 };
    enum FILETYPE { ASCII, BINARY };

public:
    UbFileOutput() : filename("") {}
    UbFileOutput(const std::string &filename) : filename(filename) {}
    virtual ~UbFileOutput() { outfile.flush(); }

    virtual bool open(const std::string &filename, CREATEOPTION opt = OUTFILE) = 0;

    virtual bool operator!() { return !(outfile); }
    virtual bool isOpen() { return !(!(outfile)); }

    virtual void flush() { outfile.flush(); }
    virtual void close() { outfile.close(); }

    virtual void writeInteger(const int &value, const int &width = 0)        = 0;
    virtual void writeDouble(const double &value, const int &width = 0)      = 0;
    virtual void writeFloat(const float &value, const int &width = 0)        = 0;
    virtual void writeBool(const bool &value, const int &width = 0)          = 0;
    virtual void writeSize_t(const std::size_t &value, const int &width = 0) = 0;
    virtual void writeChar(const char &value, const int &width = 0)          = 0;
    virtual void writeString(const std::string &value, const int &width = 0) = 0;
    virtual void writeStringOnly(const std::string &value)                   = 0;
    virtual void writeLine(const std::string &value, const int &width = 0)   = 0;
    virtual void writeLine()                                                 = 0;

    virtual void writeCommentLine(const std::string &line)                 = 0;
    virtual void writeCommentLine(char indicator, const std::string &line) = 0;
    virtual void writeCopyOfFile(const std::string &filename)              = 0;

    virtual void setCommentIndicator(char commentindicator) { this->commentindicator = commentindicator; }

    virtual void setPrecision(const int &precision) = 0;
    virtual int getPrecision()                      = 0;

    // returns "ASCII", "BINARY"
    virtual FILETYPE getFileType() = 0;

    // returns file extension:
    // e.g. "./../test/ich.inp" -> "inp", "./../test/ich" -> ""
    virtual std::string getFileExtension()
    {
        std::size_t pos1 = filename.rfind("/");
        if (pos1 == std::string::npos)
            pos1 = 0;
        std::size_t pos2 = filename.rfind(".");
        if (pos2 != std::string::npos && pos2 > pos1)
            return filename.substr(pos2 + 1);

        return "";
    }

    virtual std::string getFileName() { return this->filename; }

protected:
    std::ofstream outfile;
    std::string filename;
    char commentindicator{ 'C' };
};

#endif // UBFILEOUTPUT_H

//! \}
