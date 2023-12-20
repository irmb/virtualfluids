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
#ifndef UBFILEINPUTBINARY_H
#define UBFILEINPUTBINARY_H

#include <fstream>
#include <iostream>
#include <string>

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbFileInput.h>

/*=========================================================================*/
/*  UbFileInputBinary                                                      */
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

class UbFileInputBinary : public UbFileInput
{
public:
    UbFileInputBinary() : UbFileInput() {}
    UbFileInputBinary(std::string filename);

    bool open(std::string filename) override;

    void skipLine() override; // Springt zur naechsten Zeile
    void readLine() override;
    std::string readStringLine() override;
    std::size_t readSize_t() override;
    int readInteger() override;                   // Liest einen Int-Wert ein
    double readDouble() override;                 // Liest einen double-Wert ein
    float readFloat() override;                   // Liest einen float-Wert ein
    bool readBool() override;                     // Liest einen bool-Wert ein
    char readChar() override;                     // Liest einen char-Wert ein
    std::string readString() override;            // Liest ein Wort ein
    std::string readLineTill(char stop) override; // Liest gesamte Zeile ein bis zu einem bestimmten Zeichen
    std::string parseString() override;           // Liest

    bool containsString(const std::string &var) override;
    void setPosAfterLineWithString(const std::string &var) override;
    int readIntegerAfterString(const std::string &var) override;
    double readDoubleAfterString(const std::string &var) override;
    bool readBoolAfterString(const std::string &var) override;
    std::string readStringAfterString(const std::string &var) override;

    FILETYPE getFileType() override { return BINARY; }

    template <typename T>
    friend inline UbFileInputBinary &operator>>(UbFileInputBinary &file, T &data)
    {
        file.infile.read((char *)&data, sizeof(T));
        return file;
    }

    template <typename T>
    void readVector(std::vector<T> &v)
    {
        size_t size = v.size();
        if (size > 0) {
            infile.read((char *)&v[0], sizeof(T) * size);
        }
    }
};

#endif

//! \}
