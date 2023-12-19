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
#ifndef UBFILEINPUTASCII_H
#define UBFILEINPUTASCII_H

#include <fstream>
#include <iostream>
#include <string>



#include <basics/utilities/UbException.h>
#include <basics/utilities/UbFileInput.h>

/*=========================================================================*/
/*  UbFileInputASCII                                                       */
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

class UbFileInputASCII : public UbFileInput
{
public:
    UbFileInputASCII() : UbFileInput() {}
    UbFileInputASCII(std::string filename);

    bool open(std::string filename) override;

    std::string getFileName() override;
    void skipLine() override; // Springt zur naechsten Zeile

    void readLine() override;
    std::string readStringLine() override;
    int readInteger() override; // Liest einen Int-Wert ein
    long long readLongLong();   // Liest einen long-Wert ein

    std::size_t readSize_t() override;
    double readDouble() override;                 // Liest einen double-Wert ein
    float readFloat() override;                   // Liest einen float-Wert ein
    bool readBool() override;                     // Liest einen bool-Wert ein
    char readChar() override;                     // Liest einen char-Wert ein
    std::string readString() override;            // Liest ein Wort ein
    std::string readLineTill(char stop) override; // Liest gesamte Zeile ein bis zu einem bestimmten Zeichen
    std::string parseString() override;

    bool containsString(const std::string &var) override;
    void setPosAfterLineWithString(const std::string &var) override;
    int readIntegerAfterString(const std::string &var) override;
    double readDoubleAfterString(const std::string &var) override;
    bool readBoolAfterString(const std::string &var) override;
    std::string readStringAfterString(const std::string &var) override;

    FILETYPE getFileType() override { return ASCII; }

    template <typename T>
    friend inline UbFileInputASCII &operator>>(UbFileInputASCII &file, T &data)
    {
        file.infile >> data;
        return file;
    }
};

#endif // UBFILEINPUTASCII_H

//! \}
