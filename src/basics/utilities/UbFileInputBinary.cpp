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
#include <basics/utilities/UbFileInputBinary.h>
#include <cstring>

using namespace std;

/*==========================================================*/
UbFileInputBinary::UbFileInputBinary(string filename)
{
    this->filename = filename;
    infile.open(filename.c_str(), ios::in | ios::binary);
}
/*==========================================================*/
bool UbFileInputBinary::open(string filename)
{
    infile.close();
    infile.clear(); // setzt flags zurueck

    this->filename = filename;
    infile.open(this->filename.c_str(), ios::in | ios::binary);

    return infile.is_open();
}
/*==========================================================*/
int UbFileInputBinary::readInteger()
{
    int dummy;
    infile.read((char *)&dummy, sizeof(int));
    return dummy;
}
/*==========================================================*/
std::size_t UbFileInputBinary::readSize_t()
{
    std::size_t dummy;
    infile.read((char *)&dummy, sizeof(std::size_t));
    return dummy;
}
/*==========================================================*/
double UbFileInputBinary::readDouble()
{
    double dummy;
    infile.read((char *)&dummy, sizeof(double));
    return dummy;
}
/*==========================================================*/
float UbFileInputBinary::readFloat()
{
    float dummy;
    infile.read((char *)&dummy, sizeof(float));
    return dummy;
}
/*==========================================================*/
char UbFileInputBinary::readChar()
{
    char dummy;
    infile.read((char *)&dummy, sizeof(char));
    return dummy;
}
/*==========================================================*/
string UbFileInputBinary::readString()
{
    char c;
    infile.read(&c, sizeof(char));
    while (c == ' ' || c == '\t')
        infile.read(&c, sizeof(char));

    string dummy;
    dummy += c;

    infile.read(&c, sizeof(char));
    while (c != '\0' && c != ' ' && c != '\t' && c != '\n') {
        dummy += c;
        infile.read(&c, sizeof(char));
    }
    return dummy;
}
/*==========================================================*/
bool UbFileInputBinary::readBool()
{
    bool dummy;
    infile.read((char *)&dummy, sizeof(bool));
    return dummy;
}
/*==========================================================*/
void UbFileInputBinary::skipLine()
{
    char c;
    do {
        infile.read(&c, sizeof(char));
    } while (c != '\n');
}
/*==========================================================*/
void UbFileInputBinary::readLine()
{
    char c;
    infile.read(&c, sizeof(char));
    while (c != '\n')
        infile.read(&c, sizeof(char));
}
/*==========================================================*/
string UbFileInputBinary::readStringLine()
{
    char c;
    string dummy;
    infile.read(&c, sizeof(char));
    while (c != '\n') {
        dummy += c;
        infile.read(&c, sizeof(char));
    }
    return dummy;
}
/*==========================================================*/
string UbFileInputBinary::readLineTill(char /*stop*/)
{
    UB_THROW(UbException(UB_EXARGS, "method makes no sense for binary streams"));
}
/*==========================================================*/
string UbFileInputBinary::parseString()
{
    UB_THROW(UbException(UB_EXARGS, "method makes no sense for binary streams"));
}
/*==========================================================*/
bool UbFileInputBinary::containsString(const string & /*var*/)
{
    UB_THROW(UbException(UB_EXARGS, "method makes no sense for binary streams"));
}
/*==========================================================*/
void UbFileInputBinary::setPosAfterLineWithString(const string & /*var*/)
{
    UB_THROW(UbException(UB_EXARGS, "method makes no sense for binary streams"));
}
/*==========================================================*/
int UbFileInputBinary::readIntegerAfterString(const string & /*var*/)
{
    UB_THROW(UbException(UB_EXARGS, "method makes no sense for binary streams"));
}
/*==========================================================*/
double UbFileInputBinary::readDoubleAfterString(const string & /*var*/)
{
    UB_THROW(UbException(UB_EXARGS, "method makes no sense for binary streams"));
}
/*==========================================================*/
string UbFileInputBinary::readStringAfterString(const string & /*var*/)
{
    UB_THROW(UbException(UB_EXARGS, "method makes no sense for binary streams"));
}
/*==========================================================*/
bool UbFileInputBinary::readBoolAfterString(const string & /*var*/)
{
    UB_THROW(UbException(UB_EXARGS, "method makes no sense for binary streams"));
}

//! \}
