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
#include <basics/utilities/UbFileOutputBinary.h>
#include <basics/utilities/UbSystem.h>
#include <cstring>

using namespace std;

/*==========================================================*/
UbFileOutputBinary::UbFileOutputBinary(const string &filename, const bool &createPath)
{
    this->filename         = filename;
    this->commentindicator = 'C';

    outfile.open(filename.c_str(), ios::out | ios::binary);

    if (!outfile && createPath) {
        string path = UbSystem::getPathFromString(filename);
        if (path.size() > 0) {
            outfile.clear(); // flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
            UbSystem::makeDirectory(path);
            outfile.open(filename.c_str(), ios::out | ios::binary);
        }
    }

    if (!outfile)
        UB_THROW(UbException(UB_EXARGS, "couldn't open file:\n " + filename));
}
/*==========================================================*/
UbFileOutputBinary::UbFileOutputBinary(const string &filename, UbFileOutput::CREATEOPTION opt, const bool &createPath)
{
    this->filename         = filename;
    this->commentindicator = 'C';

    this->open(filename, opt);

    if (!this->open(filename, opt) && createPath) {
        string path = UbSystem::getPathFromString(filename);
        if (path.size() > 0) {
            outfile.clear(); // flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
            UbSystem::makeDirectory(path);

            this->open(filename, opt);
        }
    }

    if (!outfile)
        UB_THROW(UbException(UB_EXARGS, "couldn't open file:\n " + filename));
}
/*==========================================================*/
bool UbFileOutputBinary::open(const string &filename, UbFileOutput::CREATEOPTION opt)
{
    outfile.close();
    outfile.clear(); // setzt flags zurueck

    this->filename = filename;

    if (opt == UbFileOutput::OUTFILE)
        outfile.open(this->filename.c_str(), ios::out | ios::binary);
    else if (opt == UbFileOutput::APPENDFILE)
        outfile.open(this->filename.c_str(), ios::app | ios::binary);
    else if (opt == UbFileOutput::INANDOUTFILE)
        UB_THROW(UbException(UB_EXARGS, "undefined CREATEOPTION - INANDOUTFILE not possible for BINARY files"));
    else
        UB_THROW(UbException(UB_EXARGS, "undefined CREATEOPTION"));

    return outfile.is_open();
}
/*==========================================================*/
void UbFileOutputBinary::writeBool(const bool &value, const int & /*width*/)
{
    outfile.write((char *)&value, sizeof(bool));
}
/*==========================================================*/
void UbFileOutputBinary::writeDouble(const double &value, const int & /*width*/)
{
    outfile.write((char *)&value, sizeof(double));
}
/*==========================================================*/
void UbFileOutputBinary::writeFloat(const float &value, const int & /*width*/)
{
    outfile.write((char *)&value, sizeof(float));
}
/*==========================================================*/
void UbFileOutputBinary::setPrecision(const int & /*precision*/) { UB_THROW(UbException(UB_EXARGS, "no way")); }
/*==========================================================*/
int UbFileOutputBinary::getPrecision() { UB_THROW(UbException(UB_EXARGS, "no way")); }
/*==========================================================*/
void UbFileOutputBinary::writeInteger(const int &value, const int & /*width*/)
{
    outfile.write((char *)&value, sizeof(value));
}
/*==========================================================*/
void UbFileOutputBinary::writeSize_t(const std::size_t &value, const int & /*width*/)
{
    outfile.write((char *)&value, sizeof(value));
}
/*==========================================================*/
void UbFileOutputBinary::writeChar(const char &value, const int & /*width*/)
{
    outfile.write((char *)&value, sizeof(value));
}
/*==========================================================*/
void UbFileOutputBinary::writeString(const string &value, const int & /*width*/)
{
    char c              = '\0';
    unsigned int length = (unsigned)value.length();

    unsigned pos;
    // whitespaces und tabs am stringanfang uebergehen
    for (pos = 0; pos < length; pos++)
        if (value[pos] != ' ' && value[pos] != '\t')
            break;

    while (pos < length) {
        while (pos < length && value[pos] != ' ' && value[pos] != '\t' && value[pos] != '\0') {
            outfile.write((char *)&(value[pos++]), sizeof(char));
        }

        outfile.write(&c, sizeof(char));
        pos++;

        while (pos < length && (value[pos] == ' ' || value[pos] == '\t' || value[pos] == '\0')) {
            pos++;
        }
    }
}
/*==========================================================*/
void UbFileOutputBinary::writeStringOnly(const string & /*value*/)
{
    throw UbException(UB_EXARGS, "no way... causes to many errors"); // TODO: WTF?
}
/*==========================================================*/
void UbFileOutputBinary::writeLine(const std::string &value, const int & /*width*/)
{
    this->writeString(value);
    char c = '\n';
    outfile.write(&c, sizeof(char));
}
/*==========================================================*/
void UbFileOutputBinary::writeLine()
{
    char c = '\n';
    outfile.write(&c, sizeof(char));
}
/*==========================================================*/
void UbFileOutputBinary::writeCommentLine(const string &line)
{
    try {
        this->writeCommentLine(this->commentindicator, line);
    } catch (...) {
        UB_THROW(UbException(UB_EXARGS, "unknown error"));
    }
}
/*==========================================================*/
void UbFileOutputBinary::writeCommentLine(char indicator, const string &line)
{
    string dummy = indicator + line;
    this->writeLine(dummy);
}
/*==========================================================*/
void UbFileOutputBinary::writeCopyOfFile(const string &filename)
{
    ifstream infile(filename.c_str(), ios::in | ios::binary);
    if (!infile)
        UB_THROW(UbException(UB_EXARGS, "couldn't open file:\n " + filename));

    try {
        char c;
        while (infile.get(c)) {
            outfile.put(c); //=out<<c;
        }
        outfile.flush();
        infile.close();
    } catch (...) {
        UB_THROW(UbException(UB_EXARGS, "unknown error"));
    }
}

//! \}
