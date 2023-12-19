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
#include <basics/utilities/UbFileOutputASCII.h>
#include <basics/utilities/UbInfinity.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbSystem.h>
#include <cstring>

using namespace std;

UbFileOutputASCII::UbFileOutputASCII(const string &filename, const bool &createPath, const int & /*precision*/)
    : UbFileOutput(filename)
{
    this->commentindicator = 'C';
    this->setPrecision(20);

    outfile.open(filename.c_str(), ios::out);

    if (!outfile && createPath) {
        outfile.clear(); // flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
        string path = UbSystem::getPathFromString(filename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            outfile.open(filename.c_str(), ios::out);
        }
    }

    if (!outfile)
        UB_THROW(UbException(UB_EXARGS, "couldn't open file:\n " + filename));
}
/*==========================================================*/
UbFileOutputASCII::UbFileOutputASCII(const std::string &filename, CREATEOPTION opt, const bool &createPath,
                                     const int &precision)
    : UbFileOutput(filename)
{
    this->commentindicator = 'C';
    this->setPrecision(precision);

    if (!this->open(filename, opt) && createPath) {
        outfile.clear(); // flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
        string path = UbSystem::getPathFromString(filename);
        if (path.size() > 0)
            UbSystem::makeDirectory(path);

        this->open(filename, opt);
    }

    if (!outfile)
        UB_THROW(UbException(UB_EXARGS, "couldn't open file:\n " + filename));
}
/*==========================================================*/
bool UbFileOutputASCII::open(const std::string &filename, CREATEOPTION opt)
{
    outfile.close();
    outfile.clear(); // setzt flags zurueck
    this->filename = filename;

    if (opt == UbFileOutput::OUTFILE)
        outfile.open(this->filename.c_str(), ios::out);
    else if (opt == UbFileOutput::INANDOUTFILE)
        outfile.open(this->filename.c_str(), ios::out | ios::in);
    else if (opt == UbFileOutput::APPENDFILE)
        outfile.open(this->filename.c_str(), ios::app);
    else
        UB_THROW(UbException(UB_EXARGS, "undefined CREATEOPTION"));

    return outfile.is_open();
}
/*==========================================================*/
void UbFileOutputASCII::writeBool(const bool &value, const int &width)
{
    outfile.width(width);
    if (value)
        outfile << "true ";
    else
        outfile << "false ";
}
/*==========================================================*/
void UbFileOutputASCII::writeDouble(const double &value, const int &width)
{
    outfile.width(width);
    // Problem: Ub::inf wird gerundet
    //         -> beim Einlesen ist der Wert evtl zu gross und es kommt murks raus
    //         -> max Laenge darstellen und gut ist
    if (UbMath::equal(value, (double)Ub::inf)) {
        ios_base::fmtflags flags = outfile.flags();
        outfile << setprecision(std::numeric_limits<double>::digits10 + 2);
        outfile << value << " ";
        outfile.flags(flags);
        return;
    }
    outfile << value << " ";
}
/*==========================================================*/
void UbFileOutputASCII::writeFloat(const float &value, const int &width)
{
    outfile.width(width);
    // Problem: Ub::inf wird gerundet
    //         -> beim Einlesen ist der Wert evtl zu gross und es kommt murks raus
    //         -> max Laenge darstellen und gut ist
    if (UbMath::equal(value, (float)Ub::inf)) {
        ios_base::fmtflags flags = outfile.flags();
        outfile << setprecision(std::numeric_limits<float>::digits10 + 2);
        outfile << value << " ";
        outfile.flags(flags);
        return;
    }
    outfile << value << " ";
}
/*==========================================================*/
void UbFileOutputASCII::setPrecision(const int &precision) { outfile << setprecision(precision); }
/*==========================================================*/
void UbFileOutputASCII::writeInteger(const int &value, const int &width)
{
    outfile.width(width);
    outfile << value << " ";
}
/*==========================================================*/
void UbFileOutputASCII::writeSize_t(const std::size_t &value, const int &width)
{
    outfile.width(width);
    outfile << value << " ";
}
/*==========================================================*/
void UbFileOutputASCII::writeChar(const char &value, const int &width)
{
    outfile.width(width);
    outfile << (int)value << " ";
}
/*==========================================================*/
void UbFileOutputASCII::writeString(const string &value, const int &width)
{
    outfile.width(width);
    outfile << value.c_str() << " ";
}
/*==========================================================*/
void UbFileOutputASCII::writeStringOnly(const string &value) { outfile << value.c_str(); }

/*==========================================================*/
void UbFileOutputASCII::writeLine(const string &value, const int &width)
{
    outfile.width(width);
    outfile << value.c_str() << endl;
}
/*==========================================================*/
void UbFileOutputASCII::writeLine() { outfile << endl; }
/*==========================================================*/
void UbFileOutputASCII::writeCommentLine(const string &line) { this->writeCommentLine(this->commentindicator, line); }
/*==========================================================*/
void UbFileOutputASCII::writeCommentLine(char indicator, const string &line)
{
    this->outfile << indicator << line << endl;
}
/*==========================================================*/
void UbFileOutputASCII::writeCopyOfFile(const string &filename)
{
    ifstream infile(filename.c_str());
    if (!infile)
        UB_THROW(UbException(UB_EXARGS, "couldn't open file:\n " + filename));

    try {
        char c;
        while (infile.get(c)) {
            outfile.put(c); //=out<<c;
        }
        outfile.flush();
        infile.close();
    } catch (std::exception &e) {
        UB_THROW(UbException(UB_EXARGS, "catched std::exception, error: " + (std::string)e.what()));
    } catch (...) {
        UB_THROW(UbException(UB_EXARGS, "unknown error"));
    }
}

//! \}
