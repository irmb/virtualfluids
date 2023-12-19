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
#include <algorithm>
#include <basics/utilities/UbFileInputASCII.h>
#include <cstring>

using namespace std;

UbFileInputASCII::UbFileInputASCII(string filename)
{
    this->filename         = filename;
    this->commentindicator = 'C';

    infile.open(filename.c_str());

    // if(!infile) UB_THROW( UbException((string)("UbFileInputASCII::UbFileInputASCII(string filename, int how) couldn't
    // open file:\n "+filename)) );
}
/*==========================================================*/
bool UbFileInputASCII::open(string filename)
{
    infile.close();
    infile.clear(); // setzt flags zurueck

    this->filename = filename;
    infile.open(this->filename.c_str());

    return infile.is_open();
}
/*==========================================================*/
int UbFileInputASCII::readInteger()
{
    int dummy;
    infile >> dummy;
    return dummy;
}
/*==========================================================*/
long long UbFileInputASCII::readLongLong()
{
    long long dummy;
    infile >> dummy;
    return dummy;
}
/*==========================================================*/
string UbFileInputASCII::getFileName() { return this->filename; }

/*==========================================================*/
void UbFileInputASCII::skipLine()
{
    string dummy;
    getline(infile, dummy);
}
/*==========================================================*/
void UbFileInputASCII::readLine()
{
    string dummy;
    getline(infile, dummy);
}
/*==========================================================*/
string UbFileInputASCII::readStringLine()
{
    string dummy;
    getline(infile, dummy);
    return dummy;
}
/*==========================================================*/
string UbFileInputASCII::readLineTill(char stop)
{
    string dummy;
    getline(infile, dummy, stop);
    return dummy;
}
/*==========================================================*/
string UbFileInputASCII::parseString()
{
    string dummy;
    getline(infile, dummy, ' ');
    return dummy;
}
/*==========================================================*/
double UbFileInputASCII::readDouble()
{
    double dummy;
    infile >> dummy;
    return dummy;
}
/*==========================================================*/
float UbFileInputASCII::readFloat()
{
    float dummy;
    infile >> dummy;
    return dummy;
}
/*==========================================================*/
string UbFileInputASCII::readString()
{
    string dummy;
    infile >> dummy;
    return dummy;
}
/*==========================================================*/
char UbFileInputASCII::readChar()
{
    int dummy;
    infile >> dummy;
    return (char)dummy;
}
/*==========================================================*/
std::size_t UbFileInputASCII::readSize_t()
{
    std::size_t dummy;
    infile >> dummy;
    return dummy;
}
/*==========================================================*/
void UbFileInputASCII::setPosAfterLineWithString(const string &var)
{
    infile.seekg(0L, ios::beg); // Positionszeiger der Datei auf den Anfang setzen
    char line[512];
    do {
        infile.getline(line, 512);
        if (infile.eof())
            UB_THROW(UbException(UB_EXARGS, "error at reading in file \"" + filename + "\" -> string " + var +
                                                " wasn't found in " + this->filename));
    } while (strstr(line, var.c_str()) != line); // Ende Schleife, wenn varname ganz in zeile vorkommt
}
/*==========================================================*/
bool UbFileInputASCII::containsString(const string &var)
{
    infile.clear(); // setzt den EOF-Status zurueck (wird durch infile.seekg() NICHT getan!!!)

    infile.seekg(0L, ios::beg); // Positionszeiger der Datei auf den Anfang setzen
    char line[512];
    do {
        infile.getline(line, 512);
        if (infile.eof())
            return false;
    } while (strstr(line, var.c_str()) != line); // Ende Schleife, wenn varname ganz in zeile vorkommt

    return true;
}
/*==========================================================*/
int UbFileInputASCII::readIntegerAfterString(const string &var)
// last change [29.6.2021] at [13:52]
// suchts in einer Datei nach varname und gibt den dahinter stehenden int-Wert zurueck
// z.B. timesteps 9
{
    infile.clear(); // setzt den EOF-Status zurueck (wird durch infile.seekg() NICHT getan!!!)

    infile.seekg(0L, ios::beg); // Positionszeiger der Datei auf den Anfang setzen

    char line[512];

    do {
        infile.getline(line, 512);
        if (infile.eof())
            UB_THROW(UbException(UB_EXARGS, "error at reading in file \"" + filename + "\" -> " + var +
                                                " wasn't found in " + this->filename));
    } while (strstr(line, var.c_str()) != line); // Ende Schleife, wenn varname ganz in zeile vorkommt

    std::string temp{ line };
    temp = temp.substr(var.size()); // zeile um "varname" kuerzen

    temp.erase(std::remove(temp.begin(), temp.end(), ' '), temp.end());  // remove whitespace
    temp.erase(std::remove(temp.begin(), temp.end(), '\t'), temp.end()); // remove tabs

    return std::stoi(temp);
}
/*==========================================================*/
// last change [29.6.2021] at [13:52]
// sucht in einer Datei nach varname und gibt den dahinter stehenden int-Wert zurueck
// z.B. nue 9.5
double UbFileInputASCII::readDoubleAfterString(const string &var)
{
    infile.clear(); // setzt den EOF-Status zurueck (wird durch infile.seekg() NICHT getan!!!)

    infile.seekg(0L, ios::beg); // Positionszeiger der Datei auf den Anfang setzen

    char line[512];

    do {
        infile.getline(line, 512);
        if (infile.eof())
            UB_THROW(UbException(UB_EXARGS, "error at reading in file \"" + filename + "\" -> " + var +
                                                " wasn't found in " + this->filename));
    } while (/*!strncmp(varname,line,sizeof(varname))==0*/ strstr(line, var.c_str()) !=
             line); // Ende Schleife, wenn varname ganz in zeile vorkommt

    std::string temp{ line };
    temp = temp.substr(var.size()); // zeile um "varname" kuerzen

    temp.erase(std::remove(temp.begin(), temp.end(), ' '), temp.end());  // remove whitespace
    temp.erase(std::remove(temp.begin(), temp.end(), '\t'), temp.end()); // remove tabs

    return std::stod(temp);
}
/*==========================================================*/
// last change [29.6.2021] at [13:52]
// liefert string-Wert der hinter dem uebergebenen char feld in der datei infile steht
// zudem wird der wert in die uebergebene variable value uebertragen (falls man das ergebniss als char benoetig)
string UbFileInputASCII::readStringAfterString(const string &var) //,char *value)
{
    infile.clear(); // setzt den EOF-Status zurueck (wird durch infile.seekg() NICHT getan!!!)

    infile.seekg(0L, ios::beg); // Positionszeiger der Datei auf den Anfang setzen

    char line[512];
    // string line_copy[512];

    do {
        infile.getline(line, 512);
        if (infile.eof())
            UB_THROW(UbException(UB_EXARGS, "error at reading in file \"" + filename + "\" -> " + var +
                                                " wasn't found in " + this->filename));
    } while (strstr(line, var.c_str()) != line); // Ende Schleife, wenn varname ganz in zeile vorkommt

    std::string temp{ line };
    temp = temp.substr(var.size()); // zeile um "varname" kuerzen

    temp.erase(std::remove(temp.begin(), temp.end(), ' '), temp.end());  // remove whitespace
    temp.erase(std::remove(temp.begin(), temp.end(), '\t'), temp.end()); // remove tabs

    return temp;
}
/*==========================================================*/
// last change [10.3.2004] at [9:46]
// sucht in einer Datei nach varname und gibt den dahinter stehenden int-Wert zurueck
// z.B. nue 9.5
bool UbFileInputASCII::readBoolAfterString(const string &var)
{
    if (this->readStringAfterString(var) == "true")
        return true;
    else if (this->readStringAfterString(var) == "false")
        return false;
    else
        UB_THROW(UbException(UB_EXARGS, "error at reading in file \"" + filename + "\" -> expression after " + var +
                                            " is not equal to 'true' or 'false' in " + this->filename));
}
/*==========================================================*/
// last change [10.3.2004] at [9:46]
// sucht in einer Datei nach varname und gibt den dahinter stehenden int-Wert zurueck
// z.B. nue 9.5
bool UbFileInputASCII::readBool()
{
    string tmp = this->readString();
    if (tmp == "true")
        return true;
    else if (tmp == "false")
        return false;
    else
        UB_THROW(UbException(UB_EXARGS, "error at reading in file \"" + filename + "\" -> expression=\"" + tmp +
                                            "\" is not equal to 'true' or 'false' in " + this->filename));
}

//! \}
