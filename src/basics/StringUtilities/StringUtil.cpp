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
//! \addtogroup StringUtilities
//! \ingroup basics
//! \{
//! \author Konstantin Kutscher, Soeren Textor, Sebastian Geller
//=======================================================================================
#include "StringUtil.h"

#include <string.h>
#include <regex>
#include <sstream>

std::string StringUtil::findAndReplace(const std::string &source, const std::string &find, const std::string &replace)
{
    std::string output = source;
    size_t j;
    for (; (j = output.find(find)) != std::string::npos;)
        output.replace(j, find.length(), replace);
    return output;
}

std::string StringUtil::makeUpper(const std::string &instring)
{
    std::string output = instring;
    transform(output.begin(), output.end(), output.begin(), ::toupper);
    return output;
}

std::string StringUtil::makeLower(const std::string &instring)
{
    std::string output = instring;
    transform(output.begin(), output.end(), output.begin(), ::tolower);
    return output;
}

bool StringUtil::contains(const std::string &source, const char *find) { return (0 != strstr(source.c_str(), find)); }

std::string StringUtil::pad(const std::string &input, char pad, int length)
{
    std::string outstring = input;
    for (int i = (int)outstring.length(); i < length; ++i)
        outstring += pad;
    return outstring;
}

std::string StringUtil::trim(const std::string &input, const std::string &trim /*= std::string(" \t\n")*/)
{
    if (input.size() == 0)
        return input;
    std::string temp              = "";
    std::string::size_type begpos = input.find_first_not_of(trim);
    if (begpos == std::string::npos) {
        return temp;
    } else {
        std::string::size_type endpos = input.find_last_not_of(trim);
        temp                          = input.substr(begpos, endpos - begpos + 1);
    }
    return temp;
}

int StringUtil::toInt(const std::string &input) { return std::stoi(input); }

float StringUtil::toFloat(const std::string &input) { return std::stof(input); }

double StringUtil::toDouble(const std::string &input) { return std::stod(input); }

bool StringUtil::toBool(const std::string &input)
{
    bool b{ false };
    std::string trimmedInput = trim(input);
    if (!toBool(b, trimmedInput, std::boolalpha))
        throw "StringUtils::toBool() - Not a bool: " + trimmedInput;
    return b;
}

bool StringUtil::toBool(bool &t, const std::string &input, std::ios_base &(*f)(std::ios_base &))
{
    std::istringstream iss(input);
    return !(iss >> f >> t).fail();
}

std::vector<std::string> split(const std::string &str, char delim = ' ')
{
    std::stringstream ss(str);
    std::string token;
    std::vector<std::string> list;
    while (std::getline(ss, token, delim)) {
        list.push_back(token);
    }
    return list;
}

std::vector<std::string> StringUtil::split(const std::string &input, const std::string &delim /*= " "*/)
{
    std::stringstream ss;
    ss << "[" << delim << "]";

    std::regex re(ss.str());
    std::sregex_token_iterator first{ input.begin(), input.end(), re, -1 },
        last; // the '-1' is what makes the regex split (-1 := what was not matched)
    std::vector<std::string> tokens{ first, last };
    tokens.erase(std::remove_if(tokens.begin(), tokens.end(), [](std::string &token) { return token.empty(); }),
                 tokens.end());

    return tokens;
}

std::vector<int> StringUtil::toIntVector(const std::string &input)
{
    std::vector<int> v;
    std::vector<std::string> inputEntries;
    inputEntries = split(input, " \n\t");
    for (std::string entry : inputEntries)
        if (entry != "")
            v.push_back(toInt(entry));
    return v;
}

std::vector<unsigned int> StringUtil::toUintVector(const std::string &input)
{
    std::vector<unsigned int> v;
    std::vector<std::string> inputEntries;
    inputEntries = split(input, " \n\t");
    for (std::string entry : inputEntries)
        if (entry != "")
            v.push_back(toInt(entry));
    return v;
}

std::vector<bool> StringUtil::toBoolVector(const std::string &input)
{
    std::vector<bool> v;
    std::vector<std::string> inputEntries;
    inputEntries = split(input, " \n\t");
    for (std::string entry : inputEntries) {
        bool b{ false };
        std::string trimmedInput = trim(input);
        if (toBool(b, trimmedInput, std::noboolalpha))
            v.push_back(b);
    }
    return v;
}

std::vector<std::string> StringUtil::toStringVector(const std::string &input) { return split(input, " \n\t"); }

std::vector<double> StringUtil::toDoubleVector(const std::string &input)
{
    std::vector<double> v;
    std::vector<std::string> inputEntries;
    inputEntries = split(input, " \n\t");
    for (std::string entry : inputEntries)
        if (entry != "")
            v.push_back(toDouble(entry));
    return v;
}

template <typename T>
std::string StringUtil::toString(const T &t)
{
    std::ostringstream stream;
    stream << t;
    return stream.str();
}

template std::string StringUtil::toString<int>(const int &t);

bool StringUtil::endsWith(const std::string &input, const std::string &end)
{
    if (input.length() >= end.length()) {
        return (0 == input.compare(input.length() - end.length(), end.length(), end));
    } else {
        return false;
    }
}
//! \}
