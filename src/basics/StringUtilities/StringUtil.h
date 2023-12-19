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
#ifndef STRINGUTIL_H
#define STRINGUTIL_H

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>



#define SSTR(x) static_cast<std::ostringstream &>((std::ostringstream() << std::dec << x)).str()

class StringUtil
{
public:
    static std::string findAndReplace(const std::string &source, const std::string &find,
                                                    const std::string &replace);
    static std::string makeUpper(const std::string &instring);
    static std::string makeLower(const std::string &instring);
    static std::vector<std::string> split(const std::string &input, const std::string &delim = " ");
    static bool contains(const std::string &source, const char *find);
    static std::string pad(const std::string &input, char pad, int length);
    static std::string trim(const std::string &input, const std::string &trim = std::string(" \t\n"));
    static int toInt(const std::string &input);
    static float toFloat(const std::string &input);
    static double toDouble(const std::string &input);
    static bool toBool(const std::string &input);
    static std::vector<int> toIntVector(const std::string &s);
    static std::vector<unsigned int> toUintVector(const std::string &s);
    static std::vector<bool> toBoolVector(const std::string &s);
    static std::vector<std::string> toStringVector(const std::string &s);
    static std::vector<double> toDoubleVector(const std::string &s);
    template <typename T>
    static std::string toString(const T &t);

    static bool endsWith(const std::string &input, const std::string &end);


   template<class T>
   static T fromString(const std::string& s)
   {
      std::istringstream stream (s);
      T t;
      stream >> t;
      return t;
   }

private:
    StringUtil() = default;

    StringUtil(const StringUtil &) = default;
 
    virtual ~StringUtil() = default;


    static bool toBool(bool &t, const std::string &input, std::ios_base &(*f)(std::ios_base &));
};

#endif

//! \}
