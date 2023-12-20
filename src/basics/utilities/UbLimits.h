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
#ifndef UB_LIMITS_H
#define UB_LIMITS_H

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <limits>

//////////////////////////////////////////////////////////////////////////
//  CLASS DEFINITION
//////////////////////////////////////////////////////////////////////////
template <typename T>
struct UbLimits {
};

//////////////////////////////////////////////////////////////////////////
//  SPECIALIZATIONS
//////////////////////////////////////////////////////////////////////////
template <>
struct UbLimits<unsigned char> {
    // return the largest possible positive unsigned char value
    static inline unsigned char inf() { return std::numeric_limits<unsigned char>::max(); }
};
//*************************************************************************************************
template <>
struct UbLimits<char> {
    // return the largest possible positive char value. */
    static inline char inf() { return std::numeric_limits<char>::max(); }
    // return the largest possible negative char value
    static inline char ninf() { return std::numeric_limits<char>::min(); }
};
//*************************************************************************************************
template <>
struct UbLimits<signed char> {
    // return the largest possible positive signed char value
    static inline signed char inf() { return std::numeric_limits<signed char>::max(); }

    // return The largest possible negative signed char value
    static inline signed char ninf() { return std::numeric_limits<signed char>::min(); }
};
//*************************************************************************************************
template <>
struct UbLimits<wchar_t> {
    // return The largest possible positive wchar_t value
    static inline wchar_t inf() { return std::numeric_limits<wchar_t>::max(); }
    // return The largest possible negative wchar_t value
    static inline wchar_t ninf() { return std::numeric_limits<wchar_t>::min(); }
};
//*************************************************************************************************
template <>
struct UbLimits<unsigned short> {
    // return The largest possible positive unsigned short value
    static inline unsigned short inf() { return std::numeric_limits<unsigned short>::max(); }
};
//*************************************************************************************************
template <>
struct UbLimits<short> {
    // return The largest possible positive short value
    static inline short inf() { return std::numeric_limits<short>::max(); }
    // return The largest possible negative short value
    static inline short ninf() { return std::numeric_limits<short>::min(); }
};
//*************************************************************************************************
template <>
struct UbLimits<unsigned int> {
    // return The largest possible positive unsigned int value
    static inline unsigned int inf() { return std::numeric_limits<unsigned int>::max(); }
};
//*************************************************************************************************
template <>
struct UbLimits<int> {
    // return The largest possible positive int value
    static inline int inf() { return std::numeric_limits<int>::max(); }

    // return The largest possible negative int value
    static inline int ninf() { return std::numeric_limits<int>::min(); }
};
//*************************************************************************************************
template <>
struct UbLimits<unsigned long> {
    // return The largest possible positive unsigned long value
    static inline unsigned long inf() { return std::numeric_limits<unsigned long>::max(); }
};
//*************************************************************************************************
template <>
struct UbLimits<long> {
    // return The largest possible positive long value
    static inline long inf() { return std::numeric_limits<long>::max(); }

    // return The largest possible negative long value
    static inline long ninf() { return std::numeric_limits<long>::min(); }
};
//*************************************************************************************************
template <>
struct UbLimits<float> {
    // return The largest possible positive float value
    static inline float inf() { return std::numeric_limits<float>::max(); }

    // return The largest possible negative float value
    static inline float ninf() { return -std::numeric_limits<float>::max(); }
};
//*************************************************************************************************
template <>
struct UbLimits<double> {
    // return The largest possible positive double value
    static inline double inf() { return std::numeric_limits<double>::max(); }
    // return The largest possible negative double value
    static inline double ninf() { return -std::numeric_limits<double>::max(); }
};
//*************************************************************************************************
template <>
struct UbLimits<long double> {
    // return The largest possible positive long double value
    static inline long double inf() { return std::numeric_limits<long double>::max(); }
    // return The largest possible negative long double value
    static inline long double ninf() { return -std::numeric_limits<long double>::max(); }
};

#endif // UB_LIMITS_H

//! \}
