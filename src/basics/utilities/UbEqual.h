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
#ifndef UBEQUAL_H
#define UBEQUAL_H

#include <cmath>

//////////////////////////////////////////////////////////////////////////
//
//! \brief isUbEqual<T1,T2>(a,b)
//! Compares the equality of values a and b.
//!
//! By default operator== is used for this.
//!
//! Execption: floating-point variables
//! In these cases the type with higher precision is casted to the type of lower precision
//! and then the two values are compared.
//! e.g.: double d=1.2; int i=1; bool check = isUbEqual(d,i); -> true
//!
//! For classes operator== must be implemented for const objects!
//! e.g.: bool operator==(const Test&) const { if(blabla) return true; else return false; }
//
//////////////////////////////////////////////////////////////////////////

// std-trait, fuer alle nicht spezifischen typen:
template <typename T1, typename T2>
struct UbEqualTrait {
    using High = T1;
    using Low  = T1;
};

// std-trait, fuer gleiche T
template <typename T>
struct UbEqualTrait<T, T> {
    using High = T;
    using Low  = T;
};

// spezialisierung fuer diverse Typen-Tuples
template <>
struct UbEqualTrait<short, int> {
    using High = int;
    using Low  = short;
};
template <>
struct UbEqualTrait<short, long> {
    using High = long;
    using Low  = short;
};
template <>
struct UbEqualTrait<short, float> {
    using High = float;
    using Low  = short;
};
template <>
struct UbEqualTrait<short, double> {
    using High = double;
    using Low  = short;
};
template <>
struct UbEqualTrait<short, long double> {
    using High = long double;
    using Low  = short;
};

template <>
struct UbEqualTrait<int, short> {
    using High = int;
    using Low  = short;
};
template <>
struct UbEqualTrait<int, long> {
    using High = long;
    using Low  = int;
};
template <>
struct UbEqualTrait<int, float> {
    using High = float;
    using Low  = int;
};
template <>
struct UbEqualTrait<int, double> {
    using High = double;
    using Low  = int;
};
template <>
struct UbEqualTrait<int, long double> {
    using High = long double;
    using Low  = int;
};

template <>
struct UbEqualTrait<long, short> {
    using High = long;
    using Low  = short;
};
template <>
struct UbEqualTrait<long, int> {
    using High = long;
    using Low  = int;
};
template <>
struct UbEqualTrait<long, float> {
    using High = float;
    using Low  = long;
};
template <>
struct UbEqualTrait<long, double> {
    using High = double;
    using Low  = long;
};
template <>
struct UbEqualTrait<long, long double> {
    using High = long double;
    using Low  = long;
};

template <>
struct UbEqualTrait<float, short> {
    using High = float;
    using Low  = short;
};
template <>
struct UbEqualTrait<float, int> {
    using High = float;
    using Low  = int;
};
template <>
struct UbEqualTrait<float, long> {
    using High = float;
    using Low  = long;
};
template <>
struct UbEqualTrait<float, double> {
    using High = double;
    using Low  = float;
};
template <>
struct UbEqualTrait<float, long double> {
    using High = long double;
    using Low  = float;
};

template <>
struct UbEqualTrait<double, short> {
    using High = double;
    using Low  = short;
};
template <>
struct UbEqualTrait<double, int> {
    using High = double;
    using Low  = int;
};
template <>
struct UbEqualTrait<double, long> {
    using High = double;
    using Low  = long;
};
template <>
struct UbEqualTrait<double, float> {
    using High = double;
    using Low  = float;
};
template <>
struct UbEqualTrait<double, long double> {
    using High = long double;
    using Low  = double;
};

template <>
struct UbEqualTrait<long double, short> {
    using High = long double;
    using Low  = short;
};
template <>
struct UbEqualTrait<long double, int> {
    using High = long double;
    using Low  = int;
};
template <>
struct UbEqualTrait<long double, long> {
    using High = long double;
    using Low  = long;
};
template <>
struct UbEqualTrait<long double, float> {
    using High = long double;
    using Low  = float;
};
template <>
struct UbEqualTrait<long double, double> {
    using High = long double;
    using Low  = double;
};

//////////////////////////////////////////////////////////////////////////
// fuer Allgmeine-Typen ( operator== ):
template <typename T1, typename T2>
inline bool specific_equal(const T1 &a, const T2 &b)
{
    return a == b;
}

//////////////////////////////////////////////////////////////////////////
// fuer floating point build-in-type
// float.float
template </*float,float*/>
inline bool specific_equal<float, float>(const float &a, const float &b)
{
    return std::fabs(a - b) < 1E-8;
}

template </*double,double*/>
inline bool specific_equal<double, double>(const double &a, const double &b)
{
    return std::fabs(a - b) < 1E-13;
}

template </*long double,long double*/>
inline bool specific_equal<long double, long double>(const long double &a, const long double &b)
{
    return std::fabs(a - b) < 1E-16;
}

//////////////////////////////////////////////////////////////////////////
// globale isUbEqual - Funktion
template <typename T1, typename T2>
inline bool isUbEqual(const T1 &a, const T2 &b)
{
    using Low = typename UbEqualTrait<T1, T2>::Low;
    return specific_equal<Low, Low>(static_cast<Low>(a), static_cast<Low>(b));
}

//////////////////////////////////////////////////////////////////////////
// UbEqual-Functor
template <typename T1, typename T2>
struct UbEqual {
    bool operator()(const T1 &a, const T2 &b) { return isUbEqual(a, b); }
};

#endif // UBEQUAL_H

//! \}
