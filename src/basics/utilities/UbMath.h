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
#ifndef UBMATH_H
#define UBMATH_H

#include <basics/utilities/UbEqual.h>
#include <basics/utilities/UbSystem.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

namespace UbMath
{
extern const double PI;

//////////////////////////////////////////////////////////////////////////
// Hilfsfunktion fuer Genauigkeit
template <typename T>
struct Epsilon {
};

//////////////////////////////////////////////////////////////////////////
//  SPECIALIZATIONS von Epsilon
//////////////////////////////////////////////////////////////////////////
template <>
struct Epsilon<double> {
    static inline double val() { return 1.0E-11; }
};
template <>
struct Epsilon<float> {
    static inline float val() { return 1.0E-7f; }
};
template <>
struct Epsilon<long double> {
    static inline long double val() { return 1.0E-15; }
};
template <>
struct Epsilon<int> {
    static inline int val() { return 0; }
};

/*=======================================================*/
// -------------------------------------------------------------------------------------------------
// Funktion berechnet den Logarithmus einer Zahl z bzgl. der Basis b
// -------------------------------------------------------------------------------------------------
template <typename T>
inline T log(const T &z, const T &base)
{
    if (::log(base) == 0)
        return 1.0f;
    return ::log(z) / ::log(base);
}
/*=======================================================*/
// double x = UbMath::getNegativeInfinity<double>();
template <typename T>
inline T getNegativeInfinity()
{
    // assert(std::numeric_limits<T>::has_infinity);
    UB_STATIC_ASSERT(std::numeric_limits<T>::has_infinity);
    return -std::numeric_limits<T>::infinity();
}
/*=======================================================*/
// double x = UbMath::getPositiveInfinity<double>();
template <typename T>
inline T getPositiveInfinity()
{
    // assert(std::numeric_limits<T>::has_infinity);
    UB_STATIC_ASSERT(std::numeric_limits<T>::has_infinity);
    return std::numeric_limits<T>::infinity();
}
/*=======================================================*/
// double x; bool b = UbMath::isInfinity(x);
template <typename T>
inline bool isInfinity(const T &value)
{
    if (value == getNegativeInfinity<T>())
        return true;
    if (value == getPositiveInfinity<T>())
        return true;
    return false;
}
/*=======================================================*/
// double x = UbMath::getNaN<double>(x);
template <typename T>
inline T getNaN()
{
    UB_STATIC_ASSERT(std::numeric_limits<T>::has_quiet_NaN);
    return std::numeric_limits<T>::quiet_NaN();
}
/*=======================================================*/
// double x; bool b = UbMath::isNaN(x);
// x!=x liefert bei #QNAN "true"!
template <typename T>
inline bool isNaN(const T &x)
{
    UB_STATIC_ASSERT(std::numeric_limits<T>::has_quiet_NaN);
    return (x != x);
}
/*=======================================================*/
template <typename T>
inline T getEqualityEpsilon()
{
    return Epsilon<T>::val();
}
/*=======================================================*/
template <typename T>
inline bool zero(const T &value)
{
    return std::fabs(value) < Epsilon<T>::val();
    // return value >= -UbMath::EPSILON && value <= UbMath::EPSILON;
}
/*=======================================================*/
// spezialisierung fuer ints
template <>
inline bool zero(const int &value)
{
    return value == 0;
}
/*=======================================================*/
template <typename T1, typename T2>
inline bool zero(const T1 &value1, const T2 &value2)
{
    return !(!UbMath::zero(value1) || !UbMath::zero(value2));
}
/*=======================================================*/
template <typename T1, typename T2, typename T3>
inline bool zero(const T1 &value1, const T2 &value2, const T3 &value3)
{
    return !(!UbMath::zero(value1) || !UbMath::zero(value2, value3));
}
/*=======================================================*/
template <typename T>
inline bool negative(const T &value)
{
    return value < -Epsilon<T>::val();
}
/*=======================================================*/
template <typename T>
inline bool nonPositive(const T &value)
{
    return value <= Epsilon<T>::val();
}
/*=======================================================*/
template <typename T>
inline bool positive(const T &value)
{
    return value > +Epsilon<T>::val();
}
/*=======================================================*/
template <typename T>
inline bool nonNegative(const T &value)
{
    return value >= -Epsilon<T>::val();
}
/*=======================================================*/
template <typename T1, typename T2>
inline bool equal(const T1 &value, const T2 &reference)
{
    using High = typename UbEqualTrait<T1, T2>::High;
    return std::fabs(value - reference) < Epsilon<High>::val();
}
/*=======================================================*/
template <typename T1, typename T2, typename T3>
inline bool equal(const T1 &val1, const T2 &val2, const T3 &val3)
{
    return (UbMath::equal(val1, val2) && UbMath::equal(val1, val3));
}
/*=======================================================*/
template <typename T1, typename T2>
inline bool less(const T1 &value, const T2 &reference)
{
    using High = typename UbEqualTrait<T1, T2>::High;
    return value < reference - Epsilon<High>::val();
}
/*=======================================================*/
template <typename T1, typename T2>
inline bool lessEqual(const T1 &value, const T2 &reference)
{
    using High = typename UbEqualTrait<T1, T2>::High;
    return value <= reference + Epsilon<High>::val();
}
/*=======================================================*/
template <typename T1, typename T2>
inline bool greater(const T1 &value, const T2 &reference)
{
    using High = typename UbEqualTrait<T1, T2>::High;
    return value > reference + Epsilon<High>::val();
}
/*=======================================================*/
template <typename T1, typename T2>
inline bool greaterEqual(const T1 &value, const T2 &reference)
{
    using High = typename UbEqualTrait<T1, T2>::High;
    return value >= reference - Epsilon<High>::val();
}
/*=======================================================*/
template <typename T>
inline T round(const T &value, const int &decimalPlaces)
{
    return static_cast<T>(floor(value * pow(10.0, decimalPlaces) + 0.5) * pow(10.0, -decimalPlaces));
}
/*=======================================================*/
template <typename T>
inline int integerRounding(const T &value)
{
    return static_cast<int>(UbMath::zero(value) ? 0 : ((value < 0.0) ? (value - 0.5) : (value + 0.5)));
}
/*=======================================================*/
template <typename T>
inline T getRad(const T &degrees)
{
    return degrees * static_cast<T>(UbMath::PI / 180.0);
}
/*=======================================================*/
template <typename T>
inline T getDegrees(const T &rad)
{
    return rad * static_cast<T>(UbMath::PI / 180.0);
}
/*=======================================================*/
// aus wildmagic
template <typename T>
inline T ACos(const T &fValue)
{
    if (-1.0 < fValue) {
        if (fValue < 1.0)
            return static_cast<T>(acos(fValue));
        else
            return static_cast<T>(0.0);
    } else
        return static_cast<T>(PI);
}
/*=======================================================*/
template <typename T>
inline T ASin(const T &fValue)
{
    double HALF_PI = 0.5 * UbMath::PI;
    if (-1.0 < fValue) {
        if (fValue < 1.0)
            return static_cast<T>(asin(fValue));
        else
            return static_cast<T>(HALF_PI);
    } else
        return -static_cast<T>(HALF_PI);
}
/*=======================================================*/
template <typename T>
inline T invSqrt(const T &fValue)
{
    return static_cast<T>(1.0 / sqrt(fValue));
}

/*=======================================================*/
/**
 * Returns true, if specified values a and b are less both values c and d.
 * @param a the first value to check
 * @param b the second value to check
 * @param c the first value to check against
 * @param d the second value to check against
 * @return true, if specified values a and b are less both values c and d
 **/
template <typename T1, typename T2, typename T3, typename T4>
inline bool less2(const T1 &value1, const T2 &value2, T3 toBeLessAs1, T4 toBeLessAs2)
{
    return (less(value1, toBeLessAs1) && less(value1, toBeLessAs2) && less(value2, toBeLessAs1) &&
            less(value2, toBeLessAs2));
}
/*=======================================================*/
template <typename T1, typename T2, typename T3, typename T4>
inline bool greater2(const T1 &value1, const T2 &value2, T3 toBeGreaterAs1, T4 toBeGreaterAs2)
{
    return (greater(value1, toBeGreaterAs1) && greater(value1, toBeGreaterAs2) && greater(value2, toBeGreaterAs1) &&
            greater(value2, toBeGreaterAs2));
}
/*=======================================================*/
template <typename T1, typename T2, typename T3>
inline bool inClosedInterval(const T1 &value, const T2 &threshold1, const T3 &threshold2)
{
    if (threshold1 < threshold2) {
        return (greaterEqual(value, threshold1) && lessEqual(value, threshold2));
    }

    return (greaterEqual(value, threshold2) && lessEqual(value, threshold1));
}
/*=======================================================*/
template <typename T1, typename T2, typename T3>
inline bool inOpenInterval(const T1 &value, const T2 &threshold1, const T3 &threshold2)
{
    if (threshold1 < threshold2) {
        return (greater(value, threshold1) && less(value, threshold2));
    }

    return (greater(value, threshold2) && less(value, threshold1));
}
/*=======================================================*/
template <typename T1, typename T2, typename T3>
inline double adaptToClosedInterval(const T1 &value, const T2 &threshold1, const T3 &threshold2)
{
    if (threshold1 < threshold2) {
        if (less(value, threshold1))
            return threshold1;
        else if (greater(value, threshold2))
            return threshold2;
    } else {
        if (less(value, threshold2))
            return threshold2;
        else if (greater(value, threshold1))
            return threshold1;
    }
    return value;
}
/*=======================================================*/
// -------------------------------------------------------------------------------------------------
// Funktion berechnet den groessten gemeinsamen Teiler zweier Zahlen (MK)
// -------------------------------------------------------------------------------------------------
/*=======================================================*/
inline int calcGgt(int val1, int val2)
{
    if (val1 < val2)
        std::swap(val1, val2);
    int ggt = val2;
    while (ggt > 1) {
        if ((val1 % ggt) == 0 && (val2 % ggt) == 0)
            break;

        ggt -= 1;
    }
    return ggt;
}
/*=======================================================*/
// -------------------------------------------------------------------------------------------------
// Funktion berechnet den groessten gemeinsamen Teiler von drei Zahlen (MK)
// -------------------------------------------------------------------------------------------------
inline int calcGgt(int val1, const int &val2, int val3) { return UbMath::calcGgt(UbMath::calcGgt(val1, val2), val3); }
/*=======================================================*/
// returns the max of c2 values
// to avoid errors at mixed argument-types use: double myMax = max<double>(2,2.3);
template <typename T>
inline const T &max(const T &a1, const T &a2)
{
    return (a1 < a2) ? a2 : a1;
}
/*=======================================================*/
template <typename T>
inline const T &max(const T &a1, const T &a2, const T &a3)
{
    return max(max(a1, a2), a3);
}
/*=======================================================*/
template <typename T>
inline const T &max(const T &a1, const T &a2, const T &a3, const T &a4)
{
    return max(max(max(a1, a2), a3), a4);
}
/*=======================================================*/
template <typename T>
inline const T &min(const T &a1, const T &a2)
{
    return (a1 < a2) ? a1 : a2;
}
/*=======================================================*/
template <typename T>
inline const T &min(const T &a1, const T &a2, const T &a3)
{
    return min(min(a1, a2), a3);
}
/*=======================================================*/
template <typename T>
inline const T &min(const T &a1, const T &a2, const T &a3, const T &a4)
{
    return min(min(min(a1, a2), a3), a4);

    //       double tmp = a1;
    //       if(tmp>a2) tmp=a2;
    //       if(tmp>a3) tmp=a3;
    //       if(tmp>a4) tmp=a4;
    //       return tmp;
}

} // namespace UbMath

#endif

//! \}
