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
#ifndef UB_INFINITY_H
#define UB_INFINITY_H
#include <limits>

#include <basics/utilities/UbLimits.h>
#include <basics/utilities/UbSystem.h>

//////////////////////////////////////////////////////////////////////////
//!
//!  \brief UbNegInfinity
//!  \details Note: The UbNegInfinity class cannot be instantiated on its own, but works
//!        as a base class for the Infinity class.
//!
//////////////////////////////////////////////////////////////////////////

class UbNegInfinity
{
public:
    // name Conversion operators
    inline operator signed char() const { return UbLimits<signed char>::ninf(); }
    inline operator char() const { return UbLimits<char>::ninf(); }
    inline operator wchar_t() const { return UbLimits<wchar_t>::ninf(); }
    inline operator short() const { return UbLimits<short>::ninf(); }
    inline operator int() const { return UbLimits<int>::ninf(); }
    inline operator long() const { return UbLimits<long>::ninf(); }
    inline operator float() const { return UbLimits<float>::ninf(); }
    inline operator double() const { return UbLimits<double>::ninf(); }
    inline operator long double() const { return UbLimits<long double>::ninf(); }

    //! This function compares built-in data types with their largest possible value. The function
    //! only works for built-in data types. The attempt to compare user-defined class types will
    //! result in a compile time error.
    template <typename T>
    inline bool equal(const T &rhs) const
    {
        UB_STATIC_ASSERT(std::numeric_limits<T>::is_specialized);
        return UbLimits<T>::ninf() == rhs;
    }

protected:
    inline UbNegInfinity() = default;

private:
    UbNegInfinity(const UbNegInfinity &ninf);            // copy constructor (private & undefined)
    UbNegInfinity &operator=(const UbNegInfinity &ninf); // copy assignment operator (private & undefined)
    void *operator&() const;                             // address operator (private & undefined)
};

//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================
template <typename T>
inline bool operator==(const UbNegInfinity &lhs, const T &rhs)
{
    return lhs.equal(rhs);
}
//*************************************************************************************************
template <typename T>
inline bool operator==(const T &lhs, const UbNegInfinity &rhs)
{
    return rhs.equal(lhs);
}
//*************************************************************************************************
template <typename T>
inline bool operator!=(const UbNegInfinity &lhs, const T &rhs)
{
    return !lhs.equal(rhs);
}
//*************************************************************************************************
template <typename T>
inline bool operator!=(const T &lhs, const UbNegInfinity &rhs)
{
    return !rhs.equal(lhs);
}

//////////////////////////////////////////////////////////////////////////
//
//  UbInfinity
//
//////////////////////////////////////////////////////////////////////////
class UbInfinity : public UbNegInfinity // um spaeter -UbInfinity leichter zu implementieren!!!
{
public:
    inline UbInfinity() : UbNegInfinity() {}

    inline operator unsigned char() const { return UbLimits<unsigned char>::inf(); }
    inline operator signed char() const { return UbLimits<signed char>::inf(); }
    inline operator char() const { return UbLimits<char>::inf(); }
    inline operator wchar_t() const { return UbLimits<wchar_t>::inf(); }
    inline operator unsigned short() const { return UbLimits<unsigned short>::inf(); }
    inline operator short() const { return UbLimits<short>::inf(); }
    inline operator unsigned int() const { return UbLimits<unsigned int>::inf(); }
    inline operator int() const { return UbLimits<int>::inf(); }
    inline operator unsigned long() const { return UbLimits<unsigned long>::inf(); }
    inline operator long() const { return UbLimits<long>::inf(); }
    inline operator float() const { return UbLimits<float>::inf(); }
    inline operator double() const { return UbLimits<double>::inf(); }
    inline operator long double() const { return UbLimits<long double>::inf(); }

    inline const UbNegInfinity &operator-() const { return static_cast<const UbNegInfinity &>(*this); }

    /*==========================================================*/
    template <typename T>
    inline bool equal(const T &rhs) const
    {
        UB_STATIC_ASSERT(std::numeric_limits<T>::is_specialized);
        return UbLimits<T>::inf() == rhs;
    }

private:
    UbInfinity(const UbInfinity &inf);            // Copy constructor (private & undefined)
    UbInfinity &operator=(const UbInfinity &inf); // Copy assignment operator (private & undefined)
    void *operator&() const;                      // Address operator (private & undefined)
};

//////////////////////////////////////////////////////////////////////////
//  GLOBAL OPERATORS
//////////////////////////////////////////////////////////////////////////
template <typename T>
inline bool operator==(const UbInfinity &lhs, const T &rhs);

template <typename T>
inline bool operator==(const T &lhs, const UbInfinity &rhs);

template <typename T>
inline bool operator!=(const UbInfinity &lhs, const T &rhs);

template <typename T>
inline bool operator!=(const T &lhs, const UbInfinity &rhs);
//@}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Equality comparison between an Infinity object and a built-in data type.
//
// This operator works only for built-in data types. The attempt to compare user-defined class
// types will result in a compile time error.
*/
template <typename T>
inline bool operator==(const UbInfinity &lhs, const T &rhs)
{
    return lhs.equal(rhs);
}
//*************************************************************************************************
template <typename T>
inline bool operator==(const T &lhs, const UbInfinity &rhs)
{
    return rhs.equal(lhs);
}
//*************************************************************************************************
template <typename T>
inline bool operator!=(const UbInfinity &lhs, const T &rhs)
{
    return !lhs.equal(rhs);
}
//*************************************************************************************************
template <typename T>
inline bool operator!=(const T &lhs, const UbInfinity &rhs)
{
    return !rhs.equal(lhs);
}
//*************************************************************************************************

//////////////////////////////////////////////////////////////////////////
//  GLOBAL INFINITY VALUE
//////////////////////////////////////////////////////////////////////////
namespace Ub
{
// e.g. double x = UbSystem::inf;  float x = -Ub::inf;
const UbInfinity inf;
} // namespace Ub

#endif // UB_INFINITY_H

//! \}
