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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file UbEqual.h
//! \ingroup utilities
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef UBEQUAL_H
#define UBEQUAL_H

#include<cmath>

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

//std-trait, fuer alle nicht spezifischen typen:
template < typename T1, typename T2 >
struct UbEqualTrait
{
   typedef T1 High;
   typedef T1 Low;
};

//std-trait, fuer gleiche T
template < typename T >
struct UbEqualTrait< T, T >
{
   typedef T High;
   typedef T Low;
};

//spezialisierung fuer diverse Typen-Tuples
template<> struct UbEqualTrait< short, int >          { typedef int         High; typedef short  Low; };
template<> struct UbEqualTrait< short, long >         { typedef long        High; typedef short  Low; };
template<> struct UbEqualTrait< short, float >        { typedef float       High; typedef short  Low; };
template<> struct UbEqualTrait< short, double >       { typedef double      High; typedef short  Low; };
template<> struct UbEqualTrait< short, long double >  { typedef long double High; typedef short  Low; };

template<> struct UbEqualTrait< int, short >          { typedef int         High; typedef short  Low; };
template<> struct UbEqualTrait< int, long >           { typedef long        High; typedef int    Low; };
template<> struct UbEqualTrait< int, float >          { typedef float       High; typedef int    Low; };
template<> struct UbEqualTrait< int, double >         { typedef double      High; typedef int    Low; };
template<> struct UbEqualTrait< int, long double >    { typedef long double High; typedef int    Low; };

template<> struct UbEqualTrait< long, short >         { typedef long        High; typedef short  Low; };
template<> struct UbEqualTrait< long, int >           { typedef long        High; typedef int    Low; };
template<> struct UbEqualTrait< long, float >         { typedef float       High; typedef long   Low; };
template<> struct UbEqualTrait< long, double >        { typedef double      High; typedef long   Low; };
template<> struct UbEqualTrait< long, long double >   { typedef long double High; typedef long   Low; };

template<> struct UbEqualTrait< float, short >        { typedef float       High; typedef short  Low; };
template<> struct UbEqualTrait< float, int >          { typedef float       High; typedef int    Low; };
template<> struct UbEqualTrait< float, long >         { typedef float       High; typedef long   Low; };
template<> struct UbEqualTrait< float, double >       { typedef double      High; typedef float  Low; };
template<> struct UbEqualTrait< float, long double >  { typedef long double High; typedef float  Low; };

template<> struct UbEqualTrait< double, short >       { typedef double      High; typedef short  Low; };
template<> struct UbEqualTrait< double, int >         { typedef double      High; typedef int    Low; };
template<> struct UbEqualTrait< double, long >        { typedef double      High; typedef long   Low; };
template<> struct UbEqualTrait< double, float >       { typedef double      High; typedef float  Low; };
template<> struct UbEqualTrait< double, long double > { typedef long double High; typedef double Low; };

template<> struct UbEqualTrait< long double, short >  { typedef long double High; typedef short  Low; };
template<> struct UbEqualTrait< long double, int >    { typedef long double High; typedef int    Low; };
template<> struct UbEqualTrait< long double, long >   { typedef long double High; typedef long   Low; };
template<> struct UbEqualTrait< long double, float >  { typedef long double High; typedef float  Low; };
template<> struct UbEqualTrait< long double, double > { typedef long double High; typedef double Low; };

//////////////////////////////////////////////////////////////////////////
//fuer Allgmeine-Typen ( operator== ):
template< typename T1, typename T2 >
inline bool specific_equal(const T1& a, const T2& b) { return a==b; }

//////////////////////////////////////////////////////////////////////////
//fuer floating point build-in-type
//float.float
template< /*float,float*/>
inline bool specific_equal< float, float >(const float& a, const float& b) {  return std::fabs( a - b ) < 1E-8; }

template</*double,double*/>
inline bool specific_equal< double, double >(const double& a, const double& b) { return std::fabs( a - b ) < 1E-13; }

template</*long double,long double*/>
inline bool specific_equal< long double, long double >(const long double& a, const long double& b) { return std::fabs( a - b ) < 1E-16; }

//////////////////////////////////////////////////////////////////////////
//globale isUbEqual - Funktion
template< typename T1, typename T2 >
inline bool isUbEqual(const T1& a, const T2& b)
{
   typedef typename UbEqualTrait<T1,T2>::Low Low;
   return specific_equal< Low, Low >(static_cast< Low >( a ),static_cast< Low >( b ));
};

//////////////////////////////////////////////////////////////////////////
//UbEqual-Functor
template< typename T1, typename T2 >
struct UbEqual
{
   bool operator()(const T1& a, const T2& b)
   {
      return isUbEqual(a,b);
   }
};

#endif //UBEQUAL_H
