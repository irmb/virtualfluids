//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBEQUAL_H
#define UBEQUAL_H

#include<cmath>

//////////////////////////////////////////////////////////////////////////
//isUbEqual<T1,T2>(a,b)
//vergleicht die gleichtheit der beiden werte a und b
//
//std-maessig wird hierfür der operator== verwendet
//
//Ausnahme: floating-points
//hier wird jeweils der "genauere typ zum ungenaueren gecastet und dann verglichen"
//e.g.: double d=1.2; int i=1; bool check = isUbEqual(d,i); -> true
//
//bei klassen muss hier operator== fuer const objecte implementiert sein!!!
//e.g.: bool operator==(const Test&) const { if(blabla) return true; else return false; }
//
//
//author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
//version 1.0 - 25.03.2008
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

//spezialisierung für diverse Typen-Tuples
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
