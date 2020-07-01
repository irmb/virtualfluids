#ifndef UB_INFINITY_H
#define UB_INFINITY_H
#include <limits>

#include <basics/utilities/UbLimits.h>
#include <basics/utilities/UbSystem.h>


//////////////////////////////////////////////////////////////////////////
//
//  UbNegInfinity
//  Anm: keine template klasse, da man am Ende eine Instanz "inf" verwendet
//       die in "verschiedene"(!!!) Typen konvertiert werden kann und nicht 
//       nur in den template Typ!
//  Note: The UbNegInfinity class cannot be instantiated on its own, but works 
//        as a base class for the Infinity class.
//////////////////////////////////////////////////////////////////////////
class UbNegInfinity
{
 public:
   //name Conversion operators 
   inline operator signed char() const { return UbLimits<signed char>::ninf(); }
   inline operator char()        const { return UbLimits<char>::ninf();        }
   inline operator wchar_t()     const { return UbLimits<wchar_t>::ninf();     }
   inline operator short()       const { return UbLimits<short>::ninf();       }
   inline operator int()         const { return UbLimits<int>::ninf();         }
   inline operator long()        const { return UbLimits<long>::ninf();        }
   inline operator float()       const { return UbLimits<float>::ninf();       }
   inline operator double()      const { return UbLimits<double>::ninf();      }
   inline operator long double() const { return UbLimits<long double>::ninf(); }

   // This function compares built-in data types with their largest possible value. The function
   // only works for built-in data types. The attempt to compare user-defined class types will
   // result in a compile time error.
   template< typename T >
   inline bool equal( const T& rhs ) const
   {
      UB_STATIC_ASSERT( std::numeric_limits<T>::is_specialized );
      return UbLimits<T>::ninf() == rhs;
   }
 protected:
    inline UbNegInfinity() {}

 private:
   UbNegInfinity( const UbNegInfinity& ninf );             //copy constructor (private & undefined)
   UbNegInfinity& operator=( const UbNegInfinity& ninf );  //copy assignment operator (private & undefined)
   void* operator&() const;                                //address operator (private & undefined)
};

//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================
template< typename T >
inline bool operator==( const UbNegInfinity& lhs, const T& rhs )
{
   return lhs.equal( rhs );
}
//*************************************************************************************************
template< typename T >
inline bool operator==( const T& lhs, const UbNegInfinity& rhs )
{
   return rhs.equal( lhs );
}
//*************************************************************************************************
template< typename T >
inline bool operator!=( const UbNegInfinity& lhs, const T& rhs )
{
   return !lhs.equal( rhs );
}
//*************************************************************************************************
template< typename T >
inline bool operator!=( const T& lhs, const UbNegInfinity& rhs )
{
   return !rhs.equal( lhs );
}

//////////////////////////////////////////////////////////////////////////
//
//  UbInfinity
//
//////////////////////////////////////////////////////////////////////////
class UbInfinity : public UbNegInfinity //um später -UbInfinity leichter zu implementieren!!!
{
 public:
   inline UbInfinity() 
      : UbNegInfinity()
    {}
   
   inline operator unsigned char()  const  { return UbLimits<unsigned char>::inf();  }
   inline operator signed char()    const  { return UbLimits<signed char>::inf();    }
   inline operator char()           const  { return UbLimits<char>::inf();           }
   inline operator wchar_t()        const  { return UbLimits<wchar_t>::inf();        }
   inline operator unsigned short() const  { return UbLimits<unsigned short>::inf(); }
   inline operator short()          const  { return UbLimits<short>::inf();          }
   inline operator unsigned int()   const  { return UbLimits<unsigned int>::inf();   }
   inline operator int()            const  { return UbLimits<int>::inf();            }
   inline operator unsigned long()  const  { return UbLimits<unsigned long>::inf();  }
   inline operator long()           const  { return UbLimits<long>::inf();           }
   inline operator float()          const  { return UbLimits<float>::inf();          }
   inline operator double()         const  { return UbLimits<double>::inf();         }
   inline operator long double()    const  { return UbLimits<long double>::inf();    }

   inline const UbNegInfinity& operator-() const { return static_cast<const UbNegInfinity&>( *this ); }

   /*==========================================================*/
   template< typename T >
   inline bool equal( const T& rhs ) const
   {
      UB_STATIC_ASSERT( std::numeric_limits<T>::is_specialized );
      return UbLimits<T>::inf() == rhs;
   }

 private:
   UbInfinity( const UbInfinity& inf );             //Copy constructor (private & undefined)
   UbInfinity& operator=( const UbInfinity& inf );  //Copy assignment operator (private & undefined)
   void* operator&() const;                         //Address operator (private & undefined)
};

//////////////////////////////////////////////////////////////////////////
//  GLOBAL OPERATORS
//////////////////////////////////////////////////////////////////////////
template< typename T >
inline bool operator==( const UbInfinity& lhs, const T& rhs );

template< typename T >
inline bool operator==( const T& lhs, const UbInfinity& rhs );

template< typename T >
inline bool operator!=( const UbInfinity& lhs, const T& rhs );

template< typename T >
inline bool operator!=( const T& lhs, const UbInfinity& rhs );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between an Infinity object and a built-in data type.
// \ingroup util
//
// This operator works only for built-in data types. The attempt to compare user-defined class
// types will result in a compile time error.
*/
template< typename T >
inline bool operator==( const UbInfinity& lhs, const T& rhs )
{
   return lhs.equal( rhs );
}
//*************************************************************************************************
template< typename T >
inline bool operator==( const T& lhs, const UbInfinity& rhs )
{
   return rhs.equal( lhs );
}
//*************************************************************************************************
template< typename T >
inline bool operator!=( const UbInfinity& lhs, const T& rhs )
{
   return !lhs.equal( rhs );
}
//*************************************************************************************************
template< typename T >
inline bool operator!=( const T& lhs, const UbInfinity& rhs )
{
   return !rhs.equal( lhs );
}
//*************************************************************************************************

//////////////////////////////////////////////////////////////////////////
//  GLOBAL INFINITY VALUE
//////////////////////////////////////////////////////////////////////////
namespace Ub
{
   //e.g. double x = UbSystem::inf;  float x = -Ub::inf; 
   const UbInfinity inf; 
} 

#endif //UB_INFINITY_H
