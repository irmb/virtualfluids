//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBTUPLE_H
#define UBTUPLE_H

#include <iostream>
#include <string>
#include <ostream>

#ifdef CAB_RCF
   #include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif //CAB_RCF

/*=========================================================================*/
/*  UbTuple                                                             */
/*                                                                         */
/**
...
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 23.10.06
*/ 

/*
usage: ...
////Advanced UbTuple
//Bsp:
//// create and use tuple with only one field
//UbTuple<int,int,int,int,int> t1;
//val<1>(t1) += 42;
//std::cout << t1.v1() << std::endl;

//// create and use duo
//UbTuple<bool,int> t2;
//std::cout << val<1>(t2) << ", ";
//std::cout << t2.v1() << std::endl;

//// create and use triple
//UbTuple<bool,int,double> t3;
//val<1>(t3) = true;  // new values via: val< pos >(triple) = ...
//val<2>(t3) = 42;
//val<3>(t3) = 0.2;
//t3 = makeUbTuple(false, 23, 13.13);

//std::cout << val<1>(t3) << ", ";
//std::cout << val<2>(t3) << ", ";
//std::cout << val<3>(t3) << std::endl;

//// create and use quadruple
//UbType<bool,int,float,double> t4(true,42,13,1.95583);
//std::cout << val<4>(t4) << std::endl;        //<- option 2 (std)
//std::cout << t4.v2().v2().v2() << std::endl; //<- option 2
*/

//typeop.h
// primary template
/**********************************
* typeop1.hpp:
**********************************/
template <typename T>
class UbTypeOp    // primary template
{           
public:
   typedef T         ArgT;
   typedef T         BareT;
   typedef T const   ConstT;
   typedef T &       RefT;
   typedef T &       RefBareT;
   typedef T const & RefConstT;
};
/**** end of typeop1.hpp ****/

// partial specialization for const
/**********************************
* typeop2.hpp:
**********************************/
template <typename T>
class UbTypeOp <T const>  // partial specialization for const types
{
 public:
   typedef T const   ArgT;
   typedef T         BareT;
   typedef T const   ConstT;
   typedef T const & RefT;
   typedef T &       RefBareT;
   typedef T const & RefConstT;
};
/**** end of typeop2.hpp ****/

// partial specialization for references
/**********************************
* typeop3.hpp:
**********************************/
template <typename T>
class UbTypeOp <T&>        // partial specialization for references
{
public:
   typedef T &                           ArgT;
   typedef typename UbTypeOp<T>::BareT   BareT;
   typedef T const                       ConstT;
   typedef T &                           RefT;
   typedef typename UbTypeOp<T>::BareT & RefBareT;
   typedef T const &                     RefConstT;
};
/**** end of typeop3.hpp ****/

// full specialization for void
/**********************************
* typeop4.hpp:
**********************************/
template<>
class UbTypeOp <void>      // full specialization for void
{
public:
   typedef void       ArgT;
   typedef void       BareT;
   typedef void const ConstT;
   typedef void       RefT;
   typedef void       RefBareT;
   typedef void       RefConstT;
};
/**** end of typeop4.hpp ****/

template <typename T1, typename T2>
class  UbDuo;

template <typename T1, typename T2>
std::ostream& operator << (std::ostream& os, UbDuo<T1, T2> const& d1);

//duo1.hpp
template <typename T1, typename T2>
class UbDuo 
{
public:
   typedef T1 Type1;  // type of first field
   typedef T2 Type2;  // type of second field
   enum { N = 2 };    // number of fields

public:
   // constructors
   UbDuo() : value1(), value2() {  }
   UbDuo (T1 const & a, T2 const & b) : value1(a), value2(b)  {  }

   // for implicit type conversion during construction
   template <typename U1, typename U2>
   UbDuo (UbDuo<U1,U2> const & d) : value1(d.v1()), value2(d.v2()) {  }

   // for implicit type conversion during assignments
   template <typename U1, typename U2>
   UbDuo<T1, T2>& operator = (UbDuo<U1,U2> const & d) 
   {
      value1 = d.v1();//value1;
      value2 = d.v2();//value2;
      return *this;
   }

   // field access
   T1& v1()             { return value1; }
   T1 const& v1() const { return value1; }

   T2& v2()             { return value2; }
   T2 const& v2() const { return value2; }

#ifdef CAB_RCF
   template<class Archive>
   void serialize(Archive& ar, const unsigned int version)
   {
      ar & value1;
      ar & value2;
   }
#endif //CAB_RCF

private:
   T1 value1;         // value of first field
   T2 value2;         // value of second field
};

// comparison operators (allow mixed types):
template <typename T1, typename T2,typename U1, typename U2>
inline bool operator == (UbDuo<T1,T2> const& d1, UbDuo<U1,U2> const& d2)
{
   return d1.v1()==d2.v1() && d1.v2()==d2.v2();
}

template <typename T1, typename T2,typename U1, typename U2>
inline bool operator != (UbDuo<T1,T2> const& d1, UbDuo<U1,U2> const& d2)
{
   return !(d1==d2);
}

template <typename T1, typename T2,typename U1, typename U2>
inline bool operator < (UbDuo<T1,T2> const& d1, UbDuo<U1,U2> const& d2)
{
   if     (d1.v1() <  d2.v1() ) return true;
   else if(d1.v1() == d2.v1() ) return d1.v2() < d2.v2();

   return false;
}

template <typename T1, typename T2>
std::ostream& operator << (std::ostream& os, UbDuo<T1, T2> const& d1)
{
    os << d1.v1() << ", " << d1.v2();
    return os;
}

// convenience function  for creation and initialization
template <typename T1, typename T2> 
inline UbDuo<T1,T2> makeUbDuo(T1 const & a, T2 const & b)
{
   return UbDuo<T1,T2>(a,b);
}

//duo2.hpp
template <typename A, typename B, typename C>
class UbDuo<A, UbDuo<B,C> > 
{
public:
   typedef A          T1;           // type of first field
   typedef UbDuo<B,C> T2;           // type of second field
   enum { N = UbDuo<B,C>::N + 1 };  // number of fields
 
public:
   // constructors
   UbDuo() : value1(), value2() { }
   UbDuo (T1 const & a, T2 const & b) : value1(a), value2(b) { }

   // for implicit type conversion during construction
   template <typename U1, typename U2>
   UbDuo (UbDuo<U1,U2> const & d) : value1(d.v1()), value2(d.v2()) { }

   // for implicit type conversion during assignments
   template <typename U1, typename U2>
   UbDuo<T1, T2>& operator = (UbDuo<U1,U2> const & d) 
   { 
      value1 = d.v1();//value1;     
      value2 = d.v2();//value2;
      return *this;
   }

   // field access
   T1& v1()             { return value1; }
   T1 const& v1() const { return value1; }

   T2& v2()             { return value2; }
   T2 const& v2() const { return value2; }

#ifdef CAB_RCF
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & value1;
      ar & value2;
   }
#endif //CAB_RCF

private:
   T1 value1;         // value of first field
   T2 value2;         // value of second field
};

//duo3.hpp
// primary template for type of Nth field of (duo) T
template <int N, typename T>
class UbDuoT 
{
public:
   typedef void ResultT;    // in general, the result type is void
};

// specialization for 1st field of a plain duo
template <typename A, typename B>
class UbDuoT<1, UbDuo<A,B> > 
{
public:
   typedef A ResultT;
};

// specialization for 2nd field of a plain duo
template <typename A, typename B>
class UbDuoT<2, UbDuo<A,B> > 
{
public:
   typedef B ResultT;
};

// specialization for Nth field of a recursive duo
template <int N, typename A, typename B, typename C>
class UbDuoT<N, UbDuo<A, UbDuo<B,C> > > 
{
public:
   typedef typename UbDuoT<N-1, UbDuo<B,C> >::ResultT ResultT;
};

// specialization for 1st field of a recursive duo
template <typename A, typename B, typename C>
class UbDuoT<1, UbDuo<A, UbDuo<B,C> > > 
{
public:
   typedef A ResultT;
};

// specialization for 2nd field of a recursive duo
template <typename A, typename B, typename C>
class UbDuoT<2, UbDuo<A, UbDuo<B,C> > > 
{
public:
   typedef B ResultT;
};

//duo4.hpp
// primary template for value of Nth field of (duo) T
template <int N, typename T>
class DuoValue 
{
public:
   static void get(T&) {  }      // in general, we have no value
   static void get(T const&) { }
};

// specialization for 1st field of a plain duo
template <typename A, typename B>
class DuoValue<1, UbDuo<A, B> > 
{
public:
   static A& get(UbDuo<A, B> &d)             { return d.v1(); }
   static A const& get(UbDuo<A, B> const &d) { return d.v1(); }
};

// specialization for 2nd field of a plain duo
template <typename A, typename B>
class DuoValue<2, UbDuo<A, B> > 
{
public:
   static B& get(UbDuo<A, B> &d)             { return d.v2(); }
   static B const& get(UbDuo<A, B> const &d) { return d.v2(); }
};

// specialization for Nth field of recursive duo
template <int N, typename A, typename B, typename C>
struct DuoValue<N, UbDuo<A, UbDuo<B,C> > >
{
   static typename UbTypeOp<typename UbDuoT<N-1, UbDuo<B,C> >::ResultT>::RefT  get(UbDuo<A, UbDuo<B,C> > &d)
   { 
      return DuoValue<N-1, UbDuo<B,C> >::get(d.v2()); 
   }
   static typename UbTypeOp<typename UbDuoT<N-1, UbDuo<B,C> >::ResultT>::RefConstT  get(UbDuo<A, UbDuo<B,C> > const &d)
   { 
      return DuoValue<N-1, UbDuo<B,C> >::get(d.v2()); 
   }
};

// specialization for 1st field of recursive duo
template <typename A, typename B, typename C>
class DuoValue<1, UbDuo<A, UbDuo<B,C> > > 
{
public:
   static A& get(UbDuo<A, UbDuo<B,C> > &d)             { return d.v1(); }
   static A const& get(UbDuo<A, UbDuo<B,C> > const &d) { return d.v1(); }
};

// specialization for 2nd field of recursive duo
template <typename A, typename B, typename C>
class DuoValue<2, UbDuo<A, UbDuo<B,C> > > 
{
public:
   static B& get(UbDuo<A, UbDuo<B,C> > &d)             { return d.v2().v1(); }
   static B const& get(UbDuo<A, UbDuo<B,C> > const &d) { return d.v2().v1(); }
};

//duo5.hpp
// return Nth value of variable duo
template <int N, typename A, typename B> 
inline typename UbTypeOp<typename UbDuoT<N, UbDuo<A, B> >::ResultT>::RefT val(UbDuo<A, B>& d)
{
   return DuoValue<N, UbDuo<A, B> >::get(d);
}

// return Nth value of constant duo
template <int N, typename A, typename B> 
inline typename UbTypeOp<typename UbDuoT<N, UbDuo<A, B> >::ResultT>::RefConstT val(UbDuo<A, B> const& d)
{
   return DuoValue<N, UbDuo<A, B> >::get(d);
}

//duo6.hpp
// partial specialization for UbDuo<> with only one field
template <typename A>
struct UbDuo<A,void> 
{
public:
   typedef A    T1;  // type of first field
   typedef void T2;  // type of second field
   enum { N = 1 };   // number of fields

private:
   T1 value1;        // value of first field

public:
   // constructors
   UbDuo() : value1() { }
   UbDuo (T1 const & a) : value1(a) { }

   // field access
   T1& v1()             { return value1; }
   T1 const& v1() const { return value1; }

   void v2() { }
   void v2() const { }

   #ifdef CAB_RCF
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
         ar & value1;
      }
   #endif
};

//tupel1.hpp
// type that represents unused type parameters
class UbNullT 
{
};

// UbTuple<> in general derives from UbTuple<> with one more UbNullT
template <typename P1,
          typename P2 = UbNullT,
          typename P3 = UbNullT,
          typename P4 = UbNullT,
          typename P5 = UbNullT,
          typename P6 = UbNullT,
          typename P7 = UbNullT,
          typename P8 = UbNullT >
class UbTuple : public UbDuo<P1, typename UbTuple<P2,P3,P4,P5,P6,P7,P8,UbNullT>::BaseT> 
{
public:
   typedef UbDuo<P1, typename UbTuple<P2,P3,P4,P5,P6,P7,P8,UbNullT>::BaseT>  BaseT;

   // constructor:
   UbTuple() {}
   UbTuple( typename UbTypeOp<P1>::RefConstT a1,
            typename UbTypeOp<P2>::RefConstT a2,
            typename UbTypeOp<P3>::RefConstT a3 = UbNullT(),
            typename UbTypeOp<P4>::RefConstT a4 = UbNullT(),
            typename UbTypeOp<P5>::RefConstT a5 = UbNullT(),
            typename UbTypeOp<P6>::RefConstT a6 = UbNullT(),
            typename UbTypeOp<P7>::RefConstT a7 = UbNullT(),
            typename UbTypeOp<P8>::RefConstT a8 = UbNullT() )
      : BaseT(a1, UbTuple<P2,P3,P4,P5,P6,P7,P8,UbNullT>(a2,a3,a4,a5,a6,a7,a8))
   {
   }

   // for implicit type conversion during assignments
   template <typename U1,typename U2, typename U3, typename U4, typename U5, typename U6, typename U7, typename U8 >
   UbTuple<P1,P2,P3,P4,P5,P6,P7,P8>& operator = ( const UbTuple<U1,U2,U3,U4,U5,U6,U7,U8>& rhs)
   {
      this->BaseT::operator=( typename UbTuple<U1,U2,U3,U4,U5,U6,U7,U8>::BaseT(rhs) );
      return *this;
   }

};

// specialization to end deriving recursion
template <typename P1, typename P2>
class UbTuple<P1,P2,UbNullT,UbNullT,UbNullT,UbNullT,UbNullT,UbNullT> : public UbDuo<P1,P2> {
public:
   typedef UbDuo<P1,P2> BaseT;
   
   // constructor:
   UbTuple() {}
   UbTuple( typename UbTypeOp<P1>::RefConstT a1,
            typename UbTypeOp<P2>::RefConstT a2,
            typename UbTypeOp<UbNullT>::RefConstT = UbNullT(),
            typename UbTypeOp<UbNullT>::RefConstT = UbNullT(),
            typename UbTypeOp<UbNullT>::RefConstT = UbNullT(),
            typename UbTypeOp<UbNullT>::RefConstT = UbNullT(),
            typename UbTypeOp<UbNullT>::RefConstT = UbNullT(),
            typename UbTypeOp<UbNullT>::RefConstT = UbNullT() )
      : BaseT(a1, a2) 
   {
   }

   // for implicit type conversion during assignments
   template <typename U1,typename U2 >
   UbTuple<P1,P2>& operator = ( const UbTuple<U1,U2>& rhs)
   {
      this->BaseT::operator=( typename UbTuple<U1,U2>::BaseT(rhs) );
      return *this;
   }

};

// specialization for singletons
template <typename P1>
class UbTuple<P1,UbNullT,UbNullT,UbNullT,UbNullT,UbNullT,UbNullT,UbNullT> : public UbDuo<P1,void>
{
public:
   typedef UbDuo<P1,void> BaseT;

   // constructor:
   UbTuple() {}
   UbTuple( typename UbTypeOp<P1>::RefConstT a1,
            typename UbTypeOp<UbNullT>::RefConstT = UbNullT(),
            typename UbTypeOp<UbNullT>::RefConstT = UbNullT(),
            typename UbTypeOp<UbNullT>::RefConstT = UbNullT(),
            typename UbTypeOp<UbNullT>::RefConstT = UbNullT(),
            typename UbTypeOp<UbNullT>::RefConstT = UbNullT(),
            typename UbTypeOp<UbNullT>::RefConstT = UbNullT(),
            typename UbTypeOp<UbNullT>::RefConstT = UbNullT() )
      : BaseT(a1) 
   {
   }

   // for implicit type conversion during assignments
   template <typename U1 >
   UbTuple<P1>& operator = ( const UbTuple<U1>& rhs)
   {
      this->v1() = rhs.v1();
      return *this;
   }

};

// convenience function for 1 argument
template <typename T1> 
inline UbTuple<T1> makeUbTuple(T1 const &a1)
{
   return UbTuple<T1>(a1);
}

// convenience function for 2 arguments
template <typename T1, typename T2>
inline UbTuple<T1,T2> makeUbTuple(T1 const &a1, T2 const &a2)
{
   return UbTuple<T1,T2>(a1,a2);
}

// convenience function for 3 arguments
template <typename T1, typename T2, typename T3>
inline UbTuple<T1,T2,T3> makeUbTuple(T1 const &a1, T2 const &a2, T3 const &a3)
{
   return UbTuple<T1,T2,T3>(a1,a2,a3);
}

// convenience function for 4 arguments
template <typename T1, typename T2, typename T3, typename T4>
inline UbTuple<T1,T2,T3,T4> makeUbTuple(T1 const &a1, T2 const &a2, T3 const &a3, T4 const &a4)
{
   return UbTuple<T1,T2,T3,T4>(a1,a2,a3,a4);
}

// convenience function for 5 arguments
template <typename T1, typename T2, typename T3, typename T4, typename T5>
inline UbTuple<T1,T2,T3,T4,T5> makeUbTuple(T1 const &a1, T2 const &a2, T3 const &a3, T4 const &a4,T5 const &a5)
{
   return UbTuple<T1,T2,T3,T4,T5>(a1,a2,a3,a4,a5);
}

// convenience function for 6 arguments
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
inline UbTuple<T1,T2,T3,T4,T5,T6> makeUbTuple(T1 const &a1, T2 const &a2, T3 const &a3, T4 const &a4, T5 const &a5, T6 const &a6)
{
   return UbTuple<T1,T2,T3,T4,T5,T6>(a1,a2,a3,a4,a5,a6);
}

// convenience function for 7 arguments
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
inline UbTuple<T1,T2,T3,T4,T5,T6,T7> makeUbTuple(T1 const &a1, T2 const &a2, T3 const &a3, T4 const &a4, T5 const &a5, T6 const &a6, T7 const &a7)
{
   return UbTuple<T1,T2,T3,T4,T5,T6,T7>(a1,a2,a3,a4,a5,a6,a7);
}

// convenience function for 8 arguments
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
inline UbTuple<T1,T2,T3,T4,T5,T6,T7,T8> makeUbTuple(T1 const &a1, T2 const &a2,T3 const &a3, T4 const &a4,T5 const &a5, T6 const &a6,T7 const &a7, T8 const &a8 )
{
   return UbTuple<T1,T2,T3,T4,T5,T6,T7,T8>(a1,a2,a3,a4,a5,a6,a7,a8);
}

//some typedefs
typedef UbTuple<float,float>                               UbTupleFloat2;
typedef UbTuple<float,float,float>                         UbTupleFloat3;
typedef UbTuple<int,int>                                   UbTupleInt2;
typedef UbTuple<int,int,int>                               UbTupleInt3;
typedef UbTuple<int,int,int,int>                           UbTupleInt4;
typedef UbTuple<int,int,int,int,int>                       UbTupleInt5;
typedef UbTuple<int,int,int,int,int,int>                   UbTupleInt6;
typedef UbTuple<int,int,int,int,int,int,int,int>           UbTupleInt8;
typedef UbTuple<double,double>                             UbTupleDouble2;
typedef UbTuple<double,double,double>                      UbTupleDouble3;
typedef UbTuple<double,double,double,double>               UbTupleDouble4;
typedef UbTuple<double,double,double,double,double,double> UbTupleDouble6;
typedef UbTuple<std::string,double,double>                 UbTupleStringDouble2;
typedef UbTuple<std::string,double,double,double>          UbTupleStringDouble3;
typedef UbTuple<std::string,int,int,int>                   UbTupleStringInt3;
typedef UbTuple<short,short,short,short>                   UbTupleShort4;
typedef UbTuple<bool,bool,bool>                            UbTupleBool3;
typedef UbTuple<int,double,double>                         UbTupleIntDouble2;
typedef UbTuple<int, bool>                                 UbTupleIntBool;


#endif //UBTUPLE_H
