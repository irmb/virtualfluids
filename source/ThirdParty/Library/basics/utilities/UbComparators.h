#ifndef UBCOMPARATORS_H 
#define UBCOMPARATORS_H

#include <functional> 

/*=========================================================================*/
/*  UbComparators                                                             */
/*                                                                         */
/**
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 16.08.2007
*/ 

namespace UbComparators
{
   //type_traits 
   template <typename T> struct MemberInfo; //not defined for correct compiler errors!

   // specialization for MemberFunctionsPtr
   // C - class with return T method
   template <typename T, typename C> 
   struct MemberInfo<T C::*> 
   { 
      typedef T type; 
      typedef C class_type; 

      static       T& apply(       C& c, T C::* ptr ) { return c.*ptr; } 
      static const T& apply( const C& c, T C::* ptr ) { return c.*ptr; } 
   }; 
   //specialization for MemberFunctionsPtr
   //C - class with return T method
   template <typename T, typename C> 
   struct MemberInfo<T (C::*)()> 
   { 
      typedef T type; 
      typedef C class_type; 

      static T apply( C& c, T (C::*ptr)() ) { return (c.*ptr)(); } 
   }; 
   //specialization for const MemberFunctionsPtr
   //C - class with return T method
   template <typename T, typename C> 
   struct MemberInfo<T (C::*)() const> 
   { 
      typedef T type; 
      typedef C class_type; 

      static T apply( const C& c, T (C::*ptr)() const ) { return (c.*ptr)(); } 
   }; 

   //MemberComparative-Class
   template <typename Ptr, typename Comp = std::less<typename MemberInfo<Ptr>::type> > 
   class MemComp 
      : private Comp  // -> usage of Empty Base Class Optimization (EBCO) 
   { 
      typedef typename MemberInfo<Ptr>::class_type C; 

   public: 
      MemComp( Ptr ptr, Comp c = Comp() ) 
         : Comp(c), mp_(ptr) 
      {} 

      bool operator()(C& lhs, C& rhs) 
      { 
         return Comp::operator()( MemberInfo<Ptr>::apply(lhs, mp_), MemberInfo<Ptr>::apply(rhs, mp_) ); 
      } 
      bool operator()(C& lhs, C& rhs) const 
      { 
         return Comp::operator()( MemberInfo<Ptr>::apply(lhs, mp_), MemberInfo<Ptr>::apply(rhs, mp_) ); 
      } 
      bool operator()(const C& lhs, const C& rhs) 
      { 
         return Comp::operator()( MemberInfo<Ptr>::apply(lhs, mp_), MemberInfo<Ptr>::apply(rhs, mp_) ); 
      } 
      bool operator()(const C& lhs, const C& rhs) const 
      { 
         return Comp::operator()( MemberInfo<Ptr>::apply(lhs, mp_), MemberInfo<Ptr>::apply(rhs, mp_) ); 
      } 

   private: 
      Ptr mp_; 
   }; 

   // Factoryfunktionen 
   template <typename Ptr> 
   MemComp<Ptr> membercomp(Ptr p) 
   { 
      return MemComp<Ptr>(p); 
   } 

   template<typename Comp, typename Ptr> 
   MemComp<Ptr, Comp> membercomp(Ptr p, Comp c = Comp()) 
   { 
      return MemComp<Ptr, Comp>(p, c); 
   } 

   template<template<typename> class Comp, typename Ptr> 
   MemComp<Ptr, Comp<typename MemberInfo<Ptr>::type> > 
      membercomp(Ptr p, Comp<typename MemberInfo<Ptr>::type> c = Comp<typename MemberInfo<Ptr>::type>()) 
   {
      return MemComp<Ptr, Comp<typename MemberInfo<Ptr>::type> >(p, c); 
   } 

    
   //////////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////
   //andere Variante (alerdings ist hier keine Deduction moeglich!!!)
   //////////////////////////////////////////////////////////////////////////
   //Vergleichs-Templates:
   //Funktor zum "non-const" Methodenvergleich: liste.sort( compareMethods<Klasse, int, &Klasse::getVal1  );
   template<typename K/*Klasse*/, typename M /*MethodenRueckgabeTyp*/, M (K::*fct)() /*MethodenPointer*/> // Allgemeiner Fall
   struct compareMethods
   {
      bool operator()(K& r, K& l) const // da fct nicht const ist, kann auch K nicht const sein. das const hinter der deklaration besagt dass compareMethods const sein kann
      { return (r.*fct)() < (l.*fct)();  }
   };
   //////////////////////////////////////////////////////////////////////////
   //Funktor zum "const" Methodenvergleich: liste.sort( compareMethods<Klasse, int, &Klasse::getVal1  );
   template<typename K/*Klasse*/, typename M /*MethodenRueckgabeTyp*/, M (K::*fct)() const /*MethodenPointer*/> // <- hier const 
   struct compareConstMethods
   {
      bool operator()(const K& r, const K& l) const //hier koennen die K's auch const sein, muessen sie aber nicht (const hinzufuegen geht ja problemlos)
      { return (r.*fct)() < (l.*fct)();  }
   };
   //////////////////////////////////////////////////////////////////////////
   //Funktor zum Membervergleich: lise.sort( compareMember<Klasse, int, &Klasse::member>() );
   template<typename K/*Klasse*/, typename M /*MemberTyp*/, M (K::*Member) /*MemberPointer*/> // <- hier const 
   struct compareMember
   { 
      bool operator()(const K& r,const K& l) const
      { return r.*Member < l.*Member; } 
   };
   //Bsp:
   //class Klasse{ 
   //public: 
   //   Klasse(double val1, double val2 ) : val1(val1),val2(val2) {} 
   //   double getVal1()       { return val1; } 
   //   double getVal2() const { return val2; } // <- hier const
   //   double val1, val2; 
   //}; 
   //int main(int argc, char** argv){ 
   //   std::list<Klasse> l; 
   //   l.push_back( Klasse(10,10) ); 
   //   l.push_back( Klasse(1,5)   ); 
   //   l.sort( compareMember<Klasse, double,  &Klasse::val1 >() ); 
   //   l.sort( compareMethods<Klasse, double,  &Klasse::getVal1 >() ); 
   //   l.sort( compareConstMethods<Klasse, double,  &Klasse::getVal1 >() ); 
   //} 

};

#endif //UBCOMPARATOR_H

//example
// #include <basics/utilities/UbComparators.h" 
// #include <list> 
// using namespace std; 
// using namespace UbComparators; 
// 
// struct S { 
//    S(int i) :x(i) {} 
//    int x; 
//    float f() {return x;}; 
//    double g() const {return x;} 
// }; 
// 
// struct intComp { 
//    bool operator()(int l, int r) const 
//    { return l > r; } 
// }; 
// 
// struct dblComp { 
//    bool operator()(double l,  double r) const 
//    { return l > r; } 
// }; 
// 
// template <typename T> 
// struct genComp { 
//    bool operator()(const T& l, const T& r) const
//    { return l > r; } 
// }; 
// 
// 
// int main() 
// { 
//    S a(1); 
//    S b(2); 
//    list<S> sList; 
//    sList.push_back(a); 
//    sList.push_back(b); 
//    sList.sort(UbComparators::membercomp(&S::x,intComp()));  //calls overload (1) 
//    sList.sort(UbComparators::membercomp<intComp>(&S::x));   //same 
//    sList.sort(UbComparators::membercomp(&S::x));            //calls overload (5) 
//    sList.sort(UbComparators::membercomp<genComp>(&S::x));   //calls overload(3) 
//    sList.sort(UbComparators::membercomp(&S::x, genComp<int>())); //calls overload(1) 
//    //same for nonconst function 
//    sList.sort(UbComparators::membercomp(&S::f, dblComp())); //overload(2) 
//    sList.sort(UbComparators::membercomp<dblComp>(&S::f));   //same      
//    sList.sort(UbComparators::membercomp(&S::f));            //overload(6) 
//    sList.sort(UbComparators::membercomp<genComp>(&S::f));   //overload(4) 
//    //same for const function 
//    sList.sort(UbComparators::membercomp(&S::g, dblComp())); //overload(2) 
//    sList.sort(UbComparators::membercomp<dblComp>(&S::g));   //same      
//    sList.sort(UbComparators::membercomp(&S::g));            //overload(6) 
//    sList.sort(UbComparators::membercomp<genComp>(&S::g));   //overload(4) 
// } 
