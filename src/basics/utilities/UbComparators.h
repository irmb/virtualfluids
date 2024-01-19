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
#ifndef UBCOMPARATORS_H
#define UBCOMPARATORS_H

#include <functional>

namespace UbComparators
{
// type_traits
template <typename T>
struct MemberInfo; // not defined for correct compiler errors!

// specialization for MemberFunctionsPtr
// C - class with return T method
template <typename T, typename C>
struct MemberInfo<T C::*> {
    using type       = T;
    using class_type = C;

    static T &apply(C &c, T C::*ptr) { return c.*ptr; }
    static const T &apply(const C &c, T C::*ptr) { return c.*ptr; }
};
// specialization for MemberFunctionsPtr
// C - class with return T method
template <typename T, typename C>
struct MemberInfo<T (C::*)()> {
    using type       = T;
    using class_type = C;

    static T apply(C &c, T (C::*ptr)()) { return (c.*ptr)(); }
};
// specialization for const MemberFunctionsPtr
// C - class with return T method
template <typename T, typename C>
struct MemberInfo<T (C::*)() const> {
    using type       = T;
    using class_type = C;

    static T apply(const C &c, T (C::*ptr)() const) { return (c.*ptr)(); }
};

// MemberComparative-Class
template <typename Ptr, typename Comp = std::less<typename MemberInfo<Ptr>::type>>
class MemComp : private Comp // -> usage of Empty Base Class Optimization (EBCO)
{
    using C = typename MemberInfo<Ptr>::class_type;

public:
    MemComp(Ptr ptr, Comp c = Comp()) : Comp(c), mp_(ptr) {}

    bool operator()(C &lhs, C &rhs)
    {
        return Comp::operator()(MemberInfo<Ptr>::apply(lhs, mp_), MemberInfo<Ptr>::apply(rhs, mp_));
    }
    bool operator()(C &lhs, C &rhs) const
    {
        return Comp::operator()(MemberInfo<Ptr>::apply(lhs, mp_), MemberInfo<Ptr>::apply(rhs, mp_));
    }
    bool operator()(const C &lhs, const C &rhs)
    {
        return Comp::operator()(MemberInfo<Ptr>::apply(lhs, mp_), MemberInfo<Ptr>::apply(rhs, mp_));
    }
    bool operator()(const C &lhs, const C &rhs) const
    {
        return Comp::operator()(MemberInfo<Ptr>::apply(lhs, mp_), MemberInfo<Ptr>::apply(rhs, mp_));
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

template <typename Comp, typename Ptr>
MemComp<Ptr, Comp> membercomp(Ptr p, Comp c = Comp())
{
    return MemComp<Ptr, Comp>(p, c);
}

template <template <typename> class Comp, typename Ptr>
MemComp<Ptr, Comp<typename MemberInfo<Ptr>::type>>
membercomp(Ptr p, Comp<typename MemberInfo<Ptr>::type> c = Comp<typename MemberInfo<Ptr>::type>())
{
    return MemComp<Ptr, Comp<typename MemberInfo<Ptr>::type>>(p, c);
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// andere Variante (alerdings ist hier keine Deduction moeglich!!!)
//////////////////////////////////////////////////////////////////////////
// Vergleichs-Templates:
// Funktor zum "non-const" Methodenvergleich: liste.sort( compareMethods<Klasse, int, &Klasse::getVal1  );
template <typename K /*Klasse*/, typename M /*MethodenRueckgabeTyp*/,
          M (K::*fct)() /*MethodenPointer*/> // Allgemeiner Fall
struct compareMethods {
    bool operator()(K &r, K &l) const // da fct nicht const ist, kann auch K nicht const sein. das const hinter der
                                      // deklaration besagt dass compareMethods const sein kann
    {
        return (r.*fct)() < (l.*fct)();
    }
};
//////////////////////////////////////////////////////////////////////////
// Funktor zum "const" Methodenvergleich: liste.sort( compareMethods<Klasse, int, &Klasse::getVal1  );
template <typename K /*Klasse*/, typename M /*MethodenRueckgabeTyp*/,
          M (K::*fct)() const /*MethodenPointer*/> // <- hier const
struct compareConstMethods {
    bool operator()(const K &r, const K &l)
        const // hier koennen die K's auch const sein, muessen sie aber nicht (const hinzufuegen geht ja problemlos)
    {
        return (r.*fct)() < (l.*fct)();
    }
};
//////////////////////////////////////////////////////////////////////////
// Funktor zum Membervergleich: lise.sort( compareMember<Klasse, int, &Klasse::member>() );
template <typename K /*Klasse*/, typename M /*MemberTyp*/, M(K::*Member) /*MemberPointer*/> // <- hier const
struct compareMember {
    bool operator()(const K &r, const K &l) const { return r.*Member < l.*Member; }
};
// Bsp:
// class Klasse{
// public:
//   Klasse(double val1, double val2 ) : val1(val1),val2(val2) {}
//   double getVal1()       { return val1; }
//   double getVal2() const { return val2; } // <- hier const
//   double val1, val2;
//};
// int main(int argc, char** argv){
//   std::list<Klasse> l;
//   l.push_back( Klasse(10,10) );
//   l.push_back( Klasse(1,5)   );
//   l.sort( compareMember<Klasse, double,  &Klasse::val1 >() );
//   l.sort( compareMethods<Klasse, double,  &Klasse::getVal1 >() );
//   l.sort( compareConstMethods<Klasse, double,  &Klasse::getVal1 >() );
//}

} // namespace UbComparators

#endif // UBCOMPARATOR_H

// example
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

//! \}
