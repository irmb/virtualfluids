/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef A_H
#define A_H


#include "PointerDefinitions.h"


class VF_PUBLIC A
{
public:
    virtual ~A() = default;
    static SPtr<A> make();

private:
    A() = default;
    A(const A& a) = default;

};

#endif
