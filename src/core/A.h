/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef A_H
#define A_H

#include <memory>

#include <VirtualFluidsDefinitions.h>

class A;
std::shared_ptr<A> APtr;

class VF_PUBLIC A
{
public:
    virtual ~A();
    std::shared_ptr<A> make() const;

private:
    A();
    A(const A& a);

};

#endif
