#include "A.h"

A::~A()
{

}

std::shared_ptr<A> A::make() const
{
    return std::shared_ptr<A>(new A());
}

A::A()
{

}

A::A(const A& a)
{

}
