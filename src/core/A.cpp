#include "A.h"


SPtr<A> A::make()
{
    return SPtr<A>(new A());
}
