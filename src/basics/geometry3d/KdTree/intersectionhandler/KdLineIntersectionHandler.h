#ifndef KDLINEINTERSECTIONHANDLER_H
#define KDLINEINTERSECTIONHANDLER_H

#include <basics/utilities/UbTuple.h>
//#include <geometry3d/GbTriFaceMesh3D.h>

#include <set>

// #ifdef CAB_RCF
// #  include <3rdParty/rcf/RcfSerializationIncludes.h>
// #end
namespace Kd
{
template <typename T>
class Node;

template <typename T>
class LineIntersectionHandler
{
public:
    virtual bool intersectLine(const UbTuple<T, T, T> &n1, const UbTuple<T, T, T> &n2, Node<T> &parent,
                               Node<T> *&child1, Node<T> *&child2) const = 0;
    virtual ~LineIntersectionHandler()                                   = default;
};
} // namespace Kd

// #if defined(RCF_USE_SF_SERIALIZATION) && !defined(SWIG)
//    SF_NO_CTOR(Kd::LineIntersectionHandler<float>);
//    SF_NO_CTOR(Kd::LineIntersectionHandler<double>);
// #endif //RCF_USE_SF_SERIALIZATI
#endif // KDLINEINTERSECTIONHANDLER_H
