#ifndef KDRAYINTERSECTIONHANDLER_H
#define KDRAYINTERSECTIONHANDLER_H

#include <basics/utilities/UbTuple.h>
#include <geometry3d/KdTree/KdRay.h>

#include <set>

namespace Kd
{
template <typename T>
class RayIntersectionHandler
{
public:
    virtual int intersectRay(const Ray<T> &ray, Node<T> &parent, Node<T> *&child1, Node<T> *&child2,
                             std::set<UbKeys::Key3<int>> &mailbox) const = 0;
    virtual ~RayIntersectionHandler()                                    = default;
};
} // namespace Kd

#endif // KDRAYINTERSECTIONHANDLER_H
