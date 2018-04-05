#ifndef GeometryMemento_H
#define GeometryMemento_H


#include "GridGenerator/global.h"

#include <memory>
#include <string>
#include <vector>

#ifndef __CUDACC__
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#endif

#include <GridGenerator/geometries/BoundingBox/Serialization/BoundingBoxMemento.h>
#include <GridGenerator/geometries/Triangle/Serialization/TriangleMemento.h>

class VF_PUBLIC GeometryMemento
{
#ifndef __CUDACC__
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & triangles;
        ar & minmax;
    }
#endif

public:
    std::vector<TriangleMemento> triangles;
    BoundingBoxMemento minmax;
};


#endif

