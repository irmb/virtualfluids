#ifndef BoundingBoxMemento_H
#define BoundingBoxMemento_H


#include "GridGenerator/global.h"

#include <memory>
#include <string>
#include <vector>

#ifndef __CUDACC__
#include <boost/serialization/access.hpp>
#endif


class VF_PUBLIC BoundingBoxMemento
{
#ifndef __CUDACC__
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & minX;
        ar & maxX;
        ar & minY;
        ar & maxY;
        ar & minZ;
        ar & maxZ;
    }
#endif

public:
    real minX;
    real maxX;
    real minY;
    real maxY;
    real minZ;
    real maxZ;

};


#endif

