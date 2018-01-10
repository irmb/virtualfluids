#ifndef TriangleMemento_H
#define TriangleMemento_H


#include "GridGenerator/global.h"

#include <memory>
#include <string>

#ifndef __CUDACC__
#include <boost/serialization/access.hpp>
#endif

#include <GridGenerator/geometries/Vertex/Serialization/VertexMemento.h>

class VF_PUBLIC TriangleMemento
{
#ifndef __CUDACC__
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & v1 & v2 & v3 & normal;
        ar & alphaAngles;
    }
#endif

public:
    VertexMemento v1, v2, v3, normal;
    real alphaAngles[3];

};


#endif

