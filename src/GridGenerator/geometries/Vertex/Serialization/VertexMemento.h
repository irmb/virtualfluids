#ifndef VertexSerializer_H
#define VertexSerializer_H


#include "GridGenerator/global.h"

#include <memory>
#include <string>

#ifndef __CUDACC__
#include <boost/serialization/access.hpp>
#endif


class VF_PUBLIC VertexMemento
{
#ifndef __CUDACC__
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & x & y & z;
    }
#endif
public:
    doubflo x, y, z;
};


#endif

