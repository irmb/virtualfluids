#include "gmock/gmock.h"

#include <memory>
#include <fstream>

#include <GridGenerator/geometries/Geometry/Serialization/GeometryMemento.h>
#include <GridGenerator/geometries/Geometry/Geometry.cuh>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

void serialize(const GeometryMemento &memento, const std::string &filename)
{

    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    oa << memento;

}
void deserialize(GeometryMemento &memento, const std::string &filename)
{

    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    ia >> memento;

}

TEST(GeometryMementoTest, serializeAndDeserializeGeoemtry)
{
    std::string fileName = "geometry";
    Triangle t1 = Triangle(Vertex(0, 0, 0), Vertex(10, 0, 0), Vertex(0, 10, 0), Vertex(0, 0, 1));
    t1.alphaAngles[0] = 10;
    t1.alphaAngles[1] = 20;
    t1.alphaAngles[2] = 30;

    Triangle t2 = Triangle(Vertex(0, -12, 2), Vertex(10, 12, 0), Vertex(0, 10, -10), Vertex(1, 0, 0));
    t2.alphaAngles[0] = 40;
    t2.alphaAngles[1] = 50;
    t2.alphaAngles[2] = 60;

    Geometry geo;
    geo.triangleVec = { t1, t2 };
    geo.size = 2;

    GeometryMemento sut = geo.getState();
    serialize(sut, fileName);

    GeometryMemento newSut;
    deserialize(newSut, fileName);
    Geometry geoNew;
    geoNew.setState(newSut);

    EXPECT_TRUE(geo == geoNew);
}
