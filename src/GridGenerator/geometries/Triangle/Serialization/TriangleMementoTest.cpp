#include "gmock/gmock.h"

#include <memory>
#include <fstream>

#include <GridGenerator/geometries/Triangle/Serialization/TriangleMemento.h>
#include <GridGenerator/geometries/Triangle/Triangle.cuh>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

void serialize(const TriangleMemento &memento, const std::string &filename)
{

    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    oa << memento;

}
void deserialize(TriangleMemento &memento, const std::string &filename)
{

    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    ia >> memento;

}

TEST(TriangleVertexTest, serializeAndDeserializeTriangle)
{
    std::string fileName = "triangle";
    Triangle tOld = Triangle(Vertex(0, 0, 0), Vertex(10, 0, 0), Vertex(0, 10, 0), Vertex(0, 0, 1));
    tOld.alphaAngles[0] = 10;
    tOld.alphaAngles[1] = 20;
    tOld.alphaAngles[2] = 30;

    TriangleMemento sut = tOld.getState();
    serialize(sut, fileName);

    TriangleMemento newSut;
    deserialize(newSut, fileName);
    Triangle tNew;
    tNew.setState(newSut);

    EXPECT_TRUE(tOld == tNew);
}
