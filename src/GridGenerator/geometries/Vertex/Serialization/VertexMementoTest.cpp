#include "gmock/gmock.h"

#include <memory>
#include <fstream>

#include <GridGenerator/geometries/Vertex/Serialization/VertexMemento.h>
#include <GridGenerator/geometries/Vertex/Vertex.cuh>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

void serialize(const VertexMemento &memento, const std::string &filename)
{

    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    oa << memento;

}
void deserialize(VertexMemento &memento, const std::string &filename)
{

    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    ia >> memento;

}

TEST(MementoVertexTest, serializeAndDeserializeVertex)
{
    std::string fileName = "vertex";
    Vertex vOld = Vertex(1, 2, 3);
    VertexMemento sut = vOld.getState();
    serialize(sut, fileName);

    VertexMemento newSut;
    deserialize(newSut, fileName);
    Vertex vNew;
    vNew.setState(newSut);

    EXPECT_TRUE(vOld == vNew);
}
