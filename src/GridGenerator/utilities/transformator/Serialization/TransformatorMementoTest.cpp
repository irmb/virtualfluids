//#include "gmock/gmock.h"
//
//#include <memory>
//#include <fstream>
//
//#include <GridGenerator/geometries/BoundingBox/Serialization/BoundingBoxMemento.h>
//#include <GridGenerator/geometries/BoundingBox/BoundingBox.cuh>
//
//#include <boost/archive/text_iarchive.hpp>
//#include <boost/archive/text_oarchive.hpp>
//
//void serialize(const BoundingBoxMemento &memento, const std::string &filename)
//{
//    std::ofstream ofs(filename);
//    boost::archive::text_oarchive oa(ofs);
//    oa << memento;
//}
//void deserialize(BoundingBoxMemento &memento, const std::string &filename)
//{
//    std::ifstream ifs(filename);
//    boost::archive::text_iarchive ia(ifs);
//    ia >> memento;
//}
//
//TEST(BoundingBoxMementoTest, serializeAndDeserializeBB)
//{
//    std::string fileName = "boundingbox";
//
//    BoundingBox<real> box(1.2, 22.2, -23.2, 2, 0.0001, 1212122.1);
//
//
//    BoundingBoxMemento sut = box.getState();
//    serialize(sut, fileName);
//
//    BoundingBoxMemento newSut;
//    deserialize(newSut, fileName);
//    BoundingBox<real> boxNew;
//    boxNew.setState(newSut);
//
//    EXPECT_TRUE(box == boxNew);
//}
