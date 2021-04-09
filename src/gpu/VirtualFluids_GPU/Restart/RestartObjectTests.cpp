#include <gmock/gmock.h>

#include "RestartObject.h"


TEST(RestartObjectTests, saveAndLoad_ascii)
{
    std::shared_ptr<RestartObject> sut = std::make_shared<ASCIIRestartObject>();

    sut->fs = std::vector<std::vector<float>> {
                { 1,2,3 },
                { 4,5,6 }
            };

    std::string name {"test_do_check_point"};
    sut->serialize_internal(name);


    std::shared_ptr<RestartObject> obj_read = std::make_shared<ASCIIRestartObject>();
    obj_read->deserialize_internal(name);

    EXPECT_THAT(sut->fs, ::testing::ContainerEq(obj_read->fs));
}


TEST(RestartObjectTests, saveAndLoad_binary)
{
    std::shared_ptr<RestartObject> sut = std::make_shared<BinaryRestartObject>();

    sut->fs = std::vector<std::vector<float>> {
                { 1,2,3 },
                { 4,5,6 }
            };

    std::string name {"test_do_check_point"};
    sut->serialize_internal(name);


    std::shared_ptr<RestartObject> obj_read = std::make_shared<BinaryRestartObject>();
    obj_read->deserialize_internal(name);

    EXPECT_THAT(sut->fs, ::testing::ContainerEq(obj_read->fs));
}
