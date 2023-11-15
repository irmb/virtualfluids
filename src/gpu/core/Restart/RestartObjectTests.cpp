#include <gmock/gmock.h>

#include "RestartObject.h"


template <typename Type>
void saveAndLoad()
{
    std::shared_ptr<RestartObject> write_object = std::make_shared<Type>();

    write_object->fs = std::vector<std::vector<real>> {
                { 1,2,3 },
                { 4,5,6 }
            };

    const std::string name {"test_do_check_point"};
    write_object->serialize_internal(name);


    std::shared_ptr<RestartObject> read_object = std::make_shared<Type>();
    read_object->deserialize_internal(name);

    EXPECT_THAT(write_object->fs, ::testing::ContainerEq(read_object->fs));
}

TEST(RestartObjectTests, saveAndLoad_ascii)
{
    saveAndLoad<ASCIIRestartObject>();
}


TEST(RestartObjectTests, saveAndLoad_binary)
{
    saveAndLoad<BinaryRestartObject>();
}
