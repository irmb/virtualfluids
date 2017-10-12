#ifdef BUILD_JSONCPP

#include "JsonInput.h"
#include "gmock/gmock.h"


using namespace input;


TEST(JsonInputTest, configFileWithKeyValuePair_ShouldReadTheCorrectValuePair)
{
    std::istringstream stringstream(
        "{\n"
        "\"key\": \"value\",\n"
        "\"key2\": \"value2\"\n"
        "}"
        );

    JsonInput input(stringstream);

    EXPECT_THAT(input.getValue("key"), "value");
    EXPECT_THAT(input.getValue("key2"), "value2");
}

TEST(JsonInputTest, configFileWith3DKeyValuePair_ShouldReadTheCorrectValuePair)
{
    std::istringstream stringstream 
        (
        "{\n"
            "\"1Dkey\":\n"
          "{\n"
            "\"2Dkey\":\n"
            "{\n"
              "\"key\": \"value\"\n"
            "}\n"
          "}\n"
        "}\n"
        );

    JsonInput input(stringstream);

    EXPECT_THAT(input.getValue("1Dkey 2Dkey key"), "value");
}

TEST(JsonInputTest, configFileWithKeyValuePair_ShouldReadTheCorrectValuePair_WithIgnoringBlankLines)
{
    std::istringstream stringstream(
        "{\n"
        "\"key\": \"value\",\n"

        "\"key2\": \"value2\"\n"


        "}"
        );

    JsonInput input(stringstream);

    EXPECT_THAT(input.getValue("key"), "value");
    EXPECT_THAT(input.getValue("key2"), "value2");
}

TEST(JsonInputTest, configFileWithKeyValuePair_ShouldReadTheCorrectValuePair_WithIgnoringSpacesAfterAndBeforeColon)
{
    std::istringstream stringstream(
        "{\n"
        "\"key\":              \"value\",\n"

        "\"key2\"      : \"value2\"\n"


        "}"
        );

    JsonInput input(stringstream);

    EXPECT_THAT(input.getValue("key"), "value");
    EXPECT_THAT(input.getValue("key2"), "value2");
}

TEST(JsonInputTest, configFileWithKeyValuePair_getValueFromInvalidKey_ShouldReturnEmptyString)
{
    std::istringstream stringstream(
        ""
        );

    JsonInput input(stringstream);

    EXPECT_THAT(input.getValue("key"), "");
}

TEST(JsonInputTest, configFileWithKeyValuePair_HasKeyWithInvalidKey_ReturnFalse)
{
    std::istringstream stringstream(
        ""
        );

    JsonInput input(stringstream);

    EXPECT_THAT(input.hasValue("key"), false);
}

TEST(JsonInputTest, configFileWithKeyValuePair_HasKeyWithValidKey_ReturnTrue)
{
    std::istringstream stringstream(
        "{\n"
        "\"key\": \"value\",\n"
        "\"key2\": \"value2\"\n"
        "}"
        );

    JsonInput input(stringstream);

    EXPECT_THAT(input.hasValue("key"), true);
}

#endif
