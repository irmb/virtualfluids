#include "ConfigInput.h"
#include "gmock/gmock.h"



using namespace input;


TEST(ConfigInputTest, configFileWithKeyValuePair_ShouldReadTheCorrectValuePair)
{
    std::istringstream stringstream(
        "key=value\n"
        "key2=value2\n"
        );

    ConfigInput input(stringstream);

    EXPECT_THAT(input.getValue("key"), "value");
    EXPECT_THAT(input.getValue("key2"), "value2");
}

TEST(ConfigInputTest, configFileWithKeyValuePair_ShouldReadTheCorrectValuePair_WithIgnoringBlankLines)
{
    std::istringstream stringstream(
        "\n\n"
        "key=value\n"
        "\n"
        "key2=value2\n"
        "\n"
        );

    ConfigInput input(stringstream);

    EXPECT_THAT(input.getValue("key"), "value");
    EXPECT_THAT(input.getValue("key2"), "value2");
}

TEST(ConfigInputTest, configFileWithKeyValuePair_ShouldReadTheCorrectValuePair_WithIgnoringSpacesAfterEqual)
{
    std::istringstream stringstream(
        "\n\n"
        "key=  value\n"
        "\n"
        "key2=        value2\n"
        "\n"
        );

    ConfigInput input(stringstream);

    EXPECT_THAT(input.getValue("key"), "value");
    EXPECT_THAT(input.getValue("key2"), "value2");
}


TEST(DISABLED_ConfigInputTest, configFileWithKeyValuePair_ShouldReadTheCorrectValuePair_WithIgnoringSpacesBeforeEqual)
{
    std::istringstream stringstream(
        "\n\n"
        "key   =  value\n"
        "\n"
        "key2         =        value2\n"
        "\n"
        );

    ConfigInput input(stringstream);

    EXPECT_THAT(input.getValue("key"), "value");
    EXPECT_THAT(input.getValue("key2"), "value2");
}

TEST(ConfigInputTest, configFileWithKeyValuePair_ShouldReadTheCorrectValuePair_WithMissingValue)
{
    std::istringstream stringstream(
        "\n\n"
        "key=  \n"
        "\n"
        "key2=        value2\n"
        "\n"
        );

    ConfigInput input(stringstream);

    EXPECT_THAT(input.getValue("key"), "");
    EXPECT_THAT(input.getValue("key2"), "value2");
}

TEST(ConfigInputTest, configFileWithKeyValuePair_getValueFromInvalidKey_ShouldReturnEmptyString)
{
    std::istringstream stringstream(
       ""
        );

    ConfigInput input(stringstream);

    EXPECT_THAT(input.getValue("key"), "");
}

TEST(ConfigInputTest, configFileWithKeyValuePair_HasKeyWithInvalidKey_ReturnFalse)
{
    std::istringstream stringstream(
        ""
        );

    ConfigInput input(stringstream);

    EXPECT_THAT(input.hasValue("key"), false);
}

TEST(ConfigInputTest, configFileWithKeyValuePair_HasKeyWithValidKey_ReturnTrue)
{
    std::istringstream stringstream(
        "key=value\n"
        );

    ConfigInput input(stringstream);

    EXPECT_THAT(input.hasValue("key"), true);
}

TEST(ConfigInputTest, configFileWithKeyValuePair_ReturnValueAndIgnoreComments)
{
    std::istringstream stringstream(
        "# comment\n"
        "key=value\n"
        "# comment\n"
        );

    ConfigInput input(stringstream);

    EXPECT_THAT(input.getValue("key"), "value");
}

TEST(ConfigInputTest, configFileWithKeyValuePair_ReturnValue_ShouldIgnoreQuotationMarks)
{
    std::istringstream stringstream(
        "# comment\n"
        "key=\"value\"\n"
        "# comment\n"
        );

    ConfigInput input(stringstream);

    EXPECT_THAT(input.getValue("key"), "value");
}
