#include <gmock/gmock.h>

#include "ConfigurationFile.h"

using namespace vf::basics;

TEST(ConfigurationFileTest, ContainsReturnsTrueForExistingKey)
{
    ConfigurationFile config;
    config.data["key1"] = "value1";
    config.data["key2"] = "value2";

    EXPECT_TRUE(config.contains("key1"));
    EXPECT_TRUE(config.contains("key2"));
}

TEST(ConfigurationFileTest, ContainsReturnsFalseForMissingKey)
{
    ConfigurationFile config;
    config.data["key1"] = "value1";

    EXPECT_FALSE(config.contains("key2"));
}

TEST(ConfigurationFileTest, GetValueReturnsCorrectValue)
{
    ConfigurationFile config;
    config.data["key1"] = "value1";
    config.data["key2"] = "1234";

    EXPECT_EQ(config.getValue<std::string>("key1"), "value1");
    EXPECT_EQ(config.getValue<int>("key2"), 1234);
}

TEST(ConfigurationFileTest, GetValueThrowsExceptionForMissingKey)
{
    ConfigurationFile config;
    config.data["key1"] = "value1";

    EXPECT_THROW(config.getValue<std::string>("key2"), UbException);
}

TEST(ConfigurationFileTest, GetVectorReturnsCorrectValues)
{
    ConfigurationFile config;
    config.data["key1"] = "1, 2, 3";
    config.data["key2"] = "4; 5; 6";

    std::vector<int> v1 = config.getVector<int>("key1");
    std::vector<int> v2 = config.getVector<int>("key2");

    EXPECT_EQ(v1.size(), 3);
    EXPECT_EQ(v1[0], 1);
    EXPECT_EQ(v1[1], 2);
    EXPECT_EQ(v1[2], 3);

    EXPECT_EQ(v2.size(), 3);
    EXPECT_EQ(v2[0], 4);
    EXPECT_EQ(v2[1], 5);
    EXPECT_EQ(v2[2], 6);
}

TEST(ConfigurationFileTest, GetVectorThrowsExceptionForMissingKey)
{
    ConfigurationFile config;
    config.data["key1"] = "1, 2, 3";

    EXPECT_THROW(config.getVector<int>("key2"), UbException);
}

TEST(ConfigurationFileTest, GetValueReturnsDefaultValueForMissingKey)
{
    ConfigurationFile config;
    config.data["key1"] = "value1";

    int defaultValue = 42;
    int value = config.getValue<int>("key2", defaultValue);

    EXPECT_EQ(value, defaultValue);
}

TEST(ConfigurationFileTest, GetValueReturnsCorrectValueForExistingKey)
{
    ConfigurationFile config;
    config.data["key1"] = "42";

    int defaultValue = 0;
    int value = config.getValue<int>("key1", defaultValue);

    EXPECT_EQ(value, 42);
}

TEST(ConfigurationFileTest, FromStringConvertsValueToBool)
{
    ConfigurationFile config;
    config.data["key1"] = "true";
    config.data["key2"] = "false";

    bool valueTrue = config.getValue<bool>("key1");
    bool valueFalse = config.getValue<bool>("key2");

    EXPECT_TRUE(valueTrue);
    EXPECT_FALSE(valueFalse);
}

TEST(ConfigurationFileTest, GetValueThrowsExceptionForWrongTypeConversion)
{
    ConfigurationFile config;
    config.data["key1"] = "text";

    EXPECT_THROW(config.getVector<int>("key1"), UbException);
}

TEST(ConfigurationFileTest, ClearRemovesAllData)
{
    ConfigurationFile config;
    config.data["key1"] = "value1";
    config.data["key2"] = "value2";

    config.clear();

    EXPECT_TRUE(config.data.empty());
}
