#include <gmock/gmock.h>

#include "communication.h"

using namespace CommunicationDirections;

TEST(communicationTest, isNegative)
{
    EXPECT_TRUE(isNegative(CommunicationDirection::MX));
    EXPECT_TRUE(isNegative(CommunicationDirection::MY));
    EXPECT_TRUE(isNegative(CommunicationDirection::MZ));
    EXPECT_FALSE(isNegative(CommunicationDirection::PX));
    EXPECT_FALSE(isNegative(CommunicationDirection::PY));
    EXPECT_FALSE(isNegative(CommunicationDirection::PZ));
}

TEST(communicationTest, isPositive)
{
    EXPECT_TRUE(isPositive(CommunicationDirection::PX));
    EXPECT_TRUE(isPositive(CommunicationDirection::PY));
    EXPECT_TRUE(isPositive(CommunicationDirection::PZ));
    EXPECT_FALSE(isPositive(CommunicationDirection::MX));
    EXPECT_FALSE(isPositive(CommunicationDirection::MY));
    EXPECT_FALSE(isPositive(CommunicationDirection::MZ));
}

TEST(communicationTest, getNegativeDirectionAlongAxis)
{
    EXPECT_THAT(getNegativeDirectionAlongAxis(Axis::x), CommunicationDirection::MX);
    EXPECT_THAT(getNegativeDirectionAlongAxis(Axis::y), CommunicationDirection::MY);
    EXPECT_THAT(getNegativeDirectionAlongAxis(Axis::z), CommunicationDirection::MZ);
}

TEST(communicationTest, getPositiveDirectionAlongAxis)
{
    EXPECT_THAT(getPositiveDirectionAlongAxis(Axis::x), CommunicationDirection::PX);
    EXPECT_THAT(getPositiveDirectionAlongAxis(Axis::y), CommunicationDirection::PY);
    EXPECT_THAT(getPositiveDirectionAlongAxis(Axis::z), CommunicationDirection::PZ);
}