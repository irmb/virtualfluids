#include "gmock/gmock.h"
#include <cmath>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "QDebugVtkWriter.hpp"
#include <tuple>

MATCHER(DoubleNear5, "") {
    return abs(std::get<0>(arg) - std::get<1>(arg)) < 0.00001;
}

using namespace QDebugVtkWriter;

double calcVectorLength(const std::array<double, 3> coords, const std::array<double, 3> neighborCoords)
{
    return std::sqrt(std::pow((neighborCoords[0] - coords[0]), 2) + std::pow((neighborCoords[1] - coords[1]), 2) +
                     std::pow((neighborCoords[2] - coords[2]), 2));
}

TEST(QDebugVtkWriterTest, modifyLineLengthsForQsSameCoords3)
{
    const std::array<double, 3> coords = { 0, 0, 0 };
    std::array<double, 3> neighborCoords = { 1, 1, 1 };
    const real q = 0.3;
    const real initialLength = calcVectorLength(coords, neighborCoords);

    modifyLineLengthsForQs(coords, neighborCoords, q);

    std::array<double, 3> expectedNeighborCoords = { 0.3, 0.3, 0.3 };
    EXPECT_THAT(neighborCoords,testing::Pointwise(DoubleNear5(), expectedNeighborCoords));
    EXPECT_THAT(calcVectorLength(coords, neighborCoords), testing::DoubleNear(q*initialLength, 0.00001));
}

TEST(QDebugVtkWriterTest, modifyLineLengthDifferentCoords)
{
    const std::array<double, 3> coords = { 0, 0, 0 };
    std::array<double, 3> neighborCoords = { 1, 2, 3 };
    const real q = 0.3;
    const real initialLength = calcVectorLength(coords, neighborCoords);

    modifyLineLengthsForQs(coords, neighborCoords, q);

    std::array<double, 3> expectedNeighborCoords = { 0.3, 0.6, 0.9 };
    EXPECT_THAT(neighborCoords,testing::Pointwise(DoubleNear5(), expectedNeighborCoords));
    EXPECT_THAT(calcVectorLength(coords, neighborCoords), testing::DoubleNear(q*initialLength, 0.00001));
}

TEST(QDebugVtkWriterTest, modifyLineLengthNegativeCoord)
{
    const std::array<double, 3> coords = { 0, 0, 0 };
    std::array<double, 3> neighborCoords = { 1, 2, -3 };
    const real q = 0.3;
    const real initialLength = calcVectorLength(coords, neighborCoords);

    modifyLineLengthsForQs(coords, neighborCoords, q);

    std::array<double, 3> expectedNeighborCoords = { 0.3, 0.6, -0.9 };
    EXPECT_THAT(neighborCoords,testing::Pointwise(DoubleNear5(), expectedNeighborCoords));
    EXPECT_THAT(calcVectorLength(coords, neighborCoords), testing::DoubleNear(q*initialLength, 0.00001));
}