#include "gmock/gmock.h"

#include "D3Q27System.h"



TEST(D3Q27SystemTest, expectThatNormOfDirectionsIsOne)
{
    for (int fDir = D3Q27System::FSTARTDIR; fDir <= D3Q27System::FENDDIR; fDir++)
    {
        double norm = std::sqrt(std::pow(D3Q27System::cNorm[0][fDir], 2) + std::pow(D3Q27System::cNorm[1][fDir], 2) + std::pow(D3Q27System::cNorm[2][fDir], 2));
        EXPECT_THAT(1.0, testing::DoubleEq(norm));
    }
}
