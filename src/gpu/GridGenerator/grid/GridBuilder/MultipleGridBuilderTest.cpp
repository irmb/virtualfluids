#include "basics/tests/testUtilities.h"
#include "geometries/Sphere/Sphere.h"
#include "grid/GridBuilder/MultipleGridBuilder.h"

TEST(MultipleGridBuilderTest, noCoarseGrid_addFineGrid_warns)
{
    MultipleGridBuilder gridBuilder;
    SPtr<Object> gridShape;

    testingVF::captureStdOut();
    gridBuilder.addGrid(gridShape, 1);
    EXPECT_TRUE(testingVF::stdoutContainsWarning());
}

TEST(MultipleGridBuilderTest, coarseGridExist_addFineGrid_doesNotWarn)
{
    MultipleGridBuilder gridBuilder;

    SPtr<Object> gridShape = std::make_shared<Sphere>(0.0, 0.0, 0.0, 0.0);
    gridBuilder.addCoarseGrid(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.1);

    testingVF::captureStdOut();
    gridBuilder.addGrid(gridShape, 1);
    EXPECT_FALSE(testingVF::stdoutContainsWarning());
}