#include "gmock/gmock.h"

#include <GridGenerator/geometries/Triangle/Triangle.h>
#include <GridGenerator/geometries/BoundingBox/BoundingBox.h>
#include <GridGenerator/geometries/Vertex/Vertex.h>


using namespace testing;


TEST(BoundingBoxExactTest, findMinMaxFromTriangle)
{
    BoundingBox box = BoundingBox::makeInvalidMinMaxBox();

    real minX = 1.0f;
    real minY = 23.0f;
    real minZ = 1.1222f;

    real maxX = 110.0f;
    real maxY = 50.0f;
    real maxZ = 12122.23f;
    Triangle t = Triangle(Vertex(maxX, maxY - 10, minZ + 2), Vertex(minX, maxY, maxZ), Vertex(minX + 3, minY, minZ), Vertex(0.0f, 0.0f, 0.0f));

	box.setMinMax(t);

	EXPECT_THAT(box.minX, RealEq(minX));
	EXPECT_THAT(box.minY, RealEq(minY));
	EXPECT_THAT(box.minZ, RealEq(minZ));
	
	EXPECT_THAT(box.maxX, RealEq(maxX));
	EXPECT_THAT(box.maxY, RealEq(maxY));
	EXPECT_THAT(box.maxZ, RealEq(maxZ));
}

TEST(BoundingBoxTest, isInside_true)
{
    BoundingBox box = BoundingBox();

    box.minX = 0.0f;
    box.minY = 0.0f;
    box.minZ = 0.0f;

    box.maxX = 10.0f;
    box.maxY = 10.0f;
    box.maxZ = 10.0f;

    Triangle t = Triangle(Vertex(1,1,1), Vertex(2,2,2), Vertex(3,3,3), Vertex(0.0f, 0.0f, 0.0f));

    EXPECT_TRUE(box.isInside(t));
}

TEST(BoundingBoxTest, isInside_false)
{
    BoundingBox box = BoundingBox();

    box.minX = 0.0f;
    box.minY = 0.0f;
    box.minZ = 0.0f;

    box.maxX = 10.0f;
    box.maxY = 10.0f;
    box.maxZ = 10.0f;

    Triangle t = Triangle(Vertex(1, 1, 1), Vertex(2, 2, 2), Vertex(3, 3, 11), Vertex(0.0f, 0.0f, 0.0f));

    EXPECT_FALSE(box.isInside(t));
}


