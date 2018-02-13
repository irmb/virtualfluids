#include "gmock/gmock.h"

#include <GridGenerator/geometries/Triangle/Triangle.cuh>
#include <GridGenerator/geometries/BoundingBox/BoundingBox.cuh>
#include <GridGenerator/geometries/Vertex/Vertex.cuh>


using namespace testing;


TEST(BoundingBoxTest, initWithTriangle_whenTheValueIsIntegerBoxHasToPLUS_or_MINUS_ONE) 
{
    real minX = 1.0f;
    real minY = 23.0f;
    real minZ = 1.1222f;

    real maxX = 110.0f;
    real maxY = 50.0f;
    real maxZ = 12122.23f;
	BoundingBox<int> box = BoundingBox<int>::makeNodeBox(Triangle(Vertex(maxX, maxY - 10, minZ + 2), Vertex(minX, maxY, maxZ), Vertex(minX + 3, minY, minZ), Vertex(0.0f, 0.0f, 0.0f)));
    EXPECT_THAT(box.minX, Eq(minX - 1));
    EXPECT_THAT(box.minY, Eq(minY - 1));
    EXPECT_THAT(box.minZ, Eq((int)minZ));

    EXPECT_THAT(box.maxX, Eq(maxX + 1));
    EXPECT_THAT(box.maxY, Eq(maxY + 1));
    EXPECT_THAT(box.maxZ, Eq((int)maxZ + 1));
}

TEST(BoundingBoxTest, initWithTriangle2)
{
	BoundingBox<int> box = BoundingBox<int>::makeNodeBox(Triangle(Vertex(20.0f, 1.0f, 1.0f), Vertex(1.0f, 1.0f, 1 + 1e-006f), Vertex(20.0f, 20.0f, 1.0f), Vertex(1.0f, 0.0f, 0.0f)));

    EXPECT_THAT(box.minX, Eq(0));
    EXPECT_THAT(box.minY, Eq(0));
    EXPECT_THAT(box.minZ, Eq(0));

    EXPECT_THAT(box.maxX, Eq(21));
    EXPECT_THAT(box.maxY, Eq(21));
    EXPECT_THAT(box.maxZ, Eq(2));
}


TEST(BoundingBoxTest, initWithTriangle3) 
{
	BoundingBox<int> box = BoundingBox<int>::makeNodeBox(Triangle(Vertex(20.0f, 20.0f, 20.0f), Vertex(1.0f, 20.0f, 20.0f), Vertex(20.0f, 1.0f, 20.0f), Vertex(1.0f, 0.0f, 0.0f)));

    EXPECT_THAT(box.minX, Eq(0));
    EXPECT_THAT(box.minY, Eq(0));
    EXPECT_THAT(box.minZ, Eq(19));

    EXPECT_THAT(box.maxX, Eq(21));
    EXPECT_THAT(box.maxY, Eq(21));
    EXPECT_THAT(box.maxZ, Eq(21));
}

TEST(BoundingBoxTest, whenAllValueAreFloat_BoxHasToCEIL_OR_FLOOR) 
{
    real minX = 1.5f;
    real minY = 23.2f;
    real minZ = 1.1222f;

    real maxX = 110.4f;
    real maxY = 50.5f;
    real maxZ = 12122.23f;

	BoundingBox<int> box = BoundingBox<int>::makeNodeBox(Triangle(Vertex(maxX, maxY - 10, minZ + 2), Vertex(minX, maxY, maxZ), Vertex(minX + 3, minY, minZ), Vertex(0.0f, 0.0f, 0.0f)));

    EXPECT_THAT(box.minX, Eq((int)minX));
    EXPECT_THAT(box.minY, Eq((int)minY));
    EXPECT_THAT(box.minZ, Eq((int)minZ));

    EXPECT_THAT(box.maxX, Eq((int)maxX + 1));
    EXPECT_THAT(box.maxY, Eq((int)maxY + 1));
    EXPECT_THAT(box.maxZ, Eq((int)maxZ + 1));

}

TEST(BoundingBoxExactTest, findMinMaxFromTriangle)
{
    BoundingBox<real> box = BoundingBox<real>::makeInvalidMinMaxBox();

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
    BoundingBox<real> box = BoundingBox<real>();

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
    BoundingBox<real> box = BoundingBox<real>();

    box.minX = 0.0f;
    box.minY = 0.0f;
    box.minZ = 0.0f;

    box.maxX = 10.0f;
    box.maxY = 10.0f;
    box.maxZ = 10.0f;

    Triangle t = Triangle(Vertex(1, 1, 1), Vertex(2, 2, 2), Vertex(3, 3, 11), Vertex(0.0f, 0.0f, 0.0f));

    EXPECT_FALSE(box.isInside(t));
}

TEST(BoundingBoxTest, createNodeBoxWithFloastingPointValues)
{
    Triangle t = Triangle(Vertex(1.34, -2.01, 1.8), Vertex(2, 2, 1.9), Vertex(3.99, 2.1, 1.51), Vertex(0.0, 0.0, 0.0));
    real delta = 0.5;

    BoundingBox<real> box = BoundingBox<real>::makeRealNodeBox(t, delta);

    EXPECT_THAT(box.minX, RealEq(1.0));
    EXPECT_THAT(box.minY, RealEq(-2.5));
    EXPECT_THAT(box.minZ, RealEq(1.5));

    EXPECT_THAT(box.maxX, RealEq(4));
    EXPECT_THAT(box.maxY, RealEq(2.5));
    EXPECT_THAT(box.maxZ, RealEq(2.0));
}

