#include "TransformatorImp.h"

#include "gmock/gmock.h"
#include <memory>
using namespace testing;

#include <GridGenerator/geometries/Vertex/Vertex.cuh>
#include <GridGenerator/geometries/Triangle/Triangle.cuh>
#include <GridGenerator/geometries/Geometry/Geometry.cuh>

#include <GridGenerator/geometries/Arrow/ArrowMocks.h>


class TransformatorTest : public Test 
{
public:
	std::shared_ptr<Transformator> sut;
	Vertex v;

	doubflo delta;
	Vertex translater;

	void SetUp() override 
	{
		delta = 0.01f;
		v = Vertex(2.0f, -3.2f, 4.6f);
		translater = Vertex(-2.6f, 3434.0f, 0.1f);
	}
};

void expectEqual(const Vertex& actual, const Vertex& expected)
{
	EXPECT_THAT(actual.x, FloatEq(expected.x));
	EXPECT_THAT(actual.y, FloatEq(expected.y));
	EXPECT_THAT(actual.z, FloatEq(expected.z));
}

TEST(TransformatorCopyConstructorTest, copyTransformatorShouldCreateSameTransformator)
{
	doubflo delta = 0.01f;
	doubflo dx = 0.1f;
	doubflo dy = 0.2f;
	doubflo dz = 0.3f;
	TransformatorImp trafoToCopy(delta, dx, dy, dz);

	std::shared_ptr<TransformatorImp> sut = std::shared_ptr<TransformatorImp>(new TransformatorImp(trafoToCopy));

	EXPECT_TRUE(*(sut.get())==trafoToCopy);
}

TEST_F(TransformatorTest, transformVectorToViewWithStandartConstructor_ExpectnoVectorChange) 
{
	sut = std::shared_ptr<Transformator>(new TransformatorImp());
	Vertex expected = Vertex(v);

    sut->transformWorldToGrid(v);

	expectEqual(v, expected);
}

TEST_F(TransformatorTest, transformVectorToViewWithSmallDelta_ExpectVectorScales)
{
	sut = std::shared_ptr<Transformator>(new TransformatorImp(delta, Vertex()));
	Vertex expected = Vertex(v * (1.0f / delta));

	sut->transformWorldToGrid(v);

	expectEqual(v, expected);
}

TEST_F(TransformatorTest, transformVectorWithNullDelta_ExpectExcpetion) 
{
    doubflo invalidDeltaValue = 0.0f;
    ASSERT_THROW(TransformatorImp trafo(invalidDeltaValue, Vertex(0, 0, 0)), invalidDelta);
}

TEST_F(TransformatorTest, transformVectorWithNegativeDelta_ExpectExcpetion)
{
    doubflo invalidDeltaValue = -1.0f;
    ASSERT_THROW(TransformatorImp trafo(invalidDeltaValue, Vertex(0, 0, 0)), invalidDelta);
}

TEST_F(TransformatorTest, transformVectorToViewWithTranslationsAndSmallDelta)
{
	sut = std::shared_ptr<Transformator>(new TransformatorImp(delta, translater));
	Vertex expected = Vertex((v + translater) * (1.0f / delta));

	sut->transformWorldToGrid(v);

	expectEqual(v, expected);
}

TEST_F(TransformatorTest, transformVectorToWorldCoodinatesWithTranslationsAndSmallDelta)
{
	sut = std::shared_ptr<Transformator>(new TransformatorImp(delta, translater));
	Vertex expected = Vertex(v * delta - translater);

	sut->transformGridToWorld(v);

	expectEqual(v, expected);
}

TEST_F(TransformatorTest, transformTriangleToView)
{
	sut = std::shared_ptr<Transformator>(new TransformatorImp(delta, translater));
	Triangle t(v, v, v);
	Vertex expected = Vertex((v + translater) * (1.0f / delta));

	sut->transformWorldToGrid(t);

	expectEqual(t.v1, expected);
	expectEqual(t.v2, expected);
	expectEqual(t.v3, expected);
}

TEST_F(TransformatorTest, transformGeometryToView)
{
	sut = std::shared_ptr<Transformator>(new TransformatorImp(delta, translater));
	Geometry g;
	g.triangleVec.push_back(Triangle(v,v,v));
	g.size = 1;
	Vertex expected = Vertex((v + translater) * (1.0f / delta));

	sut->transformWorldToGrid(g);

	expectEqual(g.triangleVec[0].v1, expected);
	expectEqual(g.triangleVec[0].v2, expected);
	expectEqual(g.triangleVec[0].v3, expected);
}

TEST(TransformatorTestBoundingBox, transformDoubfloBoundingBoxToView)
{
	doubflo delta = 0.01f;
	Vertex translater = Vertex(-2.6f, 3434.0f, 0.1f);
	std::shared_ptr<Transformator> sut = std::shared_ptr<Transformator>(new TransformatorImp(delta, translater));

	BoundingBox<doubflo> box(0, 0, 0, 0, 0, 0);

	sut->transformWorldToGrid(box);

	EXPECT_THAT(box.minX, Eq(translater.x * (1.0f / delta)));
	EXPECT_THAT(box.minY, Eq(translater.y * (1.0f / delta)));
	EXPECT_THAT(box.minZ, Eq(translater.z * (1.0f / delta)));

	EXPECT_THAT(box.maxX, Eq(translater.x * (1.0f / delta)));
	EXPECT_THAT(box.maxY, Eq(translater.y * (1.0f / delta)));
	EXPECT_THAT(box.maxZ, Eq(translater.z * (1.0f / delta)));
}

TEST(TransformatorTestBoundingBox, transformDoubfloBoundingBoxToWorld)
{
	doubflo delta = 0.01f;
	Vertex translater = Vertex(-2.6f, 3434.0f, 0.1f);
	std::shared_ptr<Transformator> sut = std::shared_ptr<Transformator>(new TransformatorImp(delta, translater));

	BoundingBox<doubflo> box(0, 0, 0, 0, 0, 0);

	sut->transformGridToWorld(box);

	EXPECT_THAT(box.minX, Eq(-1.0f * translater.x));
	EXPECT_THAT(box.minY, Eq(-1.0f * translater.y));
	EXPECT_THAT(box.minZ, Eq(-1.0f * translater.z));
												  
	EXPECT_THAT(box.maxX, Eq(-1.0f * translater.x));
	EXPECT_THAT(box.maxY, Eq(-1.0f * translater.y));
	EXPECT_THAT(box.maxZ, Eq(-1.0f * translater.z));
}

TEST_F(TransformatorTest, transformArrowToWorld)
{
	delta = 10.0f;
	translater = Vertex(0, 0, 0);
	std::shared_ptr<ArrowTransformator> sut = std::shared_ptr<ArrowTransformator>(new TransformatorImp(delta, translater));

	doubflo x1 = 1.23f;
	doubflo y1 = 2.23f;
	doubflo z1 = 0.023f;

	doubflo x2 = -1.23f;
	doubflo y2 = 1.23f;
	doubflo z2 = -0.20233f;
	auto v1 = std::shared_ptr<Vertex>(new Vertex(x1, y1, z1));
	auto v2 = std::shared_ptr<Vertex>(new Vertex(x2, y2, z2));
	auto arrow = ArrowStub::make(v1, v2);

	sut->transformGridToWorld(arrow);

	EXPECT_THAT(arrow->getStart()->x, FloatEq(x1 * delta - translater.x));
	EXPECT_THAT(arrow->getStart()->y, FloatEq(y1 * delta - translater.y));
	EXPECT_THAT(arrow->getStart()->z, FloatEq(z1 * delta - translater.z));

	EXPECT_THAT(arrow->getEnd()->x, FloatEq(x2 * delta - translater.x));
	EXPECT_THAT(arrow->getEnd()->y, FloatEq(y2 * delta - translater.y));
	EXPECT_THAT(arrow->getEnd()->z, FloatEq(z2 * delta - translater.z));
}
