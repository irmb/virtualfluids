#include "gmock/gmock.h"


#include "Conglomerate.h"
#include "../ObjectMocks.h"

using namespace testing;

TEST(ConglomerateTest, addTwoObjectsToConglomerate_ShouldReturnMinimumMaximum)
{
    auto objectStub1 = new ObjectStub();
    auto objectStub2 = new ObjectStub();
    objectStub1->setX1Minimum(0.0);
    objectStub2->setX1Minimum(-2.0);
    objectStub1->setX2Minimum(2.0);
    objectStub2->setX2Minimum(4.0);
    objectStub1->setX3Minimum(100.0);
    objectStub2->setX3Minimum(-100.0);

    objectStub1->setX1Maximum(0.0);
    objectStub2->setX1Maximum(-2.0);
    objectStub1->setX2Maximum(-5.0);
    objectStub2->setX2Maximum(-2.0);
    objectStub1->setX3Maximum(0.0);
    objectStub2->setX3Maximum(33.0);


    SPtr<Conglomerate> conglomerate = Conglomerate::makeShared();
    conglomerate->add(objectStub1);
    conglomerate->add(objectStub2);


    const real expectedMinX1 = -2.0;
    const real expectedMinX2 = 2.0;
    const real expectedMinX3 = -100.0;
    const real actualMinX1 = conglomerate->getX1Minimum();
    const real actualMinX2 = conglomerate->getX2Minimum();
    const real actualMinX3 = conglomerate->getX3Minimum();
    EXPECT_THAT(actualMinX1, RealEq(expectedMinX1));
    EXPECT_THAT(actualMinX2, RealEq(expectedMinX2));
    EXPECT_THAT(actualMinX3, RealEq(expectedMinX3));

    const real expectedMaxX1 = 0.0;
    const real expectedMaxX2 = -2.0;
    const real expectedMaxX3 = 33.0;
    const real actualMaxX1 = conglomerate->getX1Maximum();
    const real actualMaxX2 = conglomerate->getX2Maximum();
    const real actualMaxX3 = conglomerate->getX3Maximum();
    EXPECT_THAT(actualMaxX1, RealEq(expectedMaxX1));
    EXPECT_THAT(actualMaxX2, RealEq(expectedMaxX2));
    EXPECT_THAT(actualMaxX3, RealEq(expectedMaxX3));
}

TEST(ConglomerateTest, getCenterFromConglomerate)
{
    auto objectStub1 = new ObjectStub();
    auto objectStub2 = new ObjectStub();
    objectStub1->setX1Minimum(0.0);
    objectStub2->setX1Minimum(-2.0);
    objectStub1->setX2Minimum(2.0);
    objectStub2->setX2Minimum(4.0);
    objectStub1->setX3Minimum(100.0);
    objectStub2->setX3Minimum(-100.0);

    objectStub1->setX1Maximum(0.0);
    objectStub2->setX1Maximum(-2.0);
    objectStub1->setX2Maximum(12.0);
    objectStub2->setX2Maximum(10.0);
    objectStub1->setX3Maximum(0.0);
    objectStub2->setX3Maximum(100.0);


    SPtr<Conglomerate> conglomerate = Conglomerate::makeShared();
    conglomerate->add(objectStub1);
    conglomerate->add(objectStub2);


    const real expectedCenterX1 = -1.0;
    const real expectedCenterX2 = 7.0;
    const real expectedCenterX3 = 0.0;
    const real actualCenterX1 = conglomerate->getX1Centroid();
    const real actualCenterX2 = conglomerate->getX2Centroid();
    const real actualCenterX3 = conglomerate->getX3Centroid();
    EXPECT_THAT(actualCenterX1, RealEq(expectedCenterX1));
    EXPECT_THAT(actualCenterX2, RealEq(expectedCenterX2));
    EXPECT_THAT(actualCenterX3, RealEq(expectedCenterX3));
}


TEST(ConglomerateTest, getPointInObject_ShouldReturnTrueIfOnePointIsInObject)
{
    auto objectStub1 = new ObjectStub();
    auto objectStub2 = new ObjectStub();
    objectStub1->setIsInObject(false);
    objectStub2->setIsInObject(true);


    SPtr<Conglomerate> conglomerate = Conglomerate::makeShared();
    conglomerate->add(objectStub1);
    conglomerate->add(objectStub2);


    const bool expectedIsInObject = true;
    const bool actualIsInObject = conglomerate->isPointInObject(0,0,0,0,0);
    EXPECT_THAT(actualIsInObject, Eq(expectedIsInObject));
}

TEST(ConglomerateTest, getPointInObject_ShouldReturnTrueIfOnePointIsInBothObject)
{
    auto objectStub1 = new ObjectStub();
    auto objectStub2 = new ObjectStub();
    objectStub1->setIsInObject(true);
    objectStub2->setIsInObject(true);


    SPtr<Conglomerate> conglomerate = Conglomerate::makeShared();
    conglomerate->add(objectStub1);
    conglomerate->add(objectStub2);


    const bool expectedIsInObject = true;
    const bool actualIsInObject = conglomerate->isPointInObject(0, 0, 0, 0, 0);
    EXPECT_THAT(actualIsInObject, Eq(expectedIsInObject));
}

TEST(ConglomerateTest, addAndSubtractTwoObjects_PointIsInBothObjects_ShouldNotBeInConglomerate)
{
    auto objectStub1 = new ObjectStub();
    auto objectStub2 = new ObjectStub();
    objectStub1->setIsInObject(true);
    objectStub2->setIsInObject(true);


    SPtr<Conglomerate> conglomerate = Conglomerate::makeShared();
    conglomerate->add(objectStub1);
    conglomerate->subtract(objectStub2);


    const bool expectedIsInObject = false;
    const bool actualIsInObject = conglomerate->isPointInObject(0, 0, 0, 0, 0);
    EXPECT_THAT(actualIsInObject, Eq(expectedIsInObject));
}