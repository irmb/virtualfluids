#include "gmock/gmock.h"


#include "Conglomerate.h"
#include "../ObjectMocks.h"

using namespace testing;

TEST(ConglomerateTest, addTwoObjectsToConglomerate_ShouldReturnMinimum)
{
    ObjectStub* objectStub1 = new ObjectStub();
    ObjectStub* objectStub2 = new ObjectStub();
    objectStub1->setX1Minimum(0.0);
    objectStub2->setX1Minimum(-2.0);


    Conglomerate *conglomerate = new Conglomerate();
    conglomerate->add(objectStub1);
    conglomerate->add(objectStub2);


    const real expected = -2.0;
    const real actual = conglomerate->getX1Minimum();
    EXPECT_THAT(actual, RealEq(expected));
}
