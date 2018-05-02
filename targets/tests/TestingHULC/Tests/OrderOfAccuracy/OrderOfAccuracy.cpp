#include "OrderOfAccuracy.h"

#include "Tests\TestCout\TestCout.h"
#include "Tests\DataQueue\DataQueue.h"

TEST_P(OrderOfAccuracy, Test) {
	DataQueue input = GetParam();
	if (input.expected) {
		TEST_COUT(input.testName, input.la, input.lb, input.valueName, input.valueName, "OrderOfAccuracy", input.a, input.b, input.orderOfAccuracy);
	}
	ASSERT_THAT(OrderOfAccuracy::test(input.orderOfAccuracy, input.minOrderOfAccuracy), Eq(input.expected));
}

bool OrderOfAccuracy::test(double a, double b)
{
	return a > b;
}
