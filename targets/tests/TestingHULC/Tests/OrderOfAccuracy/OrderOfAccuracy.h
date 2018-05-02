#ifndef ORDEROFACCURACY_H
#define ORDEROFACCURACY_H

#include "gmock\gmock.h"

using namespace testing;

class DataQueue;

class OrderOfAccuracy : public ::testing::TestWithParam<DataQueue> {
public:
	bool test(double a, double b);

};
#endif 