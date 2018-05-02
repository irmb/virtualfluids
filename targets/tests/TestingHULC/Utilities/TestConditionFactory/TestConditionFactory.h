#ifndef TEST_CONDITION_FACTORY_H
#define TEST_CONDITION_FACTORY_H

#include <memory>
#include <vector>

class TestCondition;
class TestParameter;

class TestConditionFactory
{
public:
	virtual std::vector<std::shared_ptr<TestCondition>> makeTestConditions() = 0;
};
#endif
