#ifndef TEST_CONDITION_FACTORY_IMP_H
#define TEST_CONDITION_FACTORY_IMP_H

#include "TestConditionFactory.h"


class TestConditionFactoryImp: public TestConditionFactory
{
public:
	static std::shared_ptr<TestConditionFactory> getNewInstance(std::vector < std::shared_ptr < TestParameter > > testPara);
	std::vector<std::shared_ptr<TestCondition>> makeTestConditions();

protected:
	TestConditionFactoryImp() {};
	TestConditionFactoryImp(std::vector < std::shared_ptr < TestParameter > > testPara);

private:
	std::vector < std::shared_ptr < TestParameter > > testPara;
};
#endif