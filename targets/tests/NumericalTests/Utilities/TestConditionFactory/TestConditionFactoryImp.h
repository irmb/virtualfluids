#ifndef TEST_CONDITION_FACTORY_IMP_H
#define TEST_CONDITION_FACTORY_IMP_H

#include "TestConditionFactory.h"


class TestConditionFactoryImp: public TestConditionFactory
{
public:
	static std::shared_ptr< TestConditionFactory> getNewInstance();
	std::vector< std::shared_ptr< TestCondition> > makeTestConditions(std::vector< std::shared_ptr< SimulationParameter> > simPara);

protected:
	TestConditionFactoryImp();

private:
	std::vector< std::shared_ptr< SimulationParameter> > simPara;
};
#endif