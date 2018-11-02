#ifndef TEST_CONDITION_FACTORY_H
#define TEST_CONDITION_FACTORY_H

#include <memory>
#include <vector>

class TestCondition;
class SimulationParameter;

class TestConditionFactory
{
public:
	virtual std::vector< std::shared_ptr< TestCondition> > makeTestConditions(std::vector< std::shared_ptr< SimulationParameter> > simPara) = 0;
};
#endif
