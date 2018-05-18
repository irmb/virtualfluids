#ifndef TEST_CONDITION_H
#define TEST_CONDITION_H

#include "VirtualFluids_GPU/LBM/LB.h"

#include <string>
#include <memory>

class Parameter;
class GridProvider;
class DataWriter;
class Calculator;
class TestResults;

class TestCondition
{
public:
	virtual std::shared_ptr<Parameter> getParameter() = 0;
	virtual std::shared_ptr<GridProvider> getGrid() = 0;
	virtual std::shared_ptr<DataWriter> getDataWriter() = 0;
	virtual std::shared_ptr<Calculator> getCalculator() = 0;
private:

};
#endif