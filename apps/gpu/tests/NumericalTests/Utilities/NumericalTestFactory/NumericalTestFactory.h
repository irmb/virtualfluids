#ifndef NUMERICAL_TEST_FACTORY_H
#define NUMERICAL_TEST_FACTORY_H

#include <memory>
#include <vector>
#include <string>

class TestSimulation;
class TestQueue;
class LogFileQueue;

class NumericalTestFactory
{
public:
	virtual std::vector<std::shared_ptr<TestSimulation> > getTestSimulations() = 0;
	virtual std::shared_ptr<TestQueue> getTestQueue() = 0;
	virtual std::shared_ptr<LogFileQueue> getLogFileQueue() = 0;
private:

};
#endif 