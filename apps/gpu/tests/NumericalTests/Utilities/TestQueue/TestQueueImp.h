#ifndef TEST_QUEUE_IMP_H
#define TEST_QUEUE_IMP_H

#include "TestQueue.h"

#include <memory>
#include <vector>

class Test;
class ColorConsoleOutput;

class TestQueueImp : public TestQueue
{
public:
	TestSuiteResult run() override;
	void makeFinalOutput() override;

	int getNumberOfFailedTests() const noexcept override;

	static std::shared_ptr<TestQueueImp> getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput);
	void addTest(std::shared_ptr<Test> test);

private:
	TestQueueImp(std::shared_ptr<ColorConsoleOutput> colorOutput);

	void calcTestNumbers();

	std::shared_ptr<ColorConsoleOutput> colorOutput;
	std::vector<std::shared_ptr<Test> > tests;
	
	int numberOfPassedTest;
	int numberOfFailedTest;
	int numberOfErrorTest;
	int numberOfExecutedTest;
	int numberOfNotExecutedTest;

	int numberOfTests;
};
#endif