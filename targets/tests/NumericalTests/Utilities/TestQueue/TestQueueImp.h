#ifndef TEST_QUEUE_IMP_H
#define TEST_QUEUE_IMP_H

#include "TestQueue.h"

#include <memory>
#include <vector>

class Test;
class TestCout;

class TestQueueImp : public TestQueue
{
public:
	void makeFinalOutput();

	static std::shared_ptr< TestQueueImp> getNewInstance();
	void addTest(std::shared_ptr< Test> test);

private:
	TestQueueImp();

	void calcNumberOfPassedTest();

	std::vector< std::shared_ptr< Test>> tests;
	std::shared_ptr< TestCout> colorOutput;
	int numberOfPassedTest;
	int numberOfTests;
};
#endif