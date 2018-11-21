#include "TestQueueImp.h"

#include "Utilities\TestCout\TestCoutImp.h"
#include "Utilities\Test\Test.h"

void TestQueueImp::makeFinalOutput()
{
	calcNumberOfPassedTest();
	colorOutput->makeFinalTestOutputHead(numberOfPassedTest, numberOfTests);
	for (int i = 0; i < tests.size(); i++)
		tests.at(i)->makeOutput();
	colorOutput->makeFinalTestOutputFoot(numberOfPassedTest, numberOfTests);
}

std::shared_ptr<TestQueueImp> TestQueueImp::getNewInstance()
{
	return std::shared_ptr<TestQueueImp>(new TestQueueImp());
}

void TestQueueImp::addTest(std::shared_ptr<Test> test)
{
	tests.push_back(test);
}

TestQueueImp::TestQueueImp()
{
	tests.resize(0);
	colorOutput = TestCoutImp::getNewInstance();
}

void TestQueueImp::calcNumberOfPassedTest()
{
	numberOfPassedTest = 0;
	numberOfTests = 0;

	for (int i = 0; i < tests.size(); i++) {
		for (int j = 0; j < tests.at(i)->getPassedTests().size(); j++) {
			if (tests.at(i)->getPassedTests().at(j)) {
				numberOfPassedTest++;
				numberOfTests++;
			}
			else
			{
				numberOfTests++;
			}
		}
			
	}
}
