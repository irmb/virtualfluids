#include "TestQueueImp.h"

#include "Utilities\ColorConsoleOutput\ColorConsoleOutput.h"
#include "Utilities\Test\Test.h"

void TestQueueImp::makeFinalOutput()
{
	calcNumberOfPassedTest();
	colorOutput->makeFinalTestOutputHead(numberOfPassedTest, numberOfTests);
	for (int i = 0; i < tests.size(); i++)
		tests.at(i)->makeConsoleOutput();
	colorOutput->makeFinalTestOutputFoot(numberOfPassedTest, numberOfTests);
}

std::shared_ptr<TestQueueImp> TestQueueImp::getNewInstance(std::shared_ptr< ColorConsoleOutput> colorOutput)
{
	return std::shared_ptr<TestQueueImp>(new TestQueueImp(colorOutput));
}

void TestQueueImp::addTest(std::shared_ptr<Test> test)
{
	tests.push_back(test);
}

TestQueueImp::TestQueueImp(std::shared_ptr< ColorConsoleOutput> colorOutput) : colorOutput(colorOutput)
{
	tests.resize(0);
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
