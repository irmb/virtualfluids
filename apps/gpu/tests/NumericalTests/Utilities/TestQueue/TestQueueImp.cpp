#include "TestQueueImp.h"
#include <algorithm>

#include "Utilities/ColorConsoleOutput/ColorConsoleOutput.h"
#include "Utilities/Test/Test.h"

#include "basics/DataTypes.h"

TestSuiteResult TestQueueImp::run()
{
    for (const auto& test : tests)
        test->run();

    makeFinalOutput();

    return TestSuiteResult(std::clamp(numberOfFailedTest, 0, 1));
}

void TestQueueImp::makeFinalOutput()
{
    calcTestNumbers();
    colorOutput->makeFinalTestOutputHead(numberOfTests, numberOfExecutedTest, numberOfPassedTest, numberOfFailedTest,
                                         numberOfErrorTest, numberOfNotExecutedTest);
    for (uint i = 0; i < tests.size(); i++)
        tests.at(i)->makeConsoleOutput();
    colorOutput->makeFinalTestOutputFoot(numberOfTests, numberOfExecutedTest, numberOfPassedTest, numberOfFailedTest,
                                         numberOfErrorTest, numberOfNotExecutedTest);
}

int TestQueueImp::getNumberOfFailedTests() const noexcept
{
    return numberOfFailedTest;
}

std::shared_ptr<TestQueueImp> TestQueueImp::getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput)
{
    return std::shared_ptr<TestQueueImp>(new TestQueueImp(colorOutput));
}

void TestQueueImp::addTest(std::shared_ptr<Test> test)
{
    tests.push_back(test);
}

TestQueueImp::TestQueueImp(std::shared_ptr<ColorConsoleOutput> colorOutput) : colorOutput(colorOutput)
{
    tests.resize(0);
}

void TestQueueImp::calcTestNumbers()
{
    numberOfTests = tests.size();
    numberOfExecutedTest = 0;
    numberOfPassedTest = 0;
    numberOfFailedTest = 0;
    numberOfErrorTest = 0;
    numberOfNotExecutedTest = 0;

    for (uint i = 0; i < tests.size(); i++) {
        switch (tests.at(i)->getTestStatus()) {
            case passed:
                numberOfPassedTest++;
                numberOfExecutedTest++;
                break;
            case failed:
                numberOfFailedTest++;
                numberOfExecutedTest++;
                break;
            case test_error:
                numberOfErrorTest++;
                break;
            case simulationCrashed:
                numberOfNotExecutedTest++;
                break;
            default:
                break;
        }
    }
}
