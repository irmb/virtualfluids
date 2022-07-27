#ifndef TEST_QUEUE_H
#define TEST_QUEUE_H

enum TestSuiteResult { PASSED, FAILED };

class TestQueue
{
public:
    virtual TestSuiteResult run() = 0;
    virtual void makeFinalOutput() = 0;
    virtual int getNumberOfFailedTests() const noexcept = 0;
};
#endif