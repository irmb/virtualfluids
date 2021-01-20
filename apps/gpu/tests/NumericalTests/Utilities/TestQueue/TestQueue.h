#ifndef TEST_QUEUE_H
#define TEST_QUEUE_H

class TestQueue
{
public:
	virtual void makeFinalOutput() = 0;
    virtual int getNumberOfFailedTests() const noexcept = 0;
};
#endif