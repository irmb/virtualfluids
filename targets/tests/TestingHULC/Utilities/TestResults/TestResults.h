#ifndef TEST_RESULTS_H
#define TEST_RESULTS_H

class TestResults
{
public:
	virtual void evaluate() = 0;
	virtual void makeFinalOutput() = 0;
};
#endif
