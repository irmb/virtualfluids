#ifndef NUMERICAL_TEST_SUITE_H
#define NUMERICAL_TEST_SUITE_H

class NumericalTestSuite
{
public:
	virtual void makeSimulationHeadOutput() = 0;
	virtual void startPostProcessing() = 0;
};
#endif