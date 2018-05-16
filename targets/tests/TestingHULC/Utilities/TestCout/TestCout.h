#ifndef TEST_COUT_H
#define TEST_COUT_H

#include <iostream>

class TestCout
{
public:
	virtual void makeTestOutput(bool testPassed, std::string testName, int l1, int l2, std::string nameWerte1, std::string nameWerte2, std::string nameWerte3, double testWert1, double testWert2, double testWert3) = 0;
	virtual void makeSimulationHeadOutput(std::string simName, int l) = 0;

private:
};
#endif 