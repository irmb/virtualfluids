#ifndef PHIANDNUTESTRESULTS_H
#define PHIANDNUTESTRESULTS_H

#include "TestResults.h"
#include <memory>
#include <vector>
#include <iostream>

class PhiAndNuTestResults : public TestResults
{
public:
	static std::shared_ptr<PhiAndNuTestResults> getNewInstance(std::string aTestName);
	void evaluate();
	void add(double phiDiff, double nuDiff, double lx);

private:
	PhiAndNuTestResults(std::string aTestName);
	void makeLastTestOutput();

	std::vector<double> phiDiff;
	std::vector<double> nuDiff;
	std::vector<double> lx;
	std::vector<double> orderOfAccuracyPhiDiff;
	std::vector<double> orderOfAccuracyNuDiff;

	std::string testName;
};
#endif // !PHIANDNUTESTRESULTS
