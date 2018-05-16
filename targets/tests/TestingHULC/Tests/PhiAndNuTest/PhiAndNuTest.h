#ifndef PHI_AND_NU_TEST_RESULTS_H
#define PHI_AND_NU_TEST_RESULTS_H

#include "Utilities/TestResults/TestResults.h"
#include "Utilities/LogFileInformation/LogFileInformation.h"

#include <memory>
#include <vector>
#include <iostream>

class TestCout;

class PhiAndNuTest : public TestResults, public LogFileInformation
{
public:
	static std::shared_ptr<PhiAndNuTest> getNewInstance(std::string aTestName, double minOrderOfAccuracy, std::shared_ptr<TestCout> testOut);
	void evaluate();
	void makeFinalOutput();
	void add(double phiDiff, double nuDiff, double lx);
	std::string getOutput();

private:
	PhiAndNuTest(std::string aTestName, double minOrderOfAccuracy, std::shared_ptr<TestCout> testOut);
	void makeLastTestOutput();
	std::vector<double> calcOrderOfAccuracy(std::vector<double> data);
	std::vector<bool> checkTestPassed(std::vector<double> orderOfAccuracy);

	std::vector<double> phiDiff;
	std::vector<double> nuDiff;
	std::vector<double> lx;
	std::vector<double> orderOfAccuracyPhiDiff;
	std::vector<double> orderOfAccuracyNuDiff;
	std::vector<bool> phiDiffTestPassed;
	std::vector<bool> nuDiffTestPassed;
		
	double minOrderOfAccuracy;
	std::string testName;

	std::shared_ptr<TestCout> testOut;
};
#endif
