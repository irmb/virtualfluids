#ifndef PHI_AND_NU_TEST_H
#define PHI_AND_NU_TEST_H

#include "Utilities\Test\TestImp.h"

#include <memory>
#include <vector>
#include <iostream>

class FFTCalculator;

class PhiAndNuTest : public TestImp 
{
public:
	static std::shared_ptr<PhiAndNuTest> getNewInstance(std::shared_ptr< ColorConsoleOutput> colorOutput, std::string dataToCalculate, double minOrderOfAccuracy, double viscosity);
	
	void update();
	void addSimulation(std::shared_ptr< TestSimulation> sim, std::shared_ptr< SimulationInfo> simInfo);
	void evaluate();
	std::string getLogFileOutput();
	std::vector< bool> getPassedTests();
	void makeConsoleOutput();

private:
	PhiAndNuTest(std::shared_ptr< ColorConsoleOutput> colorOutput, std::string dataToCalculate, double minOrderOfAccuracy, double viscosity);
	double calcOrderOfAccuracy(std::vector<double> data);
	bool checkTestPassed(double orderOfAccuracy);
	

	std::shared_ptr< FFTCalculator> calculator;
	
	std::vector<double> lx;
	std::vector<double> phiDiff;
	std::vector<double> nuDiff;
	double orderOfAccuracyPhiDiff;
	double orderOfAccuracyNuDiff;
	double minOrderOfAccuracy;
	double viscosity;

	bool phiDiffTestPassed;
	bool nuDiffTestPassed;
	
	std::string dataToCalculate;
};
#endif
