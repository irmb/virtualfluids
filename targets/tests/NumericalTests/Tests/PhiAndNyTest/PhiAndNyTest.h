#ifndef PHI_AND_NY_TEST_H
#define PHI_AND_NY_TEST_H

#include "Utilities\Test\TestImp.h"

#include <memory>
#include <vector>
#include <iostream>

class FFTCalculator;
class PhiAndNyTestPostProcessingStrategy;
struct PhiAndNyTestParameterStruct;

class PhiAndNyTest : public TestImp 
{
public:
	static std::shared_ptr<PhiAndNyTest> getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<PhiAndNyTestParameterStruct> testPara, std::string dataToCalculate);
	
	void update();
	void addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<PhiAndNyTestPostProcessingStrategy> postProStrategy);
	void evaluate();
	std::vector<bool> getPassedTests();
	void makeConsoleOutput();

	std::string getDataToCalculate();
	std::vector<int> getLx();
	std::vector<double> getNy();
	std::vector<double> getNyDiff();
	std::vector<double> getPhiDiff();
	double getOrderOfAccuracyNyDiff();
	double getOrderOfAccuracyPhiDiff();



private:
	PhiAndNyTest(std::shared_ptr<ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<PhiAndNyTestParameterStruct> testPara, std::string dataToCalculate);
	double calcOrderOfAccuracy(std::vector<double> data);
	bool checkTestPassed(double orderOfAccuracy);
	std::vector<double> calcNyDiff(std::vector<double> ny);

	unsigned int startStepCalculation, endStepCalculation;
	std::vector<double> lx;
	std::vector<double> phiDiff;
	std::vector<double> ny, nyDiff;
	double orderOfAccuracyPhiDiff;
	double orderOfAccuracyNyDiff;
	double minOrderOfAccuracy;
	double viscosity;
	bool phiDiffTestPassed;
	bool nyDiffTestPassed;
	std::string dataToCalculate;

	std::vector<std::shared_ptr<PhiAndNyTestPostProcessingStrategy> > phiAndNyPostProStrategies;

};
#endif
