#ifndef NY_TEST_H
#define NY_TEST_H

#include "Utilities/Test/TestImp.h"

#include <memory>
#include <vector>
#include <iostream>

class FFTCalculator;
class NyTestPostProcessingStrategy;
struct NyTestParameterStruct;

class NyTest : public TestImp 
{
public:
	static std::shared_ptr<NyTest> getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<NyTestParameterStruct> testPara, std::string dataToCalculate);
	
	void update();
	void addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<NyTestPostProcessingStrategy> postProStrategy);
	
	void evaluate();

	std::string getDataToCalculate();
	std::vector<int> getLx();
	std::vector<double> getNy();
	std::vector<double> getNyDiff();
	double getOrderOfAccuracyNyDiff();

private:
	NyTest(std::shared_ptr<ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<NyTestParameterStruct> testPara, std::string dataToCalculate);
	double calcOrderOfAccuracy(std::vector<double> data);
	TestStatus checkTestPassed(double orderOfAccuracy);
	bool checkNy(std::vector<double> ny);
	std::vector<double> calcNyDiff(std::vector<double> ny);
	std::vector<std::string> buildTestOutput();
	std::vector<std::string> buildBasicTestOutput();
	std::vector<std::string> buildErrorTestOutput();
	std::vector<std::string> buildSimulationFailedTestOutput();


	unsigned int startStepCalculation, endStepCalculation;
	std::vector<double> lx;
	std::vector<double> ny, nyDiff;
	double orderOfAccuracy;
	double minOrderOfAccuracy;
	double viscosity;
	std::string dataToCalculate;

	std::vector<std::shared_ptr<NyTestPostProcessingStrategy> > postProStrategies;

};
#endif
