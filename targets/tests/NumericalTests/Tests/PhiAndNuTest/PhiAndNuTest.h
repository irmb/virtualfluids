#ifndef PHI_AND_NU_TEST_H
#define PHI_AND_NU_TEST_H

#include "Utilities\Test\TestImp.h"

#include <memory>
#include <vector>
#include <iostream>

class FFTCalculator;
class PhiAndNuTestPostProcessingStrategy;
struct PhiAndNuTestParameterStruct;

class PhiAndNuTest : public TestImp 
{
public:
	static std::shared_ptr<PhiAndNuTest> getNewInstance(std::shared_ptr< ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<PhiAndNuTestParameterStruct> testPara, std::string dataToCalculate);
	
	void update();
	void addSimulation(std::shared_ptr< TestSimulation> sim, std::shared_ptr< SimulationInfo> simInfo, std::shared_ptr< PhiAndNuTestPostProcessingStrategy> postProStrategy);
	void evaluate();
	std::vector< bool> getPassedTests();
	void makeConsoleOutput();

	std::string getDataToCalculate();
	std::vector<int> getLx();
	std::vector<double> getNu();
	std::vector<double> getNuDiff();
	std::vector<double> getPhiDiff();
	double getOrderOfAccuracyNuDiff();
	double getOrderOfAccuracyPhiDiff();



private:
	PhiAndNuTest(std::shared_ptr< ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<PhiAndNuTestParameterStruct> testPara, std::string dataToCalculate);
	double calcOrderOfAccuracy(std::vector<double> data);
	bool checkTestPassed(double orderOfAccuracy);
	std::vector< double> calcNuDiff(std::vector< double> nu);

	unsigned int startStepCalculation, endStepCalculation;
	std::vector<double> lx;
	std::vector<double> phiDiff;
	std::vector<double> nu, nuDiff;
	double orderOfAccuracyPhiDiff;
	double orderOfAccuracyNuDiff;
	double minOrderOfAccuracy;
	double viscosity;
	bool phiDiffTestPassed;
	bool nuDiffTestPassed;
	std::string dataToCalculate;

	std::vector<std::shared_ptr< PhiAndNuTestPostProcessingStrategy>> phiAndNuPostProStrategies;

};
#endif
