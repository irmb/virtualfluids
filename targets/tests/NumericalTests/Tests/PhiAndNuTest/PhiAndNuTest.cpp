#include "PhiAndNuTest.h"

#include "Utilities/ColorConsoleOutput/ColorConsoleOutput.h"
#include "Utilities\PostProcessingResults\PostProcessingResults.h"
#include "Utilities\TestSimulation\TestSimulation.h"
#include "Utilities\SimulationInfo\SimulationInfo.h"

#include <iomanip>

std::shared_ptr<PhiAndNuTest> PhiAndNuTest::getNewInstance(std::shared_ptr< ColorConsoleOutput> colorOutput, std::string dataToCalculate, double minOrderOfAccuracy, double viscosity, unsigned int startStepCalculation, unsigned int endStepCalculation)
{
	return std::shared_ptr<PhiAndNuTest>(new PhiAndNuTest(colorOutput, dataToCalculate, minOrderOfAccuracy, viscosity, startStepCalculation, endStepCalculation));
}

void PhiAndNuTest::evaluate()
{
	for (int i = 0; i < postProResults.size(); i++) {
		lx.push_back(postProResults.at(i)->getNumberOfXNodes());
		if (dataToCalculate == "Vx") {
			nu.push_back(postProResults.at(i)->getNuVx());
			phiDiff.push_back(postProResults.at(i)->getPhiDiffVx());
		}
		if (dataToCalculate == "Vy") {
			nu.push_back(postProResults.at(i)->getNuVy());
			phiDiff.push_back(postProResults.at(i)->getPhiDiffVy());
		}
		if (dataToCalculate == "Vz") {
			nu.push_back(postProResults.at(i)->getNuVz());
			phiDiff.push_back(postProResults.at(i)->getPhiDiffVz());
		}
	}
	nuDiff = calcNuDiff(nu);
	orderOfAccuracyPhiDiff = calcOrderOfAccuracy(phiDiff);
	orderOfAccuracyNuDiff = calcOrderOfAccuracy(nuDiff);
	phiDiffTestPassed = checkTestPassed(orderOfAccuracyPhiDiff);
	nuDiffTestPassed = checkTestPassed(orderOfAccuracyNuDiff);

	makeConsoleOutput();
}

void PhiAndNuTest::update()
{
	TestImp::update();
}

void PhiAndNuTest::addSimulation(std::shared_ptr<TestSimulation> sim, std::shared_ptr< SimulationInfo> simInfo, std::shared_ptr< PostProcessingResults> postProResults)
{
	TestImp::addSimulation(sim, simInfo, postProResults);
}

std::vector<bool> PhiAndNuTest::getPassedTests()
{
	std::vector< bool> passed;
	passed.push_back(phiDiffTestPassed);
	passed.push_back(nuDiffTestPassed);
	return passed;
}

void PhiAndNuTest::makeConsoleOutput()
{
	colorOutput->makeNuTestOutput(nuDiffTestPassed, simInfos.at(0), simInfos.at(1), startStepCalculation, endStepCalculation, dataToCalculate, nu.at(0), nu.at(1), nuDiff.at(0), nuDiff.at(1), orderOfAccuracyNuDiff);
	colorOutput->makePhiTestOutput(nuDiffTestPassed, simInfos.at(0), simInfos.at(1), startStepCalculation, endStepCalculation, dataToCalculate, phiDiff.at(0), phiDiff.at(1), orderOfAccuracyPhiDiff);
}

std::string PhiAndNuTest::getLogFileOutput()
{
	std::ostringstream oss;

	oss << "Nu_" << lx.at(0) << "=" << nu.at(0) << std::endl;
	oss << "NuDiff_" << lx.at(0) << "=" << nuDiff.at(0) << std::endl;
	oss << "Nu_" << lx.at(1) << "=" << nu.at(1) << std::endl;
	oss << "NuDiff_" << lx.at(1) << "=" << nuDiff.at(1) << std::endl;
	oss << "OrderOfAccuracy_NuDiff_" << lx.at(0) << lx.at(1) << "=" << orderOfAccuracyNuDiff << std::endl << std::endl;

	oss << "PhiDiff_" << lx.at(0) << "=" << phiDiff.at(0) << std::endl;
	oss << "PhiDiff_" << lx.at(1) << "=" << phiDiff.at(1) << std::endl;
	oss << "OrderOfAccuracy_NuDiff_" << lx.at(0) << lx.at(1) << "=" << orderOfAccuracyPhiDiff << std::endl << std::endl;

	return oss.str();
}

PhiAndNuTest::PhiAndNuTest(std::shared_ptr< ColorConsoleOutput> colorOutput, std::string dataToCalculate, double minOrderOfAccuracy, double viscosity, unsigned int startStepCalculation, unsigned int endStepCalculation) : TestImp(colorOutput), minOrderOfAccuracy(minOrderOfAccuracy), viscosity(viscosity), dataToCalculate(dataToCalculate), startStepCalculation(startStepCalculation), endStepCalculation(endStepCalculation)
{
	lx.resize(0);
	phiDiff.resize(0);
	nuDiff.resize(0);
}

double PhiAndNuTest::calcOrderOfAccuracy(std::vector<double> data)
{
	double ooa = log(data.at(0) / data.at(1)) / log(lx.at(1) / lx.at(0));
	
	return ooa;
}

bool PhiAndNuTest::checkTestPassed(double orderOfAccuracy)
{
	return orderOfAccuracy > minOrderOfAccuracy;
}

std::vector<double> PhiAndNuTest::calcNuDiff(std::vector<double> nu)
{
	std::vector< double> results;
	for(int i = 0; i < nu.size(); i++)
		results.push_back((nu.at(i) - viscosity) / viscosity);
	return results;
}
