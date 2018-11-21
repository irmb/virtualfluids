#include "PhiAndNuTest.h"

#include "Utilities/ColorConsoleOutput/ColorConsoleOutputImp.h"
#include "Utilities/Results/SimulationResults/SimulationResults.h"
#include "Utilities\Calculator\FFTCalculator\FFTCalculator.h"
#include "Utilities\TestSimulation\TestSimulation.h"
#include "Utilities\SimulationInfo\SimulationInfo.h"

std::shared_ptr<PhiAndNuTest> PhiAndNuTest::getNewInstance(std::string dataToCalculate, double minOrderOfAccuracy, double viscosity)
{
	return std::shared_ptr<PhiAndNuTest>(new PhiAndNuTest(dataToCalculate, minOrderOfAccuracy, viscosity));
}

void PhiAndNuTest::evaluate()
{
	for (int i = 0; i < simResults.size(); i++) {
		lx.push_back(simResults.at(i)->getNumberOfXNodes());
		calculator->setSimulationResults(simResults.at(i));
		if (dataToCalculate == "Vx")
			calculator->setVectorToCalc(simResults.at(i)->getVx());
		if (dataToCalculate == "Vz")
			calculator->setVectorToCalc(simResults.at(i)->getVz());
		calculator->calc();
		phiDiff.push_back(calculator->getPhiDiff());
		nuDiff.push_back(calculator->getNuDiff());
	}
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

void PhiAndNuTest::addSimulation(std::shared_ptr<TestSimulation> sim, std::shared_ptr< SimulationInfo> simInfo)
{
	TestImp::addSimulation(sim, simInfo);
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
	testOut->makeTestOutput(nuDiffTestPassed, simInfos.at(0), simInfos.at(1), "NuDiff", "NuDiff", "OrderOfAccuracy", nuDiff.at(0), nuDiff.at(1), orderOfAccuracyNuDiff);
	testOut->makeTestOutput(nuDiffTestPassed, simInfos.at(0), simInfos.at(1), "PhiDiff", "PhiDiff", "OrderOfAccuracy", phiDiff.at(0), phiDiff.at(1), orderOfAccuracyPhiDiff);
}

std::string PhiAndNuTest::getLogFileOutput()
{
	std::ostringstream oss;
	oss << std::setfill(' ') << std::left << std::setw(4) << lx.at(0) << std::setw(45) << nuDiff.at(0) << phiDiff.at(0) << std::endl;
	oss << std::setfill(' ') << std::left << std::setw(19) << " " << std::setw(45) << orderOfAccuracyNuDiff << orderOfAccuracyPhiDiff << std::endl;
	oss << std::setfill(' ') << std::left << std::setw(4) << lx.at(1) << std::setw(45) << nuDiff.at(1) << phiDiff.at(1) << std::endl;

	return oss.str();
}

PhiAndNuTest::PhiAndNuTest(std::string dataToCalculate, double minOrderOfAccuracy, double viscosity) : TestImp(), minOrderOfAccuracy(minOrderOfAccuracy), viscosity(viscosity), dataToCalculate(dataToCalculate)
{
	lx.resize(0);
	phiDiff.resize(0);
	nuDiff.resize(0);
	calculator = FFTCalculator::getNewInstance(viscosity);
	testOut = ColorConsoleOutputImp::getNewInstance();
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