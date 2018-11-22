#include "PhiAndNuLogFileInformation.h"

#include "Tests\PhiAndNuTest\PhiAndNuTest.h"

#include <iomanip>
#include <sstream>


std::shared_ptr<PhiAndNuInformation> PhiAndNuInformation::getNewInstance(std::vector<std::shared_ptr<PhiAndNuTest>> tests, unsigned int startTimeStepCalculation, unsigned int endTimeStepCalculation)
{
	return std::shared_ptr<PhiAndNuInformation>(new PhiAndNuInformation(tests, startTimeStepCalculation, endTimeStepCalculation));
}

std::string PhiAndNuInformation::getOutput()
{
	std::ostringstream headName;
	headName << tests.at(0)->getSimulationName() <<" Phi And Nu Test";
	makeCenterHead(headName.str());

	oss << "StartTimeStepCalculation=" << startTimeStepCalculation << std::endl;
	oss << "EndTimeStepCalculation=" << endTimeStepCalculation << std::endl;
	oss << std::endl;

	oss << std::setfill(' ') << std::left << std::setw(4) << "L" << std::setw(15) << "NuDiff" << std::setw(30) << "Order of Accuracy" << std::setw(15) << "PhiDiff" << "Order of Accuracy" << std::endl;

	for(int i = 0; i < tests.size(); i++)
		oss << tests.at(i)->getLogFileOutput();

	return oss.str();
}

PhiAndNuInformation::PhiAndNuInformation(std::vector<std::shared_ptr<PhiAndNuTest>> tests, unsigned int startTimeStepCalculation, unsigned int endTimeStepCalculation) : startTimeStepCalculation(startTimeStepCalculation), endTimeStepCalculation(endTimeStepCalculation)
{
	this->tests = tests;
}