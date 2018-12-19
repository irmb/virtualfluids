#include "PhiAndNuLogFileInformation.h"

#include "Tests\PhiAndNuTest\PhiAndNuTest.h"

#include <iomanip>
#include <sstream>


std::shared_ptr<PhiAndNuInformation> PhiAndNuInformation::getNewInstance(unsigned int startTimeStepCalculation, unsigned int endTimeStepCalculation)
{
	return std::shared_ptr<PhiAndNuInformation>(new PhiAndNuInformation(startTimeStepCalculation, endTimeStepCalculation));
}

std::string PhiAndNuInformation::getOutput()
{
	std::ostringstream headName;
	headName << testGroups.at(0).at(0)->getSimulationName() <<" Phi And Nu Test";
	makeCenterHead(headName.str());

	oss << "StartTimeStepCalculation=" << startTimeStepCalculation << std::endl;
	oss << "EndTimeStepCalculation=" << endTimeStepCalculation << std::endl;
	oss << std::endl;

	for (int i = 0; i < testGroups.size(); i++) {
		fillMyData(testGroups.at(i));
		for (int j = 0; j < nu.size(); j++)
			oss << "Nu_" << lxForErase.at(j) << "_" << dataToCalc.at(j) << "=" << nu.at(j) << std::endl;
		oss << std::endl;
		for (int j = 0; j < nuDiff.size(); j++)
			oss << "NuDiff_" << lxForErase.at(j) << "_" << dataToCalc.at(j) << "=" << nuDiff.at(j) << std::endl;
		oss << std::endl;
		for (int j = 0; j < orderOfAccuracyNuDiff.size(); j++)
			oss << "OrderOfAccuracy_NuDiff_" << lx.at(2*j) << "_" << lx.at(2*j+1) << "_" << dataToCalc.at(j) << "=" << orderOfAccuracyNuDiff.at(j) << std::endl;
		oss << std::endl;
		for (int j = 0; j < phiDiff.size(); j++)
			oss << "PhiDiff_" << lxForErase.at(j) << "_" << dataToCalc.at(j) << "=" << phiDiff.at(j) << std::endl;
		oss << std::endl;
		for (int j = 0; j < orderOfAccuracyPhiDiff.size(); j++)
			oss << "OrderOfAccuracy_PhiDiff_" << lx.at(2*j) << "_" << lx.at(2*j+1) << "_" << dataToCalc.at(j) << "=" << orderOfAccuracyPhiDiff.at(j) << std::endl;
		oss << std::endl;
	}

	return oss.str();
}

void PhiAndNuInformation::addTestGroup(std::vector<std::shared_ptr<PhiAndNuTest>> tests)
{
	testGroups.push_back(tests);
}

void PhiAndNuInformation::fillMyData(std::vector<std::shared_ptr<PhiAndNuTest>> testGroup)
{
	lxForErase.resize(0);
	lx.resize(0);
	nu.resize(0);
	nuDiff.resize(0);
	phiDiff.resize(0);
	orderOfAccuracyNuDiff.resize(0);
	orderOfAccuracyPhiDiff.resize(0);
	dataToCalc.resize(0);
	for (int i = 0; i < testGroup.size(); i++) {
		std::vector<int> myLx = testGroup.at(i)->getLx();
		std::vector<double> myNu = testGroup.at(i)->getNu();
		std::vector<double> myNuDiff = testGroup.at(i)->getNuDiff();
		std::vector<double> myPhiDiff = testGroup.at(i)->getPhiDiff();

		lx.insert(lx.end(), myLx.begin(), myLx.end());
		lxForErase.insert(lxForErase.end(), myLx.begin(), myLx.end());
		nu.insert(nu.end(), myNu.begin(), myNu.end());
		nuDiff.insert(nuDiff.end(), myNuDiff.begin(), myNuDiff.end());
		phiDiff.insert(phiDiff.end(), myPhiDiff.begin(), myPhiDiff.end());
		orderOfAccuracyNuDiff.push_back(testGroup.at(i)->getOrderOfAccuracyNuDiff());
		orderOfAccuracyPhiDiff.push_back(testGroup.at(i)->getOrderOfAccuracyPhiDiff());
		dataToCalc.push_back(testGroup.at(i)->getDataToCalculate());
		dataToCalc.push_back(testGroup.at(i)->getDataToCalculate());
	}

	for (int i = 0; i < lxForErase.size(); i++) 
		for (int j = i + 1; j < lxForErase.size(); j++)
			if (lxForErase.at(i) == lxForErase.at(j))
				lxForErase.at(j) = -1;
	
	for (int i = lxForErase.size() - 1; i >= 0; i--) {
		if (lxForErase.at(i) == -1) {
			nu.erase(nu.begin() + i);
			nuDiff.erase(nuDiff.begin() + i);
			phiDiff.erase(phiDiff.begin() + i);
			lxForErase.erase(lxForErase.begin() + i);
		}
	}

	
}

PhiAndNuInformation::PhiAndNuInformation(unsigned int startTimeStepCalculation, unsigned int endTimeStepCalculation) : startTimeStepCalculation(startTimeStepCalculation), endTimeStepCalculation(endTimeStepCalculation)
{

}