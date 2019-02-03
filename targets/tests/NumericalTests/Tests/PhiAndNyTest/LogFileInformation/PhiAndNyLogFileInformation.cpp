#include "PhiAndNyLogFileInformation.h"

#include "Tests\PhiAndNyTest\PhiAndNyTest.h"
#include "Tests\PhiAndNyTest\PhiAndNyTestParameterStruct.h"

#include <iomanip>
#include <sstream>


std::shared_ptr<PhiAndNyInformation> PhiAndNyInformation::getNewInstance(std::shared_ptr<PhiAndNyTestParameterStruct> testPara)
{
	return std::shared_ptr<PhiAndNyInformation>(new PhiAndNyInformation(testPara));
}

std::string PhiAndNyInformation::getOutput()
{
	std::ostringstream headName;
	headName << testGroups.at(0).at(0)->getSimulationName() <<" Phi And Ny Test";
	makeCenterHead(headName.str());

	oss << "StartTimeStepCalculation=" << startTimeStepCalculation << std::endl;
	oss << "EndTimeStepCalculation=" << endTimeStepCalculation << std::endl;
	oss << std::endl;

	for (int i = 0; i < testGroups.size(); i++) {
		fillMyData(testGroups.at(i));
		for (int j = 0; j < ny.size(); j++)
			oss << "Ny_" << lxForErase.at(j) << "_" << dataToCalc.at(j) << "=" << ny.at(j) << std::endl;
		oss << std::endl;
		for (int j = 0; j < nyDiff.size(); j++)
			oss << "NyDiff_" << lxForErase.at(j) << "_" << dataToCalc.at(j) << "=" << nyDiff.at(j) << std::endl;
		oss << std::endl;
		for (int j = 0; j < orderOfAccuracyNyDiff.size(); j++)
			oss << "OrderOfAccuracy_NyDiff_" << lx.at(2*j) << "_" << lx.at(2*j+1) << "_" << dataToCalc.at(j) << "=" << orderOfAccuracyNyDiff.at(j) << std::endl;
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

void PhiAndNyInformation::addTestGroup(std::vector<std::shared_ptr<PhiAndNyTest> > tests)
{
	testGroups.push_back(tests);
}

void PhiAndNyInformation::fillMyData(std::vector<std::shared_ptr<PhiAndNyTest> > testGroup)
{
	lxForErase.resize(0);
	lx.resize(0);
	ny.resize(0);
	nyDiff.resize(0);
	phiDiff.resize(0);
	orderOfAccuracyNyDiff.resize(0);
	orderOfAccuracyPhiDiff.resize(0);
	dataToCalc.resize(0);
	for (int i = 0; i < testGroup.size(); i++) {
		std::vector<int> myLx = testGroup.at(i)->getLx();
		std::vector<double> myNy = testGroup.at(i)->getNy();
		std::vector<double> myNyDiff = testGroup.at(i)->getNyDiff();
		std::vector<double> myPhiDiff = testGroup.at(i)->getPhiDiff();

		lx.insert(lx.end(), myLx.begin(), myLx.end());
		lxForErase.insert(lxForErase.end(), myLx.begin(), myLx.end());
		ny.insert(ny.end(), myNy.begin(), myNy.end());
		nyDiff.insert(nyDiff.end(), myNyDiff.begin(), myNyDiff.end());
		phiDiff.insert(phiDiff.end(), myPhiDiff.begin(), myPhiDiff.end());
		orderOfAccuracyNyDiff.push_back(testGroup.at(i)->getOrderOfAccuracyNyDiff());
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
			ny.erase(ny.begin() + i);
			nyDiff.erase(nyDiff.begin() + i);
			phiDiff.erase(phiDiff.begin() + i);
			lxForErase.erase(lxForErase.begin() + i);
		}
	}

	
}

PhiAndNyInformation::PhiAndNyInformation(std::shared_ptr<PhiAndNyTestParameterStruct> testPara)
{
	startTimeStepCalculation = testPara->startTimeStepCalculation;
	endTimeStepCalculation = testPara->endTimeStepCalculation;
}