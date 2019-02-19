#include "PhiLogFileDataImp.h"

std::shared_ptr<PhiLogFileDataImp> PhiLogFileDataImp::getNewInstance()
{
	return std::shared_ptr<PhiLogFileDataImp>(new PhiLogFileDataImp());
}

std::vector<double> PhiLogFileDataImp::getBasicGridLengths()
{
	return basicGridLengths;
}

int PhiLogFileDataImp::getStartTimeStepCalculation()
{
	return startTimeStepCalculation;
}

int PhiLogFileDataImp::getEndTimeStepCalculation()
{
	return endTimeStepCalculation;
}

std::string PhiLogFileDataImp::getDataToCalc()
{
	return dataToCalc;
}

std::vector<double> PhiLogFileDataImp::getPhiDiff()
{
	return phiDiff;
}

std::vector<std::vector<double>> PhiLogFileDataImp::getOrderOfAccuracy()
{
	return orderOfAccuracy;
}

void PhiLogFileDataImp::setBasicGridLengths(std::vector<double> basicGridLengths)
{
	this->basicGridLengths = basicGridLengths;
}

void PhiLogFileDataImp::setStartTimeStepCalculation(int startTimeStepCalculation)
{
	this->startTimeStepCalculation = startTimeStepCalculation;
}

void PhiLogFileDataImp::setEndTimeStepCalculation(int endTimeStepCalculation)
{
	this->endTimeStepCalculation = endTimeStepCalculation;
}

void PhiLogFileDataImp::setDataToCalc(std::string dataToCalc)
{
	this->dataToCalc = dataToCalc;
}

void PhiLogFileDataImp::setPhiDiff(std::vector<double> phiDiff)
{
	this->phiDiff = phiDiff;
}

void PhiLogFileDataImp::setOrderOfAccuracy(std::vector<std::vector<double>> orderOfAccuracy)
{
	this->orderOfAccuracy = orderOfAccuracy;
}

PhiLogFileDataImp::PhiLogFileDataImp()
{
}
