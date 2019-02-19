#include "NyLogFileDataImp.h"

std::shared_ptr<NyLogFileDataImp> NyLogFileDataImp::getNewInstance()
{
	return std::shared_ptr<NyLogFileDataImp>(new NyLogFileDataImp());
}

std::vector<double> NyLogFileDataImp::getBasicGridLengths()
{
	return basicGridLengths;
}

int NyLogFileDataImp::getStartTimeStepCalculation()
{
	return startTimeStepCalculation;
}

int NyLogFileDataImp::getEndTimeStepCalculation()
{
	return endTimeStepCalculation;
}

std::string NyLogFileDataImp::getDataToCalc()
{
	return dataToCalc;
}

std::vector<double> NyLogFileDataImp::getNy()
{
	return ny;
}

std::vector<double> NyLogFileDataImp::getNyDiff()
{
	return nyDiff;
}

std::vector<std::vector<double>> NyLogFileDataImp::getOrderOfAccuracy()
{
	return orderOfAccuracy;
}

void NyLogFileDataImp::setBasicGridLengths(std::vector<double> basicGridLengths)
{
	this->basicGridLengths = basicGridLengths;
}

void NyLogFileDataImp::setStartTimeStepCalculation(int startTimeStepCalculation)
{
	this->startTimeStepCalculation = startTimeStepCalculation;
}

void NyLogFileDataImp::setEndTimeStepCalculation(int endTimeStepCalculation)
{
	this->endTimeStepCalculation = endTimeStepCalculation;
}

void NyLogFileDataImp::setDataToCalcPhiAndNu(std::string dataToCalc)
{
	this->dataToCalc = dataToCalc;
}

void NyLogFileDataImp::setNy(std::vector<double> ny)
{
	this->ny = ny;
}

void NyLogFileDataImp::setNyDiff(std::vector<double> nyDiff)
{
	this->nyDiff = nyDiff;
}

void NyLogFileDataImp::setOrderOfAccuracy(std::vector<std::vector<double>> orderOfAccuracy)
{
	this->orderOfAccuracy = orderOfAccuracy;
}

NyLogFileDataImp::NyLogFileDataImp()
{
}
