#include "L2NormLogFileDataImp.h"

std::shared_ptr<L2NormLogFileDataImp> L2NormLogFileDataImp::getNewInstance()
{
	return std::shared_ptr<L2NormLogFileDataImp>(new L2NormLogFileDataImp());
}

std::vector<double> L2NormLogFileDataImp::getBasicGridLengths()
{
	return basicGridLengths;
}

std::string L2NormLogFileDataImp::getDataToCalc()
{
	return dataToCalc;
}

std::string L2NormLogFileDataImp::getNormalizeData()
{
	return normalizeData;
}

int L2NormLogFileDataImp::getBasicTimeStep()
{
	return basicTimeStep;
}

int L2NormLogFileDataImp::getDivergentTimeStep()
{
	return divergentTimeStep;
}

std::vector<double> L2NormLogFileDataImp::getL2NormForBasicTimeStep()
{
	return l2NormForBasicTimeStep;
}

std::vector<double> L2NormLogFileDataImp::getL2NormForDivergentTimeStep()
{
	return l2NormForDivergentTimeStep;
}

std::vector<double> L2NormLogFileDataImp::getL2NormDiff()
{
	return l2NormDiff;
}

void L2NormLogFileDataImp::setBasicGridLengths(std::vector<double> basicGridLengths)
{
	this->basicGridLengths = basicGridLengths;
}

void L2NormLogFileDataImp::setDataToCalc(std::string dataToCalc)
{
	this->dataToCalc = dataToCalc;
}

void L2NormLogFileDataImp::setNormalizeData(std::string normalizeData)
{
	this->normalizeData = normalizeData;
}

void L2NormLogFileDataImp::setBasicTimeStep(int basicTimeStep)
{
	this->basicTimeStep = basicTimeStep;
}

void L2NormLogFileDataImp::setDivergentTimeStep(int divergentTimeStep)
{
	this->divergentTimeStep = divergentTimeStep;
}

void L2NormLogFileDataImp::setL2NormForBasicTimeStep(std::vector<double> l2NormForBasicTimeStep)
{
	this->l2NormForBasicTimeStep = l2NormForBasicTimeStep;
}

void L2NormLogFileDataImp::setL2NormForDivergentTimeStep(std::vector<double> l2NormForDivergentTimeStep)
{
	this->l2NormForDivergentTimeStep = l2NormForDivergentTimeStep;
}

void L2NormLogFileDataImp::setL2NormDiff(std::vector<double> l2NormDiff)
{
	this->l2NormDiff = l2NormDiff;
}

L2NormLogFileDataImp::L2NormLogFileDataImp()
{
}

L2NormLogFileDataImp::~L2NormLogFileDataImp()
{
}
