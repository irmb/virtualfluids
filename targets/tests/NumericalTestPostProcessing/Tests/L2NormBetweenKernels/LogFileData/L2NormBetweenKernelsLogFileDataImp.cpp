#include "L2NormBetweenKernelsLogFileDataImp.h"

std::shared_ptr<L2NormBetweenKernelsLogFileDataImp> L2NormBetweenKernelsLogFileDataImp::getNewInstance()
{
	return std::shared_ptr<L2NormBetweenKernelsLogFileDataImp>(new L2NormBetweenKernelsLogFileDataImp());
}

std::vector<double> L2NormBetweenKernelsLogFileDataImp::getBasicGridLengths()
{
	return basicGridLengths;
}

std::string L2NormBetweenKernelsLogFileDataImp::getBasicKernel()
{
	return basicKernel;
}

std::string L2NormBetweenKernelsLogFileDataImp::getDataToCalculate()
{
	return dataToCalc;
}

int L2NormBetweenKernelsLogFileDataImp::getTimeStep()
{
	return timeStep;
}

std::vector<double> L2NormBetweenKernelsLogFileDataImp::getL2NormForBasicKernel()
{
	return l2NormForBasicKernel;
}

std::vector<double> L2NormBetweenKernelsLogFileDataImp::getL2NormForDivergentKernel()
{
	return l2NormForDivergentKernel;
}

std::vector<double> L2NormBetweenKernelsLogFileDataImp::getL2NormBetweenKernels()
{
	return l2NormBetweenKernels;
}

std::string L2NormBetweenKernelsLogFileDataImp::getNormalizeData()
{
	return normalizeData;
}

void L2NormBetweenKernelsLogFileDataImp::setBasicGridLengths(std::vector<double> basicGridLengths)
{
	this->basicGridLengths = basicGridLengths;
}

void L2NormBetweenKernelsLogFileDataImp::setBasicKernel(std::string basicKernel)
{
	this->basicKernel = basicKernel;
}

void L2NormBetweenKernelsLogFileDataImp::setDataToCalculate(std::string dataToCalc)
{
	this->dataToCalc = dataToCalc;
}

void L2NormBetweenKernelsLogFileDataImp::setTimeStep(int timeStep)
{
	this->timeStep = timeStep;
}

void L2NormBetweenKernelsLogFileDataImp::setL2NormForBasicKernel(std::vector<double> l2Norm)
{
	this->l2NormForBasicKernel = l2Norm;
}

void L2NormBetweenKernelsLogFileDataImp::setL2NormForDivergentKernel(std::vector<double> l2Norm)
{
	this->l2NormForDivergentKernel = l2Norm;
}

void L2NormBetweenKernelsLogFileDataImp::setL2NormBetweenKernels(std::vector<double> l2Norm)
{
	this->l2NormBetweenKernels = l2Norm;
}

void L2NormBetweenKernelsLogFileDataImp::setNormalizeData(std::string normalizeData)
{
	this->normalizeData = normalizeData;
}

L2NormBetweenKernelsLogFileDataImp::L2NormBetweenKernelsLogFileDataImp()
{
}

L2NormBetweenKernelsLogFileDataImp::~L2NormBetweenKernelsLogFileDataImp()
{
}
