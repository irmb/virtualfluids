#include "ShearWaveLogFileDataImp.h"

std::shared_ptr<ShearWaveLogFileDataImp> ShearWaveLogFileDataImp::getNewInstance()
{
	return std::shared_ptr<ShearWaveLogFileDataImp>(new ShearWaveLogFileDataImp());
}

std::vector<int> ShearWaveLogFileDataImp::getL0()
{
	return l0;
}

std::vector<double> ShearWaveLogFileDataImp::getUx()
{
	return ux;
}

std::vector<double> ShearWaveLogFileDataImp::getUz()
{
	return uz;
}

void ShearWaveLogFileDataImp::setL0(std::vector<int> l0)
{
	this->l0 = l0;
}

void ShearWaveLogFileDataImp::setUx(std::vector<double> ux)
{
	this->ux = ux;
}

void ShearWaveLogFileDataImp::setUz(std::vector<double> uz)
{
	this->uz = uz;
}

ShearWaveLogFileDataImp::ShearWaveLogFileDataImp()
{
}

ShearWaveLogFileDataImp::~ShearWaveLogFileDataImp()
{
}
