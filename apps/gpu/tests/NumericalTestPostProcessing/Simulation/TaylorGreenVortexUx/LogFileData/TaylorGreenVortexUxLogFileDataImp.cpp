#include "TaylorGreenVortexUxLogFileDataImp.h"

std::shared_ptr<TaylorGreenVortexUxLogFileDataImp> TaylorGreenVortexUxLogFileDataImp::getNewInstance()
{
	return std::shared_ptr<TaylorGreenVortexUxLogFileDataImp>(new TaylorGreenVortexUxLogFileDataImp());
}

std::vector<int> TaylorGreenVortexUxLogFileDataImp::getL0()
{
	return l0;
}

std::vector<double> TaylorGreenVortexUxLogFileDataImp::getUx()
{
	return ux;
}

std::vector<double> TaylorGreenVortexUxLogFileDataImp::getAmplitude()
{
	return amp;
}

void TaylorGreenVortexUxLogFileDataImp::setL0(std::vector<int> l0)
{
	this->l0 = l0;
}

void TaylorGreenVortexUxLogFileDataImp::setUx(std::vector<double> ux)
{
	this->ux = ux;
}

void TaylorGreenVortexUxLogFileDataImp::setAmplitude(std::vector<double> amp)
{
	this->amp = amp;
}

TaylorGreenVortexUxLogFileDataImp::TaylorGreenVortexUxLogFileDataImp()
{
}

TaylorGreenVortexUxLogFileDataImp::~TaylorGreenVortexUxLogFileDataImp()
{
}
