#include "TaylorGreenVortexUzLogFileDataImp.h"

std::shared_ptr<TaylorGreenVortexUzLogFileDataImp> TaylorGreenVortexUzLogFileDataImp::getNewInstance()
{
    return std::shared_ptr<TaylorGreenVortexUzLogFileDataImp>(new TaylorGreenVortexUzLogFileDataImp());
}

std::vector<int> TaylorGreenVortexUzLogFileDataImp::getL0()
{
    return l0;
}

std::vector<double> TaylorGreenVortexUzLogFileDataImp::getUz()
{
    return ux;
}

std::vector<double> TaylorGreenVortexUzLogFileDataImp::getAmplitude()
{
    return amp;
}

void TaylorGreenVortexUzLogFileDataImp::setL0(std::vector<int> l0)
{
    this->l0 = l0;
}

void TaylorGreenVortexUzLogFileDataImp::setUz(std::vector<double> ux)
{
    this->ux = ux;
}

void TaylorGreenVortexUzLogFileDataImp::setAmplitude(std::vector<double> amp)
{
    this->amp = amp;
}

TaylorGreenVortexUzLogFileDataImp::TaylorGreenVortexUzLogFileDataImp()
{
}

TaylorGreenVortexUzLogFileDataImp::~TaylorGreenVortexUzLogFileDataImp()
{
}
