#include "Thixotropy.h"

SPtr<Thixotropy> Thixotropy::instance = SPtr<Thixotropy>();
LBMReal Thixotropy::tau0 = 0;
LBMReal Thixotropy::k = 0;
LBMReal Thixotropy::n = 1;
LBMReal Thixotropy::omegaMin = 0;
LBMReal Thixotropy::beta = 0;
LBMReal Thixotropy::c = 0;
LBMReal Thixotropy::mu0 = 0;

//////////////////////////////////////////////////////////////////////////
SPtr<Thixotropy> Thixotropy::getInstance()
{
   if (!instance)
      instance = SPtr<Thixotropy>(new Thixotropy());
   return instance;
}

void Thixotropy::setYieldStress(LBMReal yieldStress)
{
	tau0 = yieldStress;
}
LBMReal Thixotropy::getYieldStress() const
{
	return tau0;
}
void Thixotropy::setViscosityParameter(LBMReal kParameter)
{
	k = kParameter;
}
LBMReal Thixotropy::getViscosityParameter() const
{
	return k;
}
void Thixotropy::setPowerIndex(LBMReal index)
{
	n = index;
}
LBMReal Thixotropy::getPowerIndex() const
{
	return n;
}

void Thixotropy::setOmegaMin(LBMReal omega)
{
	omegaMin = omega;
}
LBMReal Thixotropy::getOmegaMin() const
{
	return omegaMin;
}

void Thixotropy::setBeta(LBMReal PowellEyringBeta)
{
	beta = PowellEyringBeta;
}

LBMReal Thixotropy::getBeta() const
{
	return beta;
}

void Thixotropy::setC(LBMReal PowellEyringC)
{
	c = PowellEyringC;
}

LBMReal Thixotropy::getC() const
{
	return c;
}

void Thixotropy::setMu0(LBMReal mu)
{
	mu0 = mu;
}

LBMReal Thixotropy::getMu0() const
{
	return mu0;
}

Thixotropy::Thixotropy()
{
}