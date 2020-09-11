#ifndef Thixotropy_H
#define Thixotropy_H

#include <PointerDefinitions.h>
#include <LBMSystem.h>
#include <UbMath.h>

class Thixotropy
{
public:
	Thixotropy(Thixotropy const&) = delete;
	Thixotropy& operator=(Thixotropy const&) = delete;
	static SPtr<Thixotropy> getInstance();
	void setYieldStress(LBMReal tau0);
	LBMReal getYieldStress() const;
	
	void setViscosityParameter(LBMReal k);
	LBMReal getViscosityParameter() const;

	void setPowerIndex(LBMReal n);
	LBMReal getPowerIndex() const;

	void setOmegaMin(LBMReal omegaMin);
	LBMReal getOmegaMin() const;

	static LBMReal getBinghamCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho);
	static LBMReal getHerschelBulkleyCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho);
	static LBMReal getHerschelBulkleyCollFactorBackward(LBMReal shearRate, LBMReal drho);
private:
	Thixotropy();
	
	static SPtr<Thixotropy> instance;

	static LBMReal tau0;
	static LBMReal k;
	static LBMReal n;
	static LBMReal omegaMin;
};

//////////////////////////////////////////////////////////////////////////
inline LBMReal Thixotropy::getBinghamCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho)
{
	LBMReal cs2 = UbMath::one_over_sqrt3 * UbMath::one_over_sqrt3;
	LBMReal rho = UbMath::one + drho;
	LBMReal omega = omegaInf * (UbMath::one - (omegaInf * tau0) / (shearRate * cs2 * rho + UbMath::Epsilon<LBMReal>::val()));
	return omega;
}
//////////////////////////////////////////////////////////////////////////
inline LBMReal Thixotropy::getHerschelBulkleyCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho)
{
	LBMReal cs2 = UbMath::one_over_sqrt3 * UbMath::one_over_sqrt3;
	LBMReal rho = UbMath::one + drho;
	LBMReal gammaDot = shearRate;
	LBMReal omega = omegaInf;
	LBMReal epsilon = 1;

	while (epsilon > 1e-10)
	{
		LBMReal omegaOld = omega;
		LBMReal gammaDotPowN = std::pow(gammaDot, n);
		LBMReal omegaByOmegaInfPowN = std::pow(omega / omegaInf, n);
		LBMReal numerator = (2.0 * gammaDotPowN * k * omegaByOmegaInfPowN * omegaInf + cs2 * gammaDot * (omega - 2.0) * rho + 2.0 * omegaInf * tau0);
		LBMReal denominator = (2.0 * k * n * gammaDotPowN * omegaByOmegaInfPowN * omegaInf + cs2 * gammaDot * rho * omega) + UbMath::Epsilon<LBMReal>::val();
		omega = omega - omega * numerator / denominator;
		omega = (omega < UbMath::zeroReal) ? UbMath::c1o2 * omegaOld : omega;
        //omega = (omega < omegaMin) ? UbMath::c1o2 * (omegaOld-omegaMin)+omegaMin : omega;
		epsilon = std::abs(omega - omegaOld);
	}

	return omega;
}
//////////////////////////////////////////////////////////////////////////
inline LBMReal Thixotropy::getHerschelBulkleyCollFactorBackward(LBMReal shearRate, LBMReal drho)
{
	LBMReal rho = UbMath::one + drho;
	LBMReal gamma = shearRate + UbMath::Epsilon<LBMReal>::val();
	LBMReal cs2 = UbMath::one_over_sqrt3 * UbMath::one_over_sqrt3;

	return 1.0 / ((tau0 + k * std::pow(gamma, n)) / (cs2 * rho * gamma) + UbMath::c1o2);
}
#endif
