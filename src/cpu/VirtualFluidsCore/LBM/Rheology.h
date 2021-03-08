//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file Rheology.h
//! \ingroup LBM
//! \author Konstantin Kutscher, Martin Geier
//=======================================================================================

#ifndef Rheology_H
#define Rheology_H

#include <PointerDefinitions.h>
#include <LBMSystem.h>
#include <UbMath.h>
#include <math.h> 

class Rheology
{
public:
	Rheology(Rheology const&) = delete;
	Rheology& operator=(Rheology const&) = delete;
	static SPtr<Rheology> getInstance();
	void setYieldStress(LBMReal tau0);
	LBMReal getYieldStress() const;
	
	void setViscosityParameter(LBMReal k);
	LBMReal getViscosityParameter() const;

	void setPowerIndex(LBMReal n);
	LBMReal getPowerIndex() const;

	void setOmegaMin(LBMReal omegaMin);
	LBMReal getOmegaMin() const;

	void setBeta(LBMReal PowellEyringBeta);
	LBMReal getBeta() const;

	void setC(LBMReal PowellEyringC);
	LBMReal getC() const;

	void setMu0(LBMReal mu);
	LBMReal getMu0() const;

	static LBMReal getBinghamCollFactorOld(LBMReal omegaInf, LBMReal shearRate, LBMReal drho);
	static LBMReal getBinghamCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho);
	static LBMReal getHerschelBulkleyCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho);
	static LBMReal getHerschelBulkleyCollFactorBackward(LBMReal shearRate, LBMReal drho);
	static LBMReal getPowellEyringCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho);
private:
	Rheology();
	
	static SPtr<Rheology> instance;

	static LBMReal tau0;
	static LBMReal k;
	static LBMReal n;
	static LBMReal omegaMin;
	static LBMReal beta;
	static LBMReal c;
	static LBMReal mu0;
};

//////////////////////////////////////////////////////////////////////////
inline LBMReal Rheology::getBinghamCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho)
{
	LBMReal cs2 = UbMath::one_over_sqrt3 * UbMath::one_over_sqrt3;
	LBMReal rho = UbMath::one + drho;
	//analytical solution
	LBMReal omega = omegaInf * (UbMath::one - (omegaInf * tau0) / (shearRate * cs2 * rho + UbMath::Epsilon<LBMReal>::val()));
	
	//LBMReal omega = cs2 * cs2 * shearRate * shearRate * omegaInf * rho * rho / (cs2 * cs2 * shearRate * shearRate * rho * rho + cs2 * shearRate * omegaInf * rho * tau0+omegaInf*omegaInf*tau0*tau0);
	
	// LBMReal a = omegaInf * tau0 / (cs2 * shearRate * rho);
	//10 iterations
	//LBMReal omega = omegaInf / (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a))))))))));
	
	//20 iterations
	////LBMReal omega = omegaInf / (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a * (1 + a))))))))))))))))))));
	//
	//LBMReal omega = omegaInf*cs2 * shearRate * rho / (cs2 * shearRate * rho + omegaInf * tau0);
	//LBMReal shearRateNew = shearRate * (omega / omegaInf);

	//for (int i = 0; i < 20; i++)
	//{
	//	omega = omegaInf * cs2 * shearRateNew * rho / (cs2 * shearRateNew * rho + omegaInf * tau0);
	//	shearRateNew = shearRate * (omega / omegaInf);
	//}
	//omega = omegaInf * cs2 * shearRateNew * rho / (cs2 * shearRateNew * rho + omegaInf * tau0);
	//
	////if (omega < 0.2)
	////	omega = 0.2;
	return omega;
}

inline LBMReal Rheology::getBinghamCollFactorOld(LBMReal omegaInf, LBMReal shearRate, LBMReal drho)
{
	const LBMReal cs2 = UbMath::c1o3; // UbMath::one_over_sqrt3* UbMath::one_over_sqrt3;
	LBMReal rho = UbMath::one + drho;

	if (rho * cs2 * (UbMath::c1 / omegaInf - UbMath::c1o2) * shearRate < tau0)
		return 0.0;
	else
		return omegaInf;
}
//////////////////////////////////////////////////////////////////////////
inline LBMReal Rheology::getHerschelBulkleyCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho)
{
	LBMReal cs2 = UbMath::one_over_sqrt3 * UbMath::one_over_sqrt3;
	LBMReal rho = UbMath::one + drho;
	LBMReal gammaDot = shearRate;
	LBMReal omega = omegaInf;
	LBMReal epsilon = 1;
	LBMReal gammaDotPowN = std::pow(gammaDot, n);

	while (epsilon > 1e-10)
	{
		LBMReal omegaOld = omega;
		LBMReal omegaByOmegaInfPowN = std::pow(omega / omegaInf, n);/*
		LBMReal gammaDotPowOneMinusN = std::pow(gammaDot,1- n);
		LBMReal omegaByOmegaInfPowOneMinusN = std::pow(omega / omegaInf, 1-n);
		LBMReal numeratorA = (2.0* k *  omegaInf + cs2 * gammaDotPowOneMinusN * omegaByOmegaInfPowOneMinusN *omegaInf* rho );
		LBMReal numeratorB = ( cs2 * gammaDot * ( - 2.0) * rho + 2.0 * omegaInf * tau0);
		LBMReal denominatorA = (2.0 * k * n * omegaInf + cs2 * gammaDot * rho * omegaInf* gammaDotPowOneMinusN * omegaByOmegaInfPowOneMinusN) + UbMath::Epsilon<LBMReal>::val();
		LBMReal denominatorB = (2.0 * k * n * gammaDotPowN * omegaByOmegaInfPowN * omegaInf + cs2 * gammaDot * rho * omega) + UbMath::Epsilon<LBMReal>::val();
		omega = omega - omega *( numeratorA / denominatorA+ numeratorB / denominatorB);*/
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
inline LBMReal Rheology::getHerschelBulkleyCollFactorBackward(LBMReal shearRate, LBMReal drho)
{
	LBMReal rho = UbMath::one + drho;
	LBMReal gamma = shearRate + UbMath::Epsilon<LBMReal>::val();
	LBMReal cs2 = UbMath::one_over_sqrt3 * UbMath::one_over_sqrt3;

	return 1.0 / ((tau0 + k * std::pow(gamma, n)) / (cs2 * rho * gamma) + UbMath::c1o2);
}
//////////////////////////////////////////////////////////////////////////
inline LBMReal Rheology::getPowellEyringCollFactor(LBMReal omegaInf, LBMReal shearRate, LBMReal drho)
{
	using namespace UbMath;
	LBMReal cs2 = c1o3; // UbMath::one_over_sqrt3* UbMath::one_over_sqrt3;
	LBMReal rho = c1 + drho;
	LBMReal gammaDot = shearRate;
	LBMReal omega = omegaInf;
	LBMReal epsilon = 1;

	while (epsilon > 1e-10)
	{
		LBMReal omegaOld = omega;
		epsilon = std::abs(omega - omegaOld);

		LBMReal numerator = c*sqrt(c1+(gammaDot*gammaDot*omega*omega)/(c*c*omegaInf*omegaInf))*(beta*(c2*gammaDot*mu0*omega+cs2*gammaDot*(omega-c2)*rho+c2*omegaInf*tau0)+c2*omegaInf*(asinh((gammaDot*omega)/(c*omegaInf))));

		LBMReal denominator = gammaDot*(c2+beta*c*sqrt(c1+(gammaDot*gammaDot*omega*omega)/(c*c*omegaInf*omegaInf))*(c2*mu0+cs2*rho)) + UbMath::Epsilon<LBMReal>::val();

		omega = omega - numerator / denominator;

		omega = (omega < UbMath::zeroReal) ? UbMath::c1o2 * omegaOld : omega;
	}

	return omega;
}
#endif
