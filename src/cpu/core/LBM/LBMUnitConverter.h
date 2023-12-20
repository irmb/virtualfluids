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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup cpu_LBM LBM
//! \ingroup cpu_core core
//! \{
//! \author Soeren Freudiger
//=======================================================================================

#ifndef LBMUNITCONVERTER_H
#define LBMUNITCONVERTER_H

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <basics/utilities/UbException.h>

//! \brief A class provides converter for LB units <-> SI units
//! \details
//! \code
//! LBMUnitConverter conv(  100 /*L_World*/, 1484/*cs_water*/    , 1000/*rho_water*/
//!                         , 1000/*L_LB*/   , 1./srqt(3.)/*cs_Lb*/, 1/*rho_Lb*/ );
//! cout<<conv.toString()<<endl;
//!
//! cout<<"100m       = "<< 100  * conv.getFactorLentghWToLb()   << "dx    " << std::endl;
//! cout<<"1000dx     = "<< 1000 * conv.getFactorLentghLbToW()   << "m     " << std::endl;
//!
//! cout<<"25m/s      = "<< 25   * conv.getFactorVelocityWToLb() << "dx/dt " << std::endl;
//! cout<<"0.04 dx/dt = "<< 0.04 * conv.getFactorVelocityLbToW() << "m/s   " << std::endl;
//! \endcode
//! alternative
//! \code
//! LBMUnitConverter conv(, 100 /*L_World*/, LBMUnitConverter::WATER, 1000/*L_LB*/  );
//! \endcode

class LBMUnitConverter
{
public:
    enum WORLD_MATERIAL { WATER = 0, SEAWWATER = 1, AIR_20C = 2, OIL = 3 };

    LBMUnitConverter() = default;

    LBMUnitConverter(const real &refLengthWorld, const real &csWorld, const real &rhoWorld,
                     const real &refLengthLb, const real &csLb = 1.0 / std::sqrt(3.0), const real &rhoLb = 1.0)
    {
        this->init(refLengthWorld, csWorld, rhoWorld, csWorld, refLengthLb, rhoLb, csLb);
    }

    LBMUnitConverter(const real &refLengthWorld, WORLD_MATERIAL worldMaterial, const real &refLengthLb,
                     const real &csLb = 1.0 / std::sqrt(3.0), const real &rhoLb = 1.0)
    {
        real csWorld;
        real rhoWorld;

        if (worldMaterial == WATER) {
            csWorld  = 1484 /*m/s*/;
            rhoWorld = 1000 /*kg/m^3*/;
        } else if (worldMaterial == SEAWWATER) {
            csWorld  = 1500 /*m/s*/;
            rhoWorld = 1025 /*kg/m^3*/;
        } else if (worldMaterial == AIR_20C) {
            csWorld  = 343 /*m/s*/;
            rhoWorld = 1.290 /*kg/m^3*/;
        } else if (worldMaterial == OIL) {
            csWorld  = 1740 /*m/s*/;
            rhoWorld = 830 /*kg/m^3*/;
        } else
            throw UbException(UB_EXARGS, "unknown material");

        this->init(refLengthWorld, csWorld, rhoWorld, csWorld, refLengthLb, rhoLb, csLb);
    }

    virtual ~LBMUnitConverter() = default;

    real getRefRhoLb() { return refRhoLb; }

    real getFactorLentghLbToW() { return factorLengthLbToW; }
    real getFactorLentghWToLb() { return 1.0 / this->getFactorLentghLbToW(); }

    real getFactorTimeLbToW() { return factorTimeLbToW; }
    real getFactorTimeWToLb() { return 1.0 / this->getFactorTimeLbToW(); }

    real getFactorVelocityLbToW() { return factorLengthLbToW / factorTimeLbToW; }
    real getFactorVelocityWToLb() { return 1.0 / this->getFactorVelocityLbToW(); }

    real getFactorViscosityLbToW() { return factorLengthLbToW * factorLengthLbToW / factorTimeLbToW; }
    real getFactorViscosityWToLb() { return 1.0 / this->getFactorViscosityLbToW(); }

    real getFactorDensityLbToW() { return this->factorMassLbToW / std::pow(factorLengthLbToW, 3.0); }
    real getFactorDensityWToLb() { return 1.0 / this->getFactorDensityLbToW(); }

    real getFactorPressureLbToW(){ return this->factorMassLbToW / (factorLengthLbToW * factorTimeLbToW * factorTimeLbToW); }
    real getFactorPressureWToLb() { return 1.0 / this->getFactorPressureLbToW(); }

    real getFactorMassLbToW() { return this->factorMassLbToW; }
    real getFactorMassWToLb() { return 1.0 / this->getFactorMassLbToW(); }

    real getFactorForceLbToW() { return factorMassLbToW * factorLengthLbToW / (factorTimeLbToW * factorTimeLbToW); }
    real getFactorForceWToLb() { return 1.0 / this->getFactorForceLbToW(); }

    real getFactorTorqueLbToW() { return factorMassLbToW * factorLengthLbToW * factorLengthLbToW / (factorTimeLbToW * factorTimeLbToW);}
    real getFactorTorqueWToLb() { return 1.0 / this->getFactorTorqueLbToW(); }

    real getFactorAccLbToW() { return factorLengthLbToW / (factorTimeLbToW * factorTimeLbToW); }
    real getFactorAccWToLb() { return 1.0 / this->getFactorAccLbToW(); }

    real getFactorTimeLbToW(real deltaX) const { return factorTimeWithoutDx * deltaX; }


    /*==========================================================*/
    friend inline std::ostream &operator<<(std::ostream &os, LBMUnitConverter c)
    {
        os << c.toString();
        return os;
    }
    /*==========================================================*/
    std::string toString()
    {
        std::ostringstream out;
        out << "LB --> WORLD" << std::endl;
        out << " * lentgh 1[dx  ] = " << std::setw(12) << this->getFactorLentghLbToW() << " [m   ] " << std::endl;
        out << " * time   1[dt  ] = " << std::setw(12) << this->getFactorTimeLbToW() << " [s   ] " << std::endl;
        out << " * mass   1[mass] = " << std::setw(12) << this->getFactorMassLbToW() << " [kg  ] " << std::endl;
        out << std::endl;
        out << "WORLD --> LB" << std::endl;
        out << " * lentgh 1[m   ] = " << std::setw(12) << this->getFactorLentghWToLb() << " [dx  ] " << std::endl;
        out << " * time   1[s   ] = " << std::setw(12) << this->getFactorTimeWToLb() << " [dt  ] " << std::endl;
        out << " * mass   1[kg  ] = " << std::setw(12) << this->getFactorMassWToLb() << " [mass] " << std::endl;
        out << std::endl;
        out << "LB --> WORLD (combined units)" << std::endl;
        out << " * velocity     1 [dx/dt    ] = " << std::setw(12) << this->getFactorVelocityLbToW() << " [m/s      ]"
            << std::endl;
        out << " * density      1 [mass/dx^3] = " << std::setw(12) << this->getFactorDensityLbToW() << " [kg/m^3   ]"
            << std::endl;
        out << " * pressure     1 [F_lb/dx^2] = " << std::setw(12) << this->getFactorPressureLbToW() << " [N/m^2    ]"
            << std::endl;
        out << " * viscosity    1 [dx^2/dt  ] = " << std::setw(12) << this->getFactorViscosityLbToW() << " [m^2/s    ]"
            << std::endl;
        out << " * force        1 [F_lb     ] = " << std::setw(12) << this->getFactorForceLbToW() << " [N        ]"
            << std::endl;
        out << " * acceleration 1 [dx/dt^2  ] = " << std::setw(12) << this->getFactorAccLbToW() << " [m/s^2    ]"
            << std::endl;
        out << std::endl;
        out << "WORLD --> LB (combined units)" << std::endl;
        out << " * velocity     1 [m/s      ] = " << std::setw(12) << this->getFactorVelocityWToLb() << " [dx/dt    ]"
            << std::endl;
        out << " * density      1 [kg/m^3   ] = " << std::setw(12) << this->getFactorDensityWToLb() << " [mass/dx^3]"
            << std::endl;
        out << " * pressure     1 [N/m^2    ] = " << std::setw(12) << this->getFactorPressureWToLb() << " [F_lb/dx^2]"
            << std::endl;
        out << " * viscosity    1 [m^2/s    ] = " << std::setw(12) << this->getFactorViscosityWToLb() << " [dx^2/dt  ]"
            << std::endl;
        out << " * force        1 [N        ] = " << std::setw(12) << this->getFactorForceWToLb() << " [F_lb     ]"
            << std::endl;
        out << " * acceleration 1 [m/s^2    ] = " << std::setw(12) << this->getFactorAccWToLb() << " [dx/dt^2  ]"
            << std::endl;

        return out.str();
    }

    void init(const real &refLengthWorld, const real & /*csWorld*/, const real &rhoWorld, const real &vWorld,
              const real &refLengthLb, const real &rhoLb, const real &vLb)
    {
        factorLengthLbToW   = refLengthWorld / refLengthLb;
        factorTimeLbToW     = vLb / vWorld * factorLengthLbToW;
        factorMassLbToW     = rhoWorld / rhoLb * factorLengthLbToW * factorLengthLbToW * factorLengthLbToW;
        factorTimeWithoutDx = vLb / vWorld;
        this->refRhoLb      = rhoLb;
    }

protected:
    real factorLengthLbToW{ 1.0 };
    real factorTimeLbToW{ 1.0 };
    real factorMassLbToW{ 1.0 };
    real refRhoLb{ 1.0 };
    real factorTimeWithoutDx{ 0.0 };
};

#endif // LBMUNITCONVERTER_H

//! \}
