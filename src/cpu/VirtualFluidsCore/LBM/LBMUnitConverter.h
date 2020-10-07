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
//! \file LBMUnitConverter.h
//! \ingroup LBM
//! \author Soeren Freudiger
//=======================================================================================

#ifndef LBMUNITCONVERTER_H
#define LBMUNITCONVERTER_H

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cmath>

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

   enum WORLD_MATERIAL { WATER  = 0, SEAWWATER  = 1, AIR_20C  = 2, OIL  = 3  }; 

   LBMUnitConverter()  
                         
   {

   }

   LBMUnitConverter(   const double& refLengthWorld, const double& csWorld, const double& rhoWorld
      , const double& refLengthLb   , const double& csLb = 1.0/std::sqrt(3.0)  , const double& rhoLb = 1.0   )
   {
      this->init(  refLengthWorld, csWorld, rhoWorld, csWorld, refLengthLb, rhoLb, csLb  );

   }

   LBMUnitConverter(  const double& refLengthWorld, WORLD_MATERIAL worldMaterial
      , const double& refLengthLb   , const double& csLb = 1.0/std::sqrt(3.0) , const double& rhoLb = 1.0    )
   {
      double csWorld;
      double rhoWorld;  

      if     ( worldMaterial == WATER    ) { csWorld = 1484/*m/s*/; rhoWorld =  1000/*kg/m^3*/;  }
      else if( worldMaterial == SEAWWATER) { csWorld = 1500/*m/s*/; rhoWorld =  1025/*kg/m^3*/;  }
      else if( worldMaterial == AIR_20C  ) { csWorld =  343/*m/s*/; rhoWorld = 1.290/*kg/m^3*/;  }
      else if( worldMaterial == OIL      ) { csWorld = 1740/*m/s*/; rhoWorld =  830/*kg/m^3*/;   }
      else                                  throw UbException(UB_EXARGS,"unknown material");

      this->init(  refLengthWorld, csWorld, rhoWorld, csWorld, refLengthLb, rhoLb, csLb  );

   }

   virtual ~LBMUnitConverter() = default;

   double  getRefRhoLb()             { return refRhoLb; }

   double  getFactorLentghLbToW()    { return factorLengthLbToW;                                                       }
   double  getFactorLentghWToLb()    { return 1.0/this->getFactorLentghLbToW();                                        }

   double  getFactorTimeLbToW()      { return factorTimeLbToW;                                                         }
   double  getFactorTimeWToLb()      { return 1.0/this->getFactorTimeLbToW();                                          }

   double  getFactorVelocityLbToW()  { return factorLengthLbToW/factorTimeLbToW;                                       }
   double  getFactorVelocityWToLb()  { return 1.0/this->getFactorVelocityLbToW();                                      }

   double  getFactorViscosityLbToW() { return factorLengthLbToW*factorLengthLbToW/factorTimeLbToW;                     }
   double  getFactorViscosityWToLb() { return 1.0/this->getFactorViscosityLbToW();                                     }

   double  getFactorDensityLbToW()   { return this->factorMassLbToW/std::pow(factorLengthLbToW,3.0);                   }
   double  getFactorDensityWToLb()   { return 1.0/this->getFactorDensityLbToW();                                       }

   double  getFactorPressureLbToW()  { return this->factorMassLbToW/(std::pow(factorTimeLbToW,2.0)*factorLengthLbToW); }
   double  getFactorPressureWToLb()  { return 1.0/this->getFactorPressureLbToW();                                      }

   double  getFactorMassLbToW()      { return this->factorMassLbToW;                                                   }
   double  getFactorMassWToLb()      { return 1.0/this->getFactorMassLbToW();                                          }

   double  getFactorForceLbToW()     { return factorMassLbToW*factorLengthLbToW/(factorTimeLbToW*factorTimeLbToW);     }
   double  getFactorForceWToLb()     { return 1.0/this->getFactorForceLbToW();                                         }

   double  getFactorAccLbToW()       { return factorLengthLbToW/(factorTimeLbToW*factorTimeLbToW);                     }
   double  getFactorAccWToLb()       { return 1.0/this->getFactorAccLbToW();                                           }

   double  getFactorTimeLbToW(double deltaX)        const { return factorTimeWithoutDx * deltaX;             }

   /*==========================================================*/
   friend inline std::ostream& operator << (std::ostream& os, LBMUnitConverter c) 
   {
      os<<c.toString();
      return os;
   }
   /*==========================================================*/
   std::string toString() 
   {
      std::ostringstream out;
      out<<"LB --> WORLD" << std::endl;
      out<<" * lentgh 1[dx  ] = " << std::setw(12) << this->getFactorLentghLbToW()    << " [m   ] " << std::endl;
      out<<" * time   1[dt  ] = " << std::setw(12) << this->getFactorTimeLbToW()      << " [s   ] " << std::endl;
      out<<" * mass   1[mass] = " << std::setw(12) << this->getFactorMassLbToW()      << " [kg  ] " << std::endl;
      out<<std::endl;                                                       
      out<<"WORLD --> LB" << std::endl;                                     
      out<<" * lentgh 1[m   ] = " << std::setw(12) << this->getFactorLentghWToLb()    << " [dx  ] " << std::endl;
      out<<" * time   1[s   ] = " << std::setw(12) << this->getFactorTimeWToLb()      << " [dt  ] " << std::endl;
      out<<" * mass   1[kg  ] = " << std::setw(12) << this->getFactorMassWToLb()      << " [mass] " << std::endl;
      out<<std::endl;
      out<<"LB --> WORLD (combined units)" << std::endl;
      out<<" * velocity     1 [dx/dt    ] = " << std::setw(12) << this->getFactorVelocityLbToW()  << " [m/s      ]" << std::endl;
      out<<" * density      1 [mass/dx^3] = " << std::setw(12) << this->getFactorDensityLbToW()   << " [kg/m^3   ]" << std::endl;
      out<<" * pressure     1 [F_lb/dx^2] = " << std::setw(12) << this->getFactorPressureLbToW()  << " [N/m^2    ]" << std::endl;
      out<<" * viscosity    1 [dx^2/dt  ] = " << std::setw(12) << this->getFactorViscosityLbToW() << " [m^2/s    ]" << std::endl;
      out<<" * force        1 [F_lb     ] = " << std::setw(12) << this->getFactorForceLbToW()     << " [N        ]" << std::endl;
      out<<" * acceleration 1 [dx/dt^2  ] = " << std::setw(12) << this->getFactorAccLbToW()       << " [m/s^2    ]" << std::endl;
      out<<std::endl;                                                                       
      out<<"WORLD --> LB (combined units)" << std::endl;                                    
      out<<" * velocity     1 [m/s      ] = " << std::setw(12) << this->getFactorVelocityWToLb()  << " [dx/dt    ]" << std::endl;
      out<<" * density      1 [kg/m^3   ] = " << std::setw(12) << this->getFactorDensityWToLb()   << " [mass/dx^3]" << std::endl;
      out<<" * pressure     1 [N/m^2    ] = " << std::setw(12) << this->getFactorPressureWToLb()  << " [F_lb/dx^2]" << std::endl;
      out<<" * viscosity    1 [m^2/s    ] = " << std::setw(12) << this->getFactorViscosityWToLb() << " [dx^2/dt  ]" << std::endl;
      out<<" * force        1 [N        ] = " << std::setw(12) << this->getFactorForceWToLb()     << " [F_lb     ]" << std::endl;
      out<<" * acceleration 1 [m/s^2    ] = " << std::setw(12) << this->getFactorAccWToLb()       << " [dx/dt^2  ]" << std::endl;

      return out.str();
   }

   void init(  const double& refLengthWorld, const double&  /*csWorld*/, const double& rhoWorld, const double& vWorld, 
               const double& refLengthLb, const double& rhoLb, const double& vLb  )
   {
      factorLengthLbToW = refLengthWorld / refLengthLb;
      factorTimeLbToW   = vLb / vWorld * factorLengthLbToW;
      factorMassLbToW   = rhoWorld/rhoLb*factorLengthLbToW*factorLengthLbToW*factorLengthLbToW;
      factorTimeWithoutDx=vLb/vWorld;
      this->refRhoLb = rhoLb;
   }

   protected:
   double factorLengthLbToW{1.0};
   double factorTimeLbToW{1.0};
   double factorMassLbToW{1.0};
   double refRhoLb{1.0};
   double factorTimeWithoutDx;
};

#endif //LBMUNITCONVERTER_H
