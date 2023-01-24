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
//! \file RheologyK17LBMKernel.h
//! \ingroup LBM
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef RheologyK17LBMKernel_h__
#define RheologyK17LBMKernel_h__

#include "LBMKernel.h"
#include "BCProcessor.h"
#include "D3Q27System.h"
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

//! \brief   compressible cumulant LBM kernel with rheological properties of shear and bulk viscosity for non-Newtonian fluids.
//! \details CFD solver that use Cascaded Cumulant Lattice Boltzmann method for D3Q27 model
//! \author  K. Kutscher, M. Geier
class RheologyK17LBMKernel :  public LBMKernel
{
public:
   //! This option set relaxation parameter: NORMAL  
   enum Parameter{NORMAL, MAGIC};
public:
   RheologyK17LBMKernel();
   virtual ~RheologyK17LBMKernel(void);
   virtual void calculate(int step) override;
   virtual SPtr<LBMKernel> clone() override;
   double getCalculationTime() override;
   //! The value should not be equal to a shear viscosity
   void setBulkViscosity(real value);
protected:
   virtual void initDataSet();

   virtual real getRheologyCollFactor(real omegaInf, real shearRate, real drho) const
   {
       UB_THROW(UbException("real getRheologyCollFactor() - belongs in the derived class"));
   }

   real f[D3Q27System::ENDF+1];

   UbTimer timer;

   CbArray4D<real,IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
   CbArray4D<real,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
   CbArray3D<real,IndexerX3X2X1>::CbArray3DPtr   zeroDistributions;

   mu::value_type muX1,muX2,muX3;
   mu::value_type muDeltaT;
   mu::value_type muNu;
   real forcingX1;
   real forcingX2;
   real forcingX3;
   
   // bulk viscosity
   real OxxPyyPzz; //omega2 (bulk viscosity)
   real bulkViscosity;

};
#endif // RheologyK17LBMKernel_h__


