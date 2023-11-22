//Cascaded Cumulant LBM

#ifndef K16IncompressibleNavierStokes_H
#define K16IncompressibleNavierStokes_H

#include <basics/Timer/Timer.h>

#include "LBMKernel.h"
#include "BCSet.h"
#include "D3Q27System.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

//! \brief   Cascaded Cumulant LBM kernel. 
//! \details CFD solver that use Cascaded Cumulant Lattice Boltzmann method for D3Q27 model
//! \author  K. Kutscher, M. Geier
class K16IncompressibleNavierStokes : public LBMKernel
{
public:
   //! This option set relaxation parameter: NORMAL  
   enum Parameter{NORMAL, MAGIC};
public:
   K16IncompressibleNavierStokes();
   ~K16IncompressibleNavierStokes() override;
   void calculate(int step) override;
   SPtr<LBMKernel> clone() override;
   real getCalculationTime() override;
   void setRelaxationParameter(Parameter p);
protected:
   virtual void initDataSet();
   real f[D3Q27System::ENDF+1];

   Timer timer;

   real OxyyMxzz;
   Parameter parameter;

   CbArray4D<real,IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
   CbArray4D<real,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
   CbArray3D<real,IndexerX3X2X1>::CbArray3DPtr   zeroDistributions;

   mu::value_type muX1,muX2,muX3;
   mu::value_type muDeltaT;
   mu::value_type muNu;
   real forcingX1;
   real forcingX2;
   real forcingX3;
};

#endif
