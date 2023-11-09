#ifndef CumulantLBMKernel_h__
#define CumulantLBMKernel_h__

#include "BasicLBMKernel.h"
#include "BCSet.h"
#include "D3Q27System.h"
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

//! \brief   compressible cumulant LBM kernel. 
//! \details CFD solver that use Cascaded Cumulant Lattice Boltzmann method for D3Q27 model
//! \author  K. Kutscher, M. Geier
class CumulantLBMKernel : public BasicLBMKernel
{
public:
   //! This option set relaxation parameter: NORMAL  
   enum Parameter { NORMAL, MAGIC };
public:
   CumulantLBMKernel();
   virtual ~CumulantLBMKernel() = default;
   //virtual void calculate(int step);
   SPtr<LBMKernel> clone() override;
   real getCalculationTime() override;
   void setBulkOmegaToOmega(bool value);
   void setRelaxationParameter(Parameter p);
protected:
   void initData() override;
   void nodeCollision(int step, int x1, int x2, int x3) override;
   void initDataSet();
   real f[D3Q27System::ENDF + 1];

   UbTimer timer;

   real OxyyMxzz;
   Parameter parameter;

   CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
   CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
   CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr   zeroDistributions;

   mu::value_type muX1, muX2, muX3;
   mu::value_type muDeltaT;
   mu::value_type muNu;
   real forcingX1;
   real forcingX2;
   real forcingX3;

   // bulk viscosity
   bool bulkOmegaToOmega;
   real OxxPyyPzz;

   real omega;
};
#endif // CumulantLBMKernel_h__