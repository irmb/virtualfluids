#ifndef CompressibleCumulantLBMKernel_h__
#define CompressibleCumulantLBMKernel_h__

#include "LBMKernel.h"
#include "BCProcessor.h"
#include "D3Q27System.h"
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

//! \brief   compressible cumulant LBM kernel. 
//! \details CFD solver that use Cascaded Cumulant Lattice Boltzmann method for D3Q27 model
//! \author  K. Kutscher, M. Geier
class CompressibleCumulantLBMKernel :  public LBMKernel
{
public:
   //! This option set relaxation parameter: NORMAL  
   enum Parameter{NORMAL, MAGIC};
public:
   CompressibleCumulantLBMKernel();
   ~CompressibleCumulantLBMKernel(void) override;
   void calculate(int step) override;
   SPtr<LBMKernel> clone() override;
   double getCalculationTime() override;
   void setBulkOmegaToOmega(bool value);
   void setRelaxationParameter(Parameter p);
protected:
   virtual void initDataSet();
   LBMReal f[D3Q27System::ENDF+1];

   UbTimer timer;

   LBMReal OxyyMxzz;
   Parameter parameter;

   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr   zeroDistributions;

   mu::value_type muX1,muX2,muX3;
   mu::value_type muDeltaT;
   mu::value_type muNu;
   LBMReal forcingX1;
   LBMReal forcingX2;
   LBMReal forcingX3;
   
   // bulk viscosity
   bool bulkOmegaToOmega;
   LBMReal OxxPyyPzz; 
};
#endif // CompressibleCumulantLBMKernel_h__


