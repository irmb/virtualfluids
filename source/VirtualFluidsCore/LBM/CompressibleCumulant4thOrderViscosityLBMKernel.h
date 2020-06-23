#ifndef CompressibleCumulant4thOrderViscosityLBMKernel_h__
#define CompressibleCumulant4thOrderViscosityLBMKernel_h__

#include "LBMKernel.h"
#include "BCProcessor.h"
#include "D3Q27System.h"
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

//! \brief   compressible cumulant LBM kernel. 
//! \details CFD solver that use Cascaded Cumulant Lattice Boltzmann method for D3Q27 model
//! \author  K. Kutscher, M. Geier
class CompressibleCumulant4thOrderViscosityLBMKernel :  public LBMKernel
{
public:
   //! This option set relaxation parameter: NORMAL  
   enum Parameter{NORMAL, MAGIC};
public:
   CompressibleCumulant4thOrderViscosityLBMKernel();
   virtual ~CompressibleCumulant4thOrderViscosityLBMKernel(void);
   virtual void calculate(int step);
   virtual SPtr<LBMKernel> clone();
   double getCalculationTime() override;
   //! The value should not be equal to a shear viscosity
   void setBulkViscosity(LBMReal value);
protected:
   virtual void initDataSet();
   LBMReal f[D3Q27System::ENDF+1];

   UbTimer timer;

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
   LBMReal OxxPyyPzz; //omega2 (bulk viscosity)
   LBMReal bulkViscosity;

};
#endif // CompressibleCumulant4thOrderViscosityLBMKernel_h__


