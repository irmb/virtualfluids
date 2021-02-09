#ifndef CumulantLBMKernel_h__
#define CumulantLBMKernel_h__

#include "BasicLBMKernel.h"
#include "BCProcessor.h"
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
   virtual ~CumulantLBMKernel(void);
   //virtual void calculate(int step);
   SPtr<LBMKernel> clone();
   double getCalculationTime();
   void setBulkOmegaToOmega(bool value);
   void setRelaxationParameter(Parameter p);
protected:
   void initData() override;
   void nodeCollision(int step, int x1, int x2, int x3) override;
   void initDataSet();
   LBMReal f[D3Q27System::ENDF + 1];

   UbTimer timer;

   LBMReal OxyyMxzz;
   Parameter parameter;

   CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
   CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
   CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributions;

   mu::value_type muX1, muX2, muX3;
   mu::value_type muDeltaT;
   mu::value_type muNu;
   LBMReal forcingX1;
   LBMReal forcingX2;
   LBMReal forcingX3;

   // bulk viscosity
   bool bulkOmegaToOmega;
   LBMReal OxxPyyPzz;

   LBMReal omega;
};
#endif // CumulantLBMKernel_h__