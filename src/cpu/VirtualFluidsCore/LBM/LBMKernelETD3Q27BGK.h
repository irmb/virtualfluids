#ifndef LBMKERNELETD3Q27BGK_H
#define LBMKERNELETD3Q27BGK_H

#include "LBMKernel.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"



class LBMKernelETD3Q27BGK :  public LBMKernel
{
public:
   LBMKernelETD3Q27BGK();
   ~LBMKernelETD3Q27BGK() override;
   void calculate(int step)override;
   SPtr<LBMKernel> clone()override;
   real getCalculationTime() override;

private:
   void initDataSet();
   //void collideAllCompressible();
   //void collideAllIncompressible();

   CbArray4D<real,IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
   CbArray4D<real,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
   CbArray3D<real,IndexerX3X2X1>::CbArray3DPtr   zeroDistributions;

   mu::value_type muX1,muX2,muX3;
   real forcingX1;
   real forcingX2;
   real forcingX3;


};

#endif
