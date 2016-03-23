#ifndef LBMKERNELETD3Q27BGK_H
#define LBMKERNELETD3Q27BGK_H

#include "LBMKernelETD3Q27.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"
#include <boost/serialization/export.hpp>

class LBMKernelETD3Q27BGK;
typedef boost::shared_ptr<LBMKernelETD3Q27BGK> LBMKernelETD3Q27BGKPtr;

class LBMKernelETD3Q27BGK :  public LBMKernelETD3Q27
{
public:
   LBMKernelETD3Q27BGK();
   LBMKernelETD3Q27BGK(int nx1, int nx2, int nx3, bool compressible);
   ~LBMKernelETD3Q27BGK(void);
   void calculate();
   LBMKernel3DPtr clone();
   double getCallculationTime();

private:
   void init();
   void collideAllCompressible();
   void collideAllIncompressible();

   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr   zeroDistributions;

   mu::value_type muX1,muX2,muX3;
   LBMReal forcingX1;
   LBMReal forcingX2;
   LBMReal forcingX3;

   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<LBMKernelETD3Q27>(*this);
   }

};

#endif
