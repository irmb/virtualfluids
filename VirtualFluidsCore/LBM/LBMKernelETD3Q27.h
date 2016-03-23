#ifndef LBMKERNELETD3Q27_H
#define LBMKERNELETD3Q27_H

#include "LBMKernel3D.h"
#include "BCProcessor.h"
#include "D3Q27ETBCProcessor.h"
#include "EsoTwist3D.h"

#include <boost/serialization/base_object.hpp>

class LBMKernelETD3Q27;
typedef boost::shared_ptr<LBMKernelETD3Q27> LBMKernelETD3Q27Ptr;

class LBMKernelETD3Q27 : public LBMKernel3D
{
public:
   LBMKernelETD3Q27();
   LBMKernelETD3Q27(int nx1, int nx2, int nx3); 
   virtual ~LBMKernelETD3Q27(void);
   void swapDistributions();
protected:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<LBMKernel3D>(*this);
      ar & nx1 & nx2 & nx3;
   }
   int nx1, nx2, nx3;
};
#endif
