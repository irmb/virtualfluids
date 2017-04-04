#ifndef VoidLBMKernel_h__
#define VoidLBMKernel_h__

#include "LBMKernel.h"

class VoidLBMKernel : public LBMKernel
{
public:
   VoidLBMKernel();
   VoidLBMKernel(int nx1, int nx2, int nx3);
   ~VoidLBMKernel();
   LBMKernelPtr clone();
   void calculate();
   void swapDistributions();
   double getCallculationTime();
protected:
private:
   int nx1, nx2, nx3;
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<LBMKernel>(*this);
      ar & nx1 & nx2 & nx3;
   }

};
#endif // VoidLBMKernel_h__
