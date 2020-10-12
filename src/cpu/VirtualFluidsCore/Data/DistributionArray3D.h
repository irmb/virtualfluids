#ifndef DistributionArray3D_H
#define DistributionArray3D_H

#include <LBMSystem.h>

class DistributionArray3D
{
public:
   DistributionArray3D() = default;;
   virtual ~DistributionArray3D()= default;;

   virtual size_t getNX1() const = 0;
   virtual size_t getNX2() const = 0;
   virtual size_t getNX3() const = 0;
   virtual void getDistribution(LBMReal* const f, size_t x1, size_t x2, size_t x3) = 0;
   virtual void setDistribution(const LBMReal* const f, size_t x1, size_t x2, size_t x3) = 0;
   virtual void getDistributionInv( LBMReal* const f, size_t x1, size_t x2, size_t x3) = 0;
   virtual void setDistributionInv(const LBMReal* const f, size_t x1, size_t x2, size_t x3) = 0;
   virtual void setDistributionForDirection(const LBMReal* const f, size_t x1, size_t x2, size_t x3, unsigned long int direction) = 0;
   virtual void setDistributionForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, int direction) = 0;
   virtual LBMReal getDistributionInvForDirection(size_t x1, size_t x2, size_t x3, int direction) = 0;
   virtual void setDistributionInvForDirection(const LBMReal* const f, size_t x1, size_t x2, size_t x3, unsigned long int direction) = 0;
   virtual void setDistributionInvForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, unsigned long int direction) = 0;
   virtual LBMReal getDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction) = 0;
   virtual void swap() = 0;
protected:
private:

};

#endif
