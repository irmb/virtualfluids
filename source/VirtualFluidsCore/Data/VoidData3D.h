#ifndef VoidData3D_H
#define VoidData3D_H

#include "EsoTwist3D.h"
#include <boost/serialization/serialization.hpp>

class VoidData3D;
typedef std::shared_ptr<VoidData3D> VoidData3DPtr;

class VoidData3D : public EsoTwist3D
{
public:
   VoidData3D() {};
   VoidData3D (size_t nx1, size_t nx2, size_t nx3, LBMReal value) 
   {
      this->NX1 = nx1;
      this->NX2 = nx2;
      this->NX3 = nx3;
   }
    ~VoidData3D() {};
    size_t getNX1() const { return NX1;}
    size_t getNX2() const { return NX2;}
    size_t getNX3() const { return NX3;}
    void getDistribution(LBMReal* const f, size_t x1, size_t x2, size_t x3) {}
    void setDistribution(const LBMReal* const f, size_t x1, size_t x2, size_t x3){}
    void getDistributionInv(LBMReal* const f, size_t x1, size_t x2, size_t x3) {}
    void setDistributionInv(const LBMReal* const f, size_t x1, size_t x2, size_t x3) {}
    void setDistributionForDirection(const LBMReal* const f, size_t x1, size_t x2, size_t x3, unsigned long int direction) {}
    void setDistributionForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, int direction) {}
    LBMReal getDistributionInvForDirection(size_t x1, size_t x2, size_t x3, int direction) {return 0.0;}
    void setDistributionInvForDirection(const LBMReal* const f, size_t x1, size_t x2, size_t x3, unsigned long int direction) {}
    void setDistributionInvForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, unsigned long int direction) {}
    LBMReal getDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction) {return 0.0;}
    void swap() {}
protected:
private:
   size_t NX1, NX2, NX3;
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object< EsoTwist3D >(*this);
   }
};

#endif
