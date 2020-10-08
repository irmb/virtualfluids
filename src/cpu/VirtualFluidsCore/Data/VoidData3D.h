#ifndef VoidData3D_H
#define VoidData3D_H

#include "EsoTwist3D.h"


class VoidData3D : public EsoTwist3D
{
public:
   VoidData3D() = default;;
   VoidData3D (size_t nx1, size_t nx2, size_t nx3, LBMReal  /*value*/) 
   {
      this->NX1 = nx1;
      this->NX2 = nx2;
      this->NX3 = nx3;
   }
    ~VoidData3D() override = default;;
    size_t getNX1() const override { return NX1;}
    size_t getNX2() const override { return NX2;}
    size_t getNX3() const override { return NX3;}
    void getDistribution(LBMReal* const f, size_t x1, size_t x2, size_t x3) override {}
    void setDistribution(const LBMReal* const f, size_t x1, size_t x2, size_t x3) override{}
    void getDistributionInv(LBMReal* const f, size_t x1, size_t x2, size_t x3) override {}
    void setDistributionInv(const LBMReal* const f, size_t x1, size_t x2, size_t x3) override {}
    void setDistributionForDirection(const LBMReal* const f, size_t x1, size_t x2, size_t x3, unsigned long int direction) override {}
    void setDistributionForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, int direction) override {}
    LBMReal getDistributionInvForDirection(size_t  /*x1*/, size_t  /*x2*/, size_t  /*x3*/, int  /*direction*/) override {return 0.0;}
    void setDistributionInvForDirection(const LBMReal* const f, size_t x1, size_t x2, size_t x3, unsigned long int direction) override {}
    void setDistributionInvForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, unsigned long int direction) override {}
    LBMReal getDistributionForDirection(size_t  /*x1*/, size_t  /*x2*/, size_t  /*x3*/, int  /*direction*/) override {return 0.0;}
    void swap() override {}
protected:
private:
   size_t NX1, NX2, NX3;
};

#endif
