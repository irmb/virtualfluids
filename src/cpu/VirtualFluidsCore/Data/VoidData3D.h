#ifndef VoidData3D_H
#define VoidData3D_H

#include "EsoTwist3D.h"

class VoidData3D : public EsoTwist3D
{
public:
    VoidData3D() = default;
    
    VoidData3D(size_t nx1, size_t nx2, size_t nx3, real /*value*/)
    {
        this->NX1 = nx1;
        this->NX2 = nx2;
        this->NX3 = nx3;
    }
    ~VoidData3D() override = default;
    
    size_t getNX1() const override { return NX1; }
    size_t getNX2() const override { return NX2; }
    size_t getNX3() const override { return NX3; }
    void getPreCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3) override {}
    void setPostCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3) override {}
    void getPostCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3) override {}
    void setPreCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3) override {}
    void setPostCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                     unsigned long int direction) override
    {
    }
    void setPostCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3, int direction) override {}
    real getDistributionInvForDirection(size_t /*x1*/, size_t /*x2*/, size_t /*x3*/, int /*direction*/) override
    {
        return 0.0;
    }
    void setPreCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                        unsigned long int direction) override
    {
    }
    void setPreCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3,
                                        unsigned long int direction) override
    {
    }
    real getPreCollisionDistributionForDirection(size_t /*x1*/, size_t /*x2*/, size_t /*x3*/, int /*direction*/) override
    {
        return 0.0;
    }
    void swap() override {}

protected:
private:
    size_t NX1, NX2, NX3;
};

#endif
