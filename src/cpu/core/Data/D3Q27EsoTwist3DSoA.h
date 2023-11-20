#ifndef D3Q27EsoTwist3DSoA_h
#define D3Q27EsoTwist3DSoA_h

#include "EsoTwist3D.h"
//#include "D3Q27System.h"
//#include "basics/container/CbArray4D.h"
#include <basics/container/CbArray3D.h>
//#include <boost/serialization/serialization.hpp>
//#include <boost/serialization/base_object.hpp>

struct Distributions {
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr E;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr W;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr N;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr S;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr T;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr B;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr NE;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr SW;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr SE;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr NW;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr TE;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr BW;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr BE;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr TW;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr TN;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr BS;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr BN;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr TS;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr TNE;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr TNW;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr TSE;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr TSW;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr BNE;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr BNW;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr BSE;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr BSW;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr REST;
};

class D3Q27EsoTwist3DSoA : public EsoTwist3D
{
public:
    D3Q27EsoTwist3DSoA();
    D3Q27EsoTwist3DSoA(const size_t &nx1, const size_t &nx2, const size_t &nx3, real value);
    //////////////////////////////////////////////////////////////////////////
    ~D3Q27EsoTwist3DSoA() override;
    //////////////////////////////////////////////////////////////////////////
    void swap() override;
    //////////////////////////////////////////////////////////////////////////
    void getPreCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3) override;
    //////////////////////////////////////////////////////////////////////////
    void setPostCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3) override;
    ////////////////////////////////////////////////////////////////////////
    void getPostCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3) override;
    //////////////////////////////////////////////////////////////////////////
    void setPreCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3) override;
    //////////////////////////////////////////////////////////////////////////
    void setPostCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                     unsigned long int direction) override;
    //////////////////////////////////////////////////////////////////////////
    void setPostCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3, int direction) override;
    //////////////////////////////////////////////////////////////////////////
    real getPostCollisionDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction) override;
    //////////////////////////////////////////////////////////////////////////
    void setPreCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                        unsigned long int direction) override;
    //////////////////////////////////////////////////////////////////////////
    void setPreCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3,
                                        unsigned long int direction) override;
    //////////////////////////////////////////////////////////////////////////
    real getPreCollisionDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction) override;
    //////////////////////////////////////////////////////////////////////////
    size_t getNX1() const override;
    //////////////////////////////////////////////////////////////////////////
    size_t getNX2() const override;
    //////////////////////////////////////////////////////////////////////////
    size_t getNX3() const override;
    //////////////////////////////////////////////////////////////////////////
    Distributions getDistributions();
    //////////////////////////////////////////////////////////////////////////
    void getDistributionAfterLastStep(real *const f, size_t x1, size_t x2, size_t x3);

protected:
    Distributions d;
    size_t NX1, NX2, NX3;
};

#endif
