#ifndef D3Q27EsoTwist3DSoA_h
#define D3Q27EsoTwist3DSoA_h

#include "EsoTwist3D.h"
//#include "D3Q27System.h"
//#include "basics/container/CbArray4D.h"
#include <basics/container/CbArray3D.h>
//#include <boost/serialization/serialization.hpp>
//#include <boost/serialization/base_object.hpp>


struct Distributions
{
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr E;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr W;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr N;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr S;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr T;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr B;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr NE;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr SW;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr SE;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr NW;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr TE;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr BW;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr BE;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr TW;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr TN;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr BS;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr BN;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr TS;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr TNE;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr TNW;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr TSE;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr TSW;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr BNE;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr BNW;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr BSE;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr BSW;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr ZERO;
};


class D3Q27EsoTwist3DSoA : public EsoTwist3D
{
public:
   D3Q27EsoTwist3DSoA();
   D3Q27EsoTwist3DSoA(const size_t& nx1, const size_t& nx2, const size_t& nx3, LBMReal value);
   //////////////////////////////////////////////////////////////////////////
   ~D3Q27EsoTwist3DSoA() override;
   //////////////////////////////////////////////////////////////////////////
   void swap() override;
   //////////////////////////////////////////////////////////////////////////
   void getDistribution( LBMReal* const f, size_t x1, size_t x2, size_t x3) override;
   //////////////////////////////////////////////////////////////////////////
   void setDistribution(const LBMReal* const f, size_t x1, size_t x2, size_t x3) override;
   ////////////////////////////////////////////////////////////////////////
   void getDistributionInv( LBMReal* const f, size_t x1, size_t x2, size_t x3) override;
   //////////////////////////////////////////////////////////////////////////
   void setDistributionInv(const LBMReal* const f, size_t x1, size_t x2, size_t x3) override;
   //////////////////////////////////////////////////////////////////////////
   void setDistributionForDirection(const LBMReal* const f, size_t x1, size_t x2, size_t x3, unsigned long int direction) override;
   //////////////////////////////////////////////////////////////////////////
   void setDistributionForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, int direction) override;
   //////////////////////////////////////////////////////////////////////////
   LBMReal getDistributionInvForDirection(size_t x1, size_t x2, size_t x3, int direction) override;
   //////////////////////////////////////////////////////////////////////////
   void setDistributionInvForDirection(const LBMReal* const f, size_t x1, size_t x2, size_t x3, unsigned long int direction) override;
   //////////////////////////////////////////////////////////////////////////
   void setDistributionInvForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, unsigned long int direction) override;
   //////////////////////////////////////////////////////////////////////////
   LBMReal getDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction) override;
   //////////////////////////////////////////////////////////////////////////
   size_t getNX1() const override;
   //////////////////////////////////////////////////////////////////////////
   size_t getNX2() const override;
   //////////////////////////////////////////////////////////////////////////
   size_t getNX3() const override;
   //////////////////////////////////////////////////////////////////////////
   Distributions getDistributions();
   //////////////////////////////////////////////////////////////////////////
   void getDistributionAfterLastStep(LBMReal* const f, size_t x1, size_t x2, size_t x3);

protected:
   Distributions d;
   size_t NX1, NX2, NX3;

};

#endif

