#ifndef D3Q27EsoTwist3DSplittedVector_h
#define D3Q27EsoTwist3DSplittedVector_h

#include "EsoTwist3D.h"
#include "D3Q27System.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>

class D3Q27EsoTwist3DSplittedVector : public EsoTwist3D
{
public:
   D3Q27EsoTwist3DSplittedVector();
   D3Q27EsoTwist3DSplittedVector(size_t nx1, size_t nx2, size_t nx3, LBMReal value);
   //////////////////////////////////////////////////////////////////////////
   ~D3Q27EsoTwist3DSplittedVector();
   //////////////////////////////////////////////////////////////////////////
   void swap();
   //////////////////////////////////////////////////////////////////////////
   virtual void getDistribution( LBMReal* const f, size_t x1, size_t x2, size_t x3);
   //////////////////////////////////////////////////////////////////////////
   virtual void setDistribution(const LBMReal* const f, size_t x1, size_t x2, size_t x3);
   ////////////////////////////////////////////////////////////////////////
   virtual void getDistributionInv( LBMReal* const f, size_t x1, size_t x2, size_t x3);
   //////////////////////////////////////////////////////////////////////////
   virtual void setDistributionInv(const LBMReal* const f, size_t x1, size_t x2, size_t x3);
   //////////////////////////////////////////////////////////////////////////
   virtual void setDistributionForDirection(const LBMReal* const f, size_t x1, size_t x2, size_t x3, unsigned long int direction);
   //////////////////////////////////////////////////////////////////////////
   virtual void setDistributionForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, int direction);
   //////////////////////////////////////////////////////////////////////////
   virtual LBMReal getDistributionInvForDirection(size_t x1, size_t x2, size_t x3, int direction);
   //////////////////////////////////////////////////////////////////////////
   virtual void setDistributionInvForDirection(const LBMReal* const f, size_t x1, size_t x2, size_t x3, unsigned long int direction);
   //////////////////////////////////////////////////////////////////////////
   virtual void setDistributionInvForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, unsigned long int direction);
   //////////////////////////////////////////////////////////////////////////
   virtual LBMReal getDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction);
   //////////////////////////////////////////////////////////////////////////
   size_t getNX1() const;
   //////////////////////////////////////////////////////////////////////////
   size_t getNX2() const;
   //////////////////////////////////////////////////////////////////////////
   size_t getNX3() const;
   //////////////////////////////////////////////////////////////////////////
   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr getLocalDistributions();
   //////////////////////////////////////////////////////////////////////////
   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr getNonLocalDistributions();
   //////////////////////////////////////////////////////////////////////////
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr getZeroDistributions();
   //////////////////////////////////////////////////////////////////////////
   void setNX1(size_t newNX1);
   void setNX2(size_t newNX2);
   void setNX3(size_t newNX3);
   void setLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr array);
   void setNonLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr array);
   void setZeroDistributions(CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr array);
   
protected:
   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr   zeroDistributions;
   size_t NX1, NX2, NX3;

   friend class MPIIORestartCoProcessor;

   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object< EsoTwist3D >(*this);
      ar & NX1; 
      ar & NX2; 
      ar & NX3;
      ar &  localDistributions;
      ar &  nonLocalDistributions;
      ar &  zeroDistributions;
   }
};

#endif
