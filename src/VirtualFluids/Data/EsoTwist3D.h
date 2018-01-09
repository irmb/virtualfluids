#ifndef ESOTWIST3D_H
#define ESOTWIST3D_H

#include "DistributionArray3D.h"
#include <LBMSystem.h>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>

class EsoTwist3D;
typedef std::shared_ptr<EsoTwist3D> EsoTwist3DPtr;

class EsoTwistD3Q27UnrollArray{};
class EsoTwistPlusD3Q27UnrollArray{};
class EsoTwistPlusD3Q19UnrollArray{};

class EsoTwist3D : public DistributionArray3D
{
public:
   EsoTwist3D(){};
   virtual ~EsoTwist3D(){};
   //////////////////////////////////////////////////////////////////////////
   virtual void swap() = 0;
   //////////////////////////////////////////////////////////////////////////
   virtual void getDistribution(LBMReal* const f, size_t x1, size_t x2, size_t x3) = 0;
   //////////////////////////////////////////////////////////////////////////
   virtual void setDistribution(const LBMReal* const f, size_t x1, size_t x2, size_t x3) = 0;
   ////////////////////////////////////////////////////////////////////////
   virtual void getDistributionInv( LBMReal* const f, size_t x1, size_t x2, size_t x3) = 0;
   //////////////////////////////////////////////////////////////////////////
   virtual void setDistributionInv(const LBMReal* const f, size_t x1, size_t x2, size_t x3) = 0;
   //////////////////////////////////////////////////////////////////////////
   virtual void setDistributionForDirection(const LBMReal* const f, size_t x1, size_t x2, size_t x3, unsigned long int direction) = 0;
   //////////////////////////////////////////////////////////////////////////
   virtual void setDistributionForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, int direction) = 0;
   //////////////////////////////////////////////////////////////////////////
   //virtual void getDistributionInvForDirection(LBMReal* const& f, const size_t& x1, const size_t& x2, const size_t& x3, const unsigned long int& direction) = 0;
   //////////////////////////////////////////////////////////////////////////
   virtual LBMReal getDistributionInvForDirection(size_t x1, size_t x2, size_t x3, int direction) = 0;
   //////////////////////////////////////////////////////////////////////////
   virtual void setDistributionInvForDirection(const LBMReal* const f, size_t x1, size_t x2, size_t x3, unsigned long int direction) = 0;
   //////////////////////////////////////////////////////////////////////////
   virtual void setDistributionInvForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, unsigned long int direction) = 0;
   //////////////////////////////////////////////////////////////////////////
   virtual LBMReal getDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction) = 0;
   //////////////////////////////////////////////////////////////////////////
   virtual size_t getNX1() const = 0;
   //////////////////////////////////////////////////////////////////////////
   virtual size_t getNX2() const = 0;
   //////////////////////////////////////////////////////////////////////////
   virtual size_t getNX3() const = 0;
   //////////////////////////////////////////////////////////////////////////
  
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object< DistributionArray3D >(*this);
   }
};

#endif
