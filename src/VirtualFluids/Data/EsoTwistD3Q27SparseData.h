#ifndef EsoTwistD3Q27SparseData_h
#define EsoTwistD3Q27SparseData_h

#include "EsoTwist3D.h"
#include <D3Q27System.h>

#include <boost/unordered_map.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>

class EsoTwistD3Q27SparseData;
typedef std::shared_ptr<EsoTwistD3Q27SparseData> EsoTwistD3Q27SparseDataPtr;

class EsoTwistD3Q27SparseData : public EsoTwist3D
{
public:
   typedef boost::unordered_map<size_t, LBMReal> SparseData;
   typedef std::shared_ptr<SparseData> SparseDataPtr;
   //typedef std::map<size_t, LBMReal> SparseData;
   //typedef std::shared_ptr<SparseData> SparseDataPtr;
public:
   EsoTwistD3Q27SparseData();
   EsoTwistD3Q27SparseData(size_t ibx[3], size_t nx[3], size_t level, double value );
   //////////////////////////////////////////////////////////////////////////
   ~EsoTwistD3Q27SparseData();
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
   ////////////////////////////////////////////////////////////////////////
   SparseDataPtr getLocalDistributions();
   //////////////////////////////////////////////////////////////////////////
   SparseDataPtr getNonLocalDistributions();
   //////////////////////////////////////////////////////////////////////////
   SparseDataPtr getZeroDistributions();
   //////////////////////////////////////////////////////////////////////////
   static void setSize(int nx[]);
   //////////////////////////////////////////////////////////////////////////
   inline size_t index4D(size_t x1, size_t x2, size_t x3, size_t x4)
   {
      //return nx1*(nx2*(nx3*x4+x3)+x2)+x1;
      return nx1*(nx2*(nx3*(nx4*level+(x4+ox3))+(x3+ox2))+(x2+ox1))+x1;
   }
   //////////////////////////////////////////////////////////////////////////
   inline size_t index3D(size_t x1, size_t x2, size_t x3)
   {
      //return  nx2 * ( nx3 * x3 + x2) + x1 ;
      return nx2*(nx3*(nx4*level+(x3+ox3))+(x2+ox2))+(x1+ox1);
   }

protected:
   static SparseData localDistributions;
   static SparseData nonLocalDistributions;
   static SparseData zeroDistributions;
   SparseDataPtr ld;  //local distributions;
   SparseDataPtr nld; //non local distributions;
   SparseDataPtr zd;  //zero distributions;
   size_t NX1, NX2, NX3;
   static size_t nx1, nx2, nx3, nx4, nx5;
   size_t ox1, ox2, ox3;
   size_t level;




   //friend class boost::serialization::access;
   //template<class Archive>
   //void serialize(Archive & ar, const unsigned int version)
   //{
   //   ar & boost::serialization::base_object< EsoTwist3D >(*this);
   //   ar & NX1; 
   //   ar & NX2; 
   //   ar & NX3;
   //   ar &  localDistributions;
   //   ar &  nonLocalDistributions;
   //   ar &  zeroDistributions;
   //}
};

#endif
