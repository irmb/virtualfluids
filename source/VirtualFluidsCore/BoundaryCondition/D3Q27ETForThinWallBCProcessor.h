#ifndef D3Q27ETFORTHINWALLBCPROCESSSOR_H
#define D3Q27ETFORTHINWALLBCPROCESSSOR_H

#include "D3Q27ETBCProcessor.h"
#include "D3Q27BoundaryCondition.h"
#include "EsoTwist3D.h"
#include "BCArray3D.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

#include <boost/serialization/base_object.hpp>

class D3Q27ETForThinWallBCProcessor;
typedef boost::shared_ptr<D3Q27ETForThinWallBCProcessor> D3Q27ETForThinWallBCProcessorPtr;

class D3Q27ETForThinWallBCProcessor : public D3Q27ETBCProcessor
{
public:
   D3Q27ETForThinWallBCProcessor();
   D3Q27ETForThinWallBCProcessor(LBMKernel3DPtr kernel);
   virtual ~D3Q27ETForThinWallBCProcessor();
   virtual BCProcessorPtr clone(LBMKernel3DPtr kernel);
   void applyPostCollisionBC();
protected:
   //EsoTwist3DPtr distributionsTemp;
private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<D3Q27ETBCProcessor>(*this);
      //ar & distributionsTemp;
   }
};

#endif
