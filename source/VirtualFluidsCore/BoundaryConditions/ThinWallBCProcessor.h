#ifndef ThinWallBCProcessor_H
#define ThinWallBCProcessor_H

#include "BCProcessor.h"
#include "BoundaryConditions.h"
#include "EsoTwist3D.h"
#include "BCArray3D.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

#include <boost/serialization/base_object.hpp>

class ThinWallBCProcessor;
typedef boost::shared_ptr<ThinWallBCProcessor> ThinWallBCProcessorPtr;

class ThinWallBCProcessor : public BCProcessor
{
public:
   ThinWallBCProcessor();
   ThinWallBCProcessor(LBMKernelPtr kernel);
   ~ThinWallBCProcessor();
   BCProcessorPtr clone(LBMKernelPtr kernel);
   void applyPostCollisionBC();
protected:
private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<BCProcessor>(*this);
   }
};

#endif
