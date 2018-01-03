#ifndef ThinWallBCProcessor_H
#define ThinWallBCProcessor_H

#include <memory>

#include "BCProcessor.h"

#include <boost/serialization/base_object.hpp>

class ThinWallBCProcessor;
typedef std::shared_ptr<ThinWallBCProcessor> ThinWallBCProcessorPtr;

class LBMKernel;

class ThinWallBCProcessor : public BCProcessor
{
public:
   ThinWallBCProcessor();
   ThinWallBCProcessor(std::shared_ptr<LBMKernel> kernel);
   ~ThinWallBCProcessor();
   std::shared_ptr<BCProcessor> clone(std::shared_ptr<LBMKernel> kernel);
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
