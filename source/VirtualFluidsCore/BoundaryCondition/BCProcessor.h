#ifndef BCPROCESSOR_H
#define BCPROCESSOR_H

#include "EsoTwist3D.h"


#include <boost/serialization/serialization.hpp>
#include <boost/shared_ptr.hpp>

class BCProcessor;
typedef boost::shared_ptr<BCProcessor> BCProcessorPtr;

#include "LBMKernel3D.h"
#include "BoundaryCondition.h"

class BCProcessor
{
public:
   BCProcessor();
   virtual ~BCProcessor();
   virtual void applyPreCollisionBC() = 0;
   virtual void applyPostCollisionBC() = 0;
   virtual void addBC(BoundaryConditionPtr bc) = 0;
   virtual BoundaryConditionPtr getBC(BoundaryCondition::Type type) = 0;
   virtual BCProcessorPtr clone(LBMKernel3DPtr kernel) = 0;
protected:
private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
   }
};

#endif
