#ifndef BCPROCESSSOR_H
#define BCPROCESSSOR_H

#include "BCAlgorithm.h"
#include "Data/EsoTwist3D.h"
#include "BCArray3D.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

#include <vector>

#include <boost/serialization/base_object.hpp>

class BCProcessor;
typedef boost::shared_ptr<BCProcessor> BCProcessorPtr;

#include "LBM/LBMKernel.h"

class BCProcessor
{
public:
   BCProcessor();
   BCProcessor(LBMKernelPtr kernel);
   virtual ~BCProcessor();
   virtual BCArray3D& getBCArray();
   virtual BCProcessorPtr clone(LBMKernelPtr kernel);

   void addBC(BCAlgorithmPtr bc);
   void applyPreCollisionBC();
   void applyPostCollisionBC();
   void clearBC();
protected:
   std::vector<BCAlgorithmPtr> preBC;
   std::vector<BCAlgorithmPtr> postBC;
   BCArray3D bcArray;

private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & bcArray;
      //ar & preBC;
      //ar & postBC;
   }
};

#endif
