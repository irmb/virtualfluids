#ifndef D3Q27ETBCPROCESSSOR_H
#define D3Q27ETBCPROCESSSOR_H

#include "BCProcessor.h"
#include "D3Q27BoundaryCondition.h"
#include "EsoTwist3D.h"
#include "BCArray3D.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"
//#include "BoundaryCondition.h"
#include <vector>

#include <boost/serialization/base_object.hpp>

class D3Q27ETBCProcessor;
typedef boost::shared_ptr<D3Q27ETBCProcessor> D3Q27ETBCProcessorPtr;

class D3Q27ETBCProcessor : public BCProcessor
{
public:
   D3Q27ETBCProcessor();
   D3Q27ETBCProcessor(LBMKernel3DPtr kernel);
   virtual ~D3Q27ETBCProcessor();
   //virtual void applyBC();
   virtual BCArray3D<D3Q27BoundaryCondition>& getBCArray();
   virtual BCProcessorPtr clone(LBMKernel3DPtr kernel);

   void addBC(BoundaryConditionPtr bc);
   BoundaryConditionPtr getBC(BoundaryCondition::Type type);
   void applyPreCollisionBC();
   void applyPostCollisionBC();
   //void init();
protected:
   std::vector<BoundaryConditionPtr> preBC;
   std::vector<BoundaryConditionPtr> postBC;
   BCArray3D<D3Q27BoundaryCondition> bcArray;

private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<BCProcessor>(*this);
      ar & bcArray;
      ar & preBC;
      ar & postBC;
   }
};

#endif
