#ifndef BoundaryConditionProcessor_h__
#define BoundaryConditionProcessor_h__

#include "BoundaryCondition.h"
#include <vector>

class BoundaryConditionProcessor;
typedef boost::shared_ptr<BoundaryConditionProcessor> BoundaryConditionProcessorPtr;

class BoundaryConditionProcessor
{
public:
   BoundaryConditionProcessor();
   virtual ~BoundaryConditionProcessor() {}

   void addBC(BoundaryConditionPtr bc);
   BoundaryConditionPtr getBC(BoundaryCondition::Type type);
   void applyPreCollisionBC();
   void applyPostCollisionBC();
protected:
private:
   std::vector<BoundaryConditionPtr> preBC;
   std::vector<BoundaryConditionPtr> postBC;
};

#endif // BoundaryConditionProcessor_h__
