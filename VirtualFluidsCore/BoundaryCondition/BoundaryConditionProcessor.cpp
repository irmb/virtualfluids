#include "BoundaryConditionProcessor.h"
#include <boost/foreach.hpp>

BoundaryConditionProcessor::BoundaryConditionProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
void BoundaryConditionProcessor::addBC(BoundaryConditionPtr bc)
{
   BoundaryCondition::Type type = bc->getType();
   if (bc->isPreCollision())
   {
      preBC.push_back(bc);
   }
   else
   {
      postBC.push_back(bc);
   }
}
//////////////////////////////////////////////////////////////////////////
BoundaryConditionPtr BoundaryConditionProcessor::getBC(BoundaryCondition::Type type)
{
   BOOST_FOREACH(BoundaryConditionPtr bc, preBC)
   {
     if (bc->getType() == type)
     {
        return bc;
     }
   }

   BOOST_FOREACH(BoundaryConditionPtr bc, postBC)
   {
      if (bc->getType() == type)
      {
         return bc;
      }
   }

   return BoundaryConditionPtr();
}
//////////////////////////////////////////////////////////////////////////
void BoundaryConditionProcessor::applyPreCollisionBC()
{
   BOOST_FOREACH(BoundaryConditionPtr bc, preBC)
   {
      bc->apply();
   }
}
//////////////////////////////////////////////////////////////////////////
void BoundaryConditionProcessor::applyPostCollisionBC()
{
   BOOST_FOREACH(BoundaryConditionPtr bc, postBC)
   {
      bc->apply();
   }
}

