#include "BoundaryCondition.h"

BoundaryCondition::BoundaryCondition() : compressible(false)
{

}
//////////////////////////////////////////////////////////////////////////
void BoundaryCondition::apply()
{
   if (this->compressible)
   {
      calcFeqsForDirFct = &D3Q27System::getCompFeqForDirection;
      calcMacrosFct = &D3Q27System::calcCompMacroscopicValues;
      calcFeqFct = &D3Q27System::calcCompFeq;
   }
   else
   {
      calcFeqsForDirFct = &D3Q27System::getIncompFeqForDirection;
      calcMacrosFct = &D3Q27System::calcIncompMacroscopicValues;
      calcFeqFct = &D3Q27System::calcIncompFeq;
   }

   int cbc = 0;
   int cn = 0;
   int j;

   int nsize = (int)nodeVector.size();

      for (j = cn; j < nsize;)
      {
         x1 = nodeVector[j++];
         x2 = nodeVector[j++];
         x3 = nodeVector[j++];

         bcPtr = bcVector[cbc];
         cbc++;

         applyBC();
      }
      cn = j;
}
////////////////////////////////////////////////////////////////////////////
//void BoundaryCondition::addDistributions(EsoTwist3DPtr distributions)
//{
//   this->distributions = distributions;
//}
//////////////////////////////////////////////////////////////////////////
void BoundaryCondition::addNode(int x1, int x2, int x3)
{
   nodeVector.push_back(x1);
   nodeVector.push_back(x2);
   nodeVector.push_back(x3);
}
//////////////////////////////////////////////////////////////////////////
void BoundaryCondition::addBcPointer(D3Q27BoundaryConditionPtr bcPtr)
{
   bcVector.push_back(bcPtr);
}
//////////////////////////////////////////////////////////////////////////
void BoundaryCondition::setCompressible(bool c)
{
   compressible = c;
}
//////////////////////////////////////////////////////////////////////////
void BoundaryCondition::setCollFactor(LBMReal cf)
{
   collFactor = cf;
}

//////////////////////////////////////////////////////////////////////////
BoundaryCondition::Type BoundaryCondition::getType()
{
   return type;
}
//////////////////////////////////////////////////////////////////////////
bool BoundaryCondition::isPreCollision()
{
   return preCollision;
}




