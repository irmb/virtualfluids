#include "BCAlgorithm.h"

BCAlgorithm::BCAlgorithm() : compressible(false)
{

}
//////////////////////////////////////////////////////////////////////////
void BCAlgorithm::setNodeIndex(int x1, int x2, int x3)
{
   this->x1 = x1;
   this->x2 = x2;
   this->x3 = x3;
}
//////////////////////////////////////////////////////////////////////////
void BCAlgorithm::setBcPointer(BoundaryConditionsPtr bcPtr)
{
   this->bcPtr = bcPtr;
}
//////////////////////////////////////////////////////////////////////////
void BCAlgorithm::setCompressible(bool c)
{
   compressible = c;

   if (this->compressible)
   {
      calcFeqsForDirFct = &D3Q27System::getCompFeqForDirection;
      calcMacrosFct = &D3Q27System::calcCompMacroscopicValues;
      calcFeqFct = &D3Q27System::calcCompFeq;
      compressibleFactor = 1.0;
   }
   else
   {
      calcFeqsForDirFct = &D3Q27System::getIncompFeqForDirection;
      calcMacrosFct = &D3Q27System::calcIncompMacroscopicValues;
      calcFeqFct = &D3Q27System::calcIncompFeq;
      compressibleFactor = 0.0;
   }
}
//////////////////////////////////////////////////////////////////////////
void BCAlgorithm::setCollFactor(LBMReal cf)
{
   collFactor = cf;
}

//////////////////////////////////////////////////////////////////////////
char BCAlgorithm::getType()
{
   return type;
}
//////////////////////////////////////////////////////////////////////////
bool BCAlgorithm::isPreCollision()
{
   return preCollision;
}
//////////////////////////////////////////////////////////////////////////
BCArray3DPtr BCAlgorithm::getBcArray()
{
   return bcArray;
}
//////////////////////////////////////////////////////////////////////////
void BCAlgorithm::setBcArray(BCArray3DPtr bcarray)
{
   bcArray = bcarray;
}

