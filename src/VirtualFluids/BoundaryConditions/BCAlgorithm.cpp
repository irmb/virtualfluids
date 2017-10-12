#include "BCAlgorithm.h"

BCAlgorithm::BCAlgorithm() : compressible(false)
{

}
//////////////////////////////////////////////////////////////////////////
void BCAlgorithm::apply()
{
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
//////////////////////////////////////////////////////////////////////////
void BCAlgorithm::addNode(int x1, int x2, int x3)
{
   nodeVector.push_back(x1);
   nodeVector.push_back(x2);
   nodeVector.push_back(x3);
}
//////////////////////////////////////////////////////////////////////////
void BCAlgorithm::addBcPointer(BoundaryConditionsPtr bcPtr)
{
   bcVector.push_back(bcPtr);
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
void BCAlgorithm::clearData()
{
   nodeVector.clear();
   bcVector.clear();
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

