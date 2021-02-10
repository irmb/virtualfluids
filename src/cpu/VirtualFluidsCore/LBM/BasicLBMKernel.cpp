#include "BasicLBMKernel.h"
#include "D3Q27System.h"
#include "BCArray3D.h"
#include "BCProcessor.h"

BasicLBMKernel::BasicLBMKernel()
{

}

BasicLBMKernel::~BasicLBMKernel(void)
{
}

void BasicLBMKernel::calculate(int step)
{
   using namespace D3Q27System;
   using namespace std;

   //timer.resetAndStart();


   /////////////////////////////////////

   //localDistributions = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getLocalDistributions();
   //nonLocalDistributions = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getNonLocalDistributions();
   //zeroDistributions = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getZeroDistributions();

   initData();

   SPtr<BCArray3D> bcArray = this->getBCProcessor()->getBCArray();

   const int bcArrayMaxX1 = (int)bcArray->getNX1();
   const int bcArrayMaxX2 = (int)bcArray->getNX2();
   const int bcArrayMaxX3 = (int)bcArray->getNX3();

   minX1 = ghostLayerWidth;
   minX2 = ghostLayerWidth;
   minX3 = ghostLayerWidth;
   maxX1 = bcArrayMaxX1 - ghostLayerWidth;
   maxX2 = bcArrayMaxX2 - ghostLayerWidth;
   maxX3 = bcArrayMaxX3 - ghostLayerWidth;

   for (int x3 = minX3; x3 < maxX3; x3++)
   {
      for (int x2 = minX2; x2 < maxX2; x2++)
      {
         for (int x1 = minX1; x1 < maxX1; x1++)
         {
            if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3))
            {
               nodeCollision(step, x1, x2, x3);
            }
         }
      }
   }
}
