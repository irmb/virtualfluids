#include "BCProcessor.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "DataSet3D.h"
#include "ILBMKernel.h"
#include "BCArray3D.h"
#include "BCAlgorithm.h"

BCProcessor::BCProcessor()
{
   
}
//////////////////////////////////////////////////////////////////////////
BCProcessor::BCProcessor(SPtr<ILBMKernel> kernel)
{
   SPtr<DistributionArray3D> distributions = dynamicPointerCast<EsoTwist3D>(kernel->getDataSet()->getFdistributions());
   bcArray = SPtr<BCArray3D>(new BCArray3D( distributions->getNX1(), distributions->getNX2(), distributions->getNX3(), BCArray3D::FLUID));
}
//////////////////////////////////////////////////////////////////////////
BCProcessor::~BCProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCProcessor> BCProcessor::clone(SPtr<ILBMKernel> kernel)
{
   SPtr<BCProcessor> bcProcessor(new BCProcessor(kernel));
   return bcProcessor;
}
//////////////////////////////////////////////////////////////////////////
SPtr<BCArray3D> BCProcessor::getBCArray()
{ 
   return bcArray; 
}
//////////////////////////////////////////////////////////////////////////
void BCProcessor::setBCArray(SPtr<BCArray3D> bcarray)
{
   bcArray = bcarray;
}
//////////////////////////////////////////////////////////////////////////
void BCProcessor::addBC(SPtr<BCAlgorithm> bc)
{
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
void BCProcessor::applyPreCollisionBC()
{
   for(SPtr<BCAlgorithm> bc : preBC)
      bc->applyBC();
}
//////////////////////////////////////////////////////////////////////////
void BCProcessor::applyPostCollisionBC()
{
    for (SPtr<BCAlgorithm> bc : postBC)
        bc->applyBC();
}
//////////////////////////////////////////////////////////////////////////
void BCProcessor::clearBC()
{
   preBC.clear();
   postBC.clear();
}

