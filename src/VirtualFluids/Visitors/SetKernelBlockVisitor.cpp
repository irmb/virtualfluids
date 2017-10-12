#include "SetKernelBlockVisitor.h"
#include "Grid3DSystem.h"
#include "LBMSystem.h"

//SetKernelBlockVisitor::SetKernelBlockVisitor(LBMKernel3DPtr kernel, LBMReal nue) : 
//                        Block3DVisitor(0, Grid3DSystem::MAXLEVEL), kernel(kernel), nue(nue)
//{
//
//}
//////////////////////////////////////////////////////////////////////////
//SetKernelBlockVisitor::SetKernelBlockVisitor( LBMKernel3DPtr kernel, LBMReal nue, double availMem, double needMem ) : 
//                                              Block3DVisitor(0, Grid3DSystem::MAXLEVEL), kernel(kernel), nue(nue)
//{
//   if (needMem > availMem)
//   {
//      throw UbException(UB_EXARGS,"SetKernelBlockVisitor: Not enough memory!!!");
//   }
//}
//////////////////////////////////////////////////////////////////////////
SetKernelBlockVisitor::SetKernelBlockVisitor(LBMKernelPtr kernel, LBMReal nue, double availMem, double needMem, SetKernelBlockVisitor::Action action /*= SetKernelBlockVisitor::New*/) :
                                             Block3DVisitor(0, Grid3DSystem::MAXLEVEL), kernel(kernel), nue(nue), action(action), dataSetFlag(true)
{
   if (needMem > availMem)
   {
      throw UbException(UB_EXARGS, "SetKernelBlockVisitor: Not enough memory!!!");
   }
}
//////////////////////////////////////////////////////////////////////////
void SetKernelBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
{
   if(kernel && (block->getRank() == grid->getRank()))
   {
      LBMReal collFactor = LBMSystem::calcCollisionFactor(nue, block->getLevel());
      kernel->setCollisionFactor(collFactor);
      kernel->setIndex(block->getX1(), block->getX2(), block->getX3());
      kernel->setDeltaT(LBMSystem::getDeltaT(block->getLevel()));
      kernel->setBlock(block);
      LBMKernelPtr newKernel = kernel->clone();

      switch (action)
      {
      case SetKernelBlockVisitor::NewKernel:
         block->setKernel(newKernel);
         break;
      case SetKernelBlockVisitor::ChangeKernel:
      {
         DataSet3DPtr dataSet = block->getKernel()->getDataSet();
         if (!dataSet)
         {
            UB_THROW(UbException(UB_EXARGS, "It is not possible to change a DataSet in kernel! Old DataSet is not exist!"));
         }

         newKernel->setDataSet(dataSet);
         
         BCProcessorPtr bcProc = block->getKernel()->getBCProcessor();
         if (!bcProc)
         {
            UB_THROW(UbException(UB_EXARGS, "It is not possible to change a BCProcessor in kernel! Old BCProcessor is not exist!"));
         }
         newKernel->setBCProcessor(bcProc);
         block->setKernel(newKernel);
      }
         break;

      case SetKernelBlockVisitor::ChangeKernelWithData:
      {
         BCProcessorPtr bcProc = block->getKernel()->getBCProcessor();
         if (!bcProc)
         {
            UB_THROW(UbException(UB_EXARGS, "It is not possible to change a BCProcessor in kernel! Old BCProcessor is not exist!"));
         }
         newKernel->setBCProcessor(bcProc);
         block->setKernel(newKernel);
      }
         break;
      }

   }
}

void SetKernelBlockVisitor::setNoDataSetFlag(bool flag)
{
   dataSetFlag = flag;
}

