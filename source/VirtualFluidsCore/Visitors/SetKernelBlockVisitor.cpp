#include "SetKernelBlockVisitor.h"
#include "Grid3DSystem.h"
#include "LBMSystem.h"
#include "DataSet3D.h"
#include "BCProcessor.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "LBMKernel.h"

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
SetKernelBlockVisitor::SetKernelBlockVisitor(SPtr<LBMKernel> kernel, LBMReal nue, double availMem, double needMem, SetKernelBlockVisitor::Action action /*= SetKernelBlockVisitor::New*/) :
                                             Block3DVisitor(0, Grid3DSystem::MAXLEVEL), kernel(kernel), nue(nue), action(action), dataSetFlag(true)
{
   if (needMem > availMem)
   {
      throw UbException(UB_EXARGS, "SetKernelBlockVisitor: Not enough memory!!!");
   }
}
//////////////////////////////////////////////////////////////////////////
void SetKernelBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
   if(kernel && (block->getRank() == grid->getRank()))
   {
      LBMReal collFactor = LBMSystem::calcCollisionFactor(nue, block->getLevel());
      kernel->setCollisionFactor(collFactor);
      kernel->setIndex(block->getX1(), block->getX2(), block->getX3());
      kernel->setDeltaT(LBMSystem::getDeltaT(block->getLevel()));
      kernel->setBlock(block);
      SPtr<LBMKernel> newKernel = kernel->clone();

      switch (action)
      {
      case SetKernelBlockVisitor::NewKernel:
         block->setKernel(newKernel);
         break;
      case SetKernelBlockVisitor::ChangeKernel:
      {
         SPtr<DataSet3D> dataSet = block->getKernel()->getDataSet();
         if (!dataSet)
         {
            UB_THROW(UbException(UB_EXARGS, "It is not possible to change a DataSet in kernel! Old DataSet is not exist!"));
         }

         newKernel->setDataSet(dataSet);
         
         SPtr<BCProcessor> bcProc = block->getKernel()->getBCProcessor();
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
         SPtr<BCProcessor> bcProc = block->getKernel()->getBCProcessor();
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

