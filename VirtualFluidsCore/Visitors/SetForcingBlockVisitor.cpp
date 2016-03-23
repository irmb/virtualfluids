#include "SetForcingBlockVisitor.h"
#include "Grid3DSystem.h"
#include "LBMSystem.h"

SetForcingBlockVisitor::SetForcingBlockVisitor(LBMReal forcingX1, LBMReal forcingX2, LBMReal forcingX3) : 
                        Block3DVisitor(0, Grid3DSystem::MAXLEVEL), forcingX1(forcingX1), 
                                                                   forcingX2(forcingX2),
                                                                   forcingX3(forcingX3)
{
   ftype = 0;
}
//////////////////////////////////////////////////////////////////////////
SetForcingBlockVisitor::SetForcingBlockVisitor(const mu::Parser& muForcingX1, const mu::Parser& muForcingX2, const mu::Parser& muForcingX3) : 
                                              Block3DVisitor(0, Grid3DSystem::MAXLEVEL), muForcingX1(muForcingX1),
                                                                                         muForcingX2(muForcingX2),
                                                                                         muForcingX3(muForcingX3)

{
   ftype = 1;
}
//////////////////////////////////////////////////////////////////////////
SetForcingBlockVisitor::SetForcingBlockVisitor(const std::string& sForcingX1, const std::string& sForcingX2, const std::string& sForcingX3) : 
                                             Block3DVisitor(0, Grid3DSystem::MAXLEVEL), sForcingX1(sForcingX1),
                                                                                        sForcingX2(sForcingX2),
                                                                                        sForcingX3(sForcingX3)

{
   ftype = 2;
}
//////////////////////////////////////////////////////////////////////////
void SetForcingBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
{
   if(block->getRank() == grid->getRank())
   {
      switch (ftype)
      {
      case 0:
         block->getKernel()->setForcingX1(forcingX1);
         block->getKernel()->setForcingX2(forcingX2);
         block->getKernel()->setForcingX3(forcingX3);
         block->getKernel()->setWithForcing(true);
         break;
      case 1:
         block->getKernel()->setForcingX1(muForcingX1);
         block->getKernel()->setForcingX2(muForcingX2);
         block->getKernel()->setForcingX3(muForcingX3);
         block->getKernel()->setWithForcing(true);
         break;
      case 2:
         block->getKernel()->setForcingX1(sForcingX1);
         block->getKernel()->setForcingX2(sForcingX2);
         block->getKernel()->setForcingX3(sForcingX3);
         block->getKernel()->setWithForcing(true);
         break;
      default:
         block->getKernel()->setForcingX1(0.0);
         block->getKernel()->setForcingX2(0.0);
         block->getKernel()->setForcingX3(0.0);
         block->getKernel()->setWithForcing(false);
         break;
      }
   }
}

