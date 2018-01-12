#include "SetForcingBlockVisitor.h"
#include "Grid3DSystem.h"
#include "LBMSystem.h"
#include "Grid3D.h"
#include "Block3D.h"

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
void SetForcingBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    SPtr<LBMKernel> kernel = dynamicPointerCast<LBMKernel>(block->getKernel());
    if (!kernel)
        throw std::runtime_error("SetForcingBlockVisitor: Kernel is not a LBMKernel");

   if(block->getRank() == grid->getRank())
   {
      switch (ftype)
      {
      case 0:
         kernel->setForcingX1(forcingX1);
         kernel->setForcingX2(forcingX2);
         kernel->setForcingX3(forcingX3);
         kernel->setWithForcing(true);
         break;
      case 1:
         kernel->setForcingX1(muForcingX1);
         kernel->setForcingX2(muForcingX2);
         kernel->setForcingX3(muForcingX3);
         kernel->setWithForcing(true);
         break;
      case 2:
         kernel->setForcingX1(sForcingX1);
         kernel->setForcingX2(sForcingX2);
         kernel->setForcingX3(sForcingX3);
         kernel->setWithForcing(true);
         break;
      default:
         kernel->setForcingX1(0.0);
         kernel->setForcingX2(0.0);
         kernel->setForcingX3(0.0);
         kernel->setWithForcing(false);
         break;
      }
   }
}

