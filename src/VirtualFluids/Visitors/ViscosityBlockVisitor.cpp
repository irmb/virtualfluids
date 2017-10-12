#include "ViscosityBlockVisitor.h"
#include "Grid3DSystem.h"
#include "LBMSystem.h"

ViscosityBlockVisitor::ViscosityBlockVisitor(LBMReal nu) :
Block3DVisitor(0, Grid3DSystem::MAXLEVEL), nu(nu)
{

}
//////////////////////////////////////////////////////////////////////////
void ViscosityBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
{
   if (block->getRank() == grid->getRank())
   {
      LBMReal collFactor = LBMSystem::calcCollisionFactor(nu, block->getLevel());
      block->getKernel()->setCollisionFactor(collFactor);
   }
}

