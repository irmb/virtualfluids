#include "ViscosityBlockVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"
#include "Grid3DSystem.h"
#include "ILBMKernel.h"
#include "LBMSystem.h"

ViscosityBlockVisitor::ViscosityBlockVisitor(LBMReal nu) : Block3DVisitor(0, Grid3DSystem::MAXLEVEL), nu(nu) {}
//////////////////////////////////////////////////////////////////////////
void ViscosityBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    if (block->getRank() == grid->getRank()) {
        LBMReal collFactor = LBMSystem::calcCollisionFactor(nu, block->getLevel());
        block->getKernel()->setCollisionFactor(collFactor);
    }
}
