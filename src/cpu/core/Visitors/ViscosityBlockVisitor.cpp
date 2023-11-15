#include "ViscosityBlockVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"
#include "D3Q27System.h"
#include "ILBMKernel.h"
#include "LBMSystem.h"

ViscosityBlockVisitor::ViscosityBlockVisitor(real nu) : Block3DVisitor(0, D3Q27System::MAXLEVEL), nu(nu) {}
//////////////////////////////////////////////////////////////////////////
void ViscosityBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    if (block->getRank() == grid->getRank()) {
        real collFactor = LBMSystem::calcCollisionFactor(nu, block->getLevel());
        block->getKernel()->setCollisionFactor(collFactor);
    }
}
