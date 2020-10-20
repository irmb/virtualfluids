#include "RenumberBlockVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"
#include "Grid3DSystem.h"
#include "LBMSystem.h"

int RenumberBlockVisitor::counter = 0;

RenumberBlockVisitor::RenumberBlockVisitor() : Block3DVisitor(0, Grid3DSystem::MAXLEVEL) {}
//////////////////////////////////////////////////////////////////////////
void RenumberBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    block->setGlobalID(counter);
    Grid3D::BlockIDMap blockIdMap = grid->getBlockIDs();
    blockIdMap.insert(std::make_pair(counter, block));
    counter++;
}
