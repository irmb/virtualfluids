#include "RenumberBlockVisitor.h"
#include "Grid3DSystem.h"
#include "LBMSystem.h"
#include "Grid3D.h"
#include "Block3D.h"

int RenumberBlockVisitor::counter = 0;

RenumberBlockVisitor::RenumberBlockVisitor() :
Block3DVisitor(0, Grid3DSystem::MAXLEVEL)
{

}
//////////////////////////////////////////////////////////////////////////
void RenumberBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
{
   block->setGlobalID(counter);
   Grid3D::BlockIDMap blockIdMap = grid->getBlockIDs();
   blockIdMap.insert(std::make_pair(counter, block));
   counter++;
}

