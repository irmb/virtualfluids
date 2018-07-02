#include "SetBcBlocksBlockVisitor.h"

#include "Interactor3D.h"
#include "Grid3DSystem.h"
#include "Grid3D.h"
#include "Block3D.h"

SetBcBlocksBlockVisitor::SetBcBlocksBlockVisitor(SPtr<Interactor3D> interactor) : 
   Block3DVisitor(0, Grid3DSystem::MAXLEVEL), interactor(interactor)
{

}

void SetBcBlocksBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
   if(block->getRank() == grid->getRank())
   {
      if (block->isActive())
      {
         interactor->setBCBlock(block);
      }
   }
}

