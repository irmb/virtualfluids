#include "SetSolidOrBoundaryBlockVisitor.h"

#include "Interactor3D.h"
#include "Grid3DSystem.h"
#include "Grid3D.h"
#include "Block3D.h"

SetSolidOrBoundaryBlockVisitor::SetSolidOrBoundaryBlockVisitor(SPtr<Interactor3D> interactor, BlockType type) : 
   Block3DVisitor(0, Grid3DSystem::MAXLEVEL), interactor(interactor),
   type(type)
{

}

void SetSolidOrBoundaryBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
   if(block->getRank() == grid->getRank())
   {
      if (block->isActive())
      {
         switch (type)
         {
         case SOLID:
            interactor->setSolidBlock(block);
            break;
         case BC:
            interactor->setBCBlock(block);
            break;
         }
      }
   }
}

