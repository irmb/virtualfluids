#include "SetSolidBlockVisitor.h"
#include "Grid3DSystem.h"

SetSolidBlockVisitor::SetSolidBlockVisitor(Interactor3DPtr interactor, BlockType type) : 
   Block3DVisitor(0, Grid3DSystem::MAXLEVEL), interactor(interactor),
   type(type)
{


}
//////////////////////////////////////////////////////////////////////////
void SetSolidBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
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

