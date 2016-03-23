#include "SetSolidOrTransBlockVisitor.h"
#include "Grid3DSystem.h"

SetSolidOrTransBlockVisitor::SetSolidOrTransBlockVisitor(Interactor3DPtr interactor, BlockType type) : 
   Block3DVisitor(0, Grid3DSystem::MAXLEVEL), interactor(interactor),
   type(type)
{


}
//////////////////////////////////////////////////////////////////////////
void SetSolidOrTransBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
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
         case TRANS:
            interactor->setTransBlock(block);
            break;
         }
      }
   }
}

