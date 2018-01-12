#include "SetSolidBlockVisitor.h"

#include "Interactor3D.h"
#include "Grid3DSystem.h"
#include "Grid3D.h"
#include "Block3D.h"

SetSolidBlockVisitor::SetSolidBlockVisitor(SPtr<Interactor3D> interactor, BlockType type) : 
   Block3DVisitor(0, Grid3DSystem::MAXLEVEL), interactor(interactor),
   type(type)
{

}

void SetSolidBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
   if(block->getRank() == grid->getRank())
   {
      if (block->isActive())
      {
         switch (type)
         {
         case BlockType::SOLID:
            interactor->setSolidBlock(block);
            break;
         case BlockType::BC:
            interactor->setBCBlock(block);
            break;
         }
      }
   }
}

