#include "OverlapBlockVisitor.h"
#include "Grid3DSystem.h"

OverlapBlockVisitor::OverlapBlockVisitor(int levelDepth/*shut be maxGridLevel*/, bool includeNotActiveBlocks)
   :   Block3DVisitor(0, Grid3DSystem::MAXLEVEL)
   , levelDepth(levelDepth)
   , includeNotActiveBlocks(includeNotActiveBlocks)
{
}
//////////////////////////////////////////////////////////////////////////
void OverlapBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
{
   int ix1, ix2, ix3, level;
   ix1 = block->getX1();
   ix2 = block->getX2();
   ix3 = block->getX3();
   level = block->getLevel();

   int nix1, nix2,nix3, nlev;
   int neighix1, neighix2, neighix3, neighlev;
   std::vector<Block3DPtr> neighbors;
   grid->getAllNeighbors(ix1, ix2, ix3, level, this->levelDepth, neighbors);
   bool hasAdded = false;
   for(size_t i=0; i<neighbors.size(); i++)
   {
      if(   ( neighbors[i]->isActive() || includeNotActiveBlocks )
         && neighbors[i]->getLevel() > level) 
      {
         neighix1 = neighbors[i]->getX1();
         neighix2 = neighbors[i]->getX2();
         neighix3 = neighbors[i]->getX3();
         neighlev = neighbors[i]->getLevel();
         nix1 = neighix1>>1;
         nix2 = neighix2>>1;
         nix3 = neighix3>>1;
         nlev = neighlev-1;

         if(nlev != level) 
         {
            throw UbException(UB_EXARGS, "OverlapBlockVisitor::adaptBlock - leveldifferenz passt nicht, block: " + block->toString());
         }

         Block3DPtr newBlock = grid->getBlock(nix1,nix2,nix3,nlev);
         if(!newBlock)
         {
            newBlock = Block3DPtr(new Block3D(nix1,nix2,nix3,nlev));
            grid->addBlock(newBlock);
            hasAdded=true;
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
std::string OverlapBlockVisitor::getSpecificDescription()
{
   std::string str("Overlap:");
   return str;
}
