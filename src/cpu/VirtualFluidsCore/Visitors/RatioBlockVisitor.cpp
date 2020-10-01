#include "RatioBlockVisitor.h"
#include "Grid3DSystem.h"
#include "Grid3D.h"
#include "Block3D.h"

RatioBlockVisitor::RatioBlockVisitor(int levelDepth, bool includeNotActiveBlocks)
   :   Block3DVisitor(0, Grid3DSystem::MAXLEVEL)
   , maxLevelRatio(1)
   , expandBlocks(true)
   , levelDepth(levelDepth)
   , includeNotActiveBlocks(includeNotActiveBlocks)
{

}
//////////////////////////////////////////////////////////////////////////
void RatioBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
   int ix1, ix2, ix3, level;
   ix1 = block->getX1();
   ix2 = block->getX2();
   ix3 = block->getX3();
   level = block->getLevel();

   if( block->isActive()  || includeNotActiveBlocks )
   {
      if(this->expandBlocks)
      {
         if(this->lookForExpand(grid,ix1, ix2, ix3, level))
         {
            grid->expandBlock(ix1, ix2, ix3, level);
         }
      }
      else
      {
         if(this->lookForCollapse(grid,ix1, ix2, ix3, level))
         {
            grid->collapseBlock(ix1, ix2, ix3, level, levelDepth);
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
bool RatioBlockVisitor::lookForExpand(SPtr<Grid3D> grid, const int& ix1, const int& ix2, const int& ix3, const int& level)
{
   std::vector<SPtr<Block3D>> neighbors;
   grid->getAllNeighbors(ix1, ix2, ix3, level, this->levelDepth, neighbors);
   for(size_t i=0; i<neighbors.size(); i++)
   {
      if(   ( neighbors[i]->isActive() || includeNotActiveBlocks )
         && neighbors[i]->getLevel() > level+this->maxLevelRatio) 
      {
         return true;
      }
   }

   return false;
}
//////////////////////////////////////////////////////////////////////////
bool RatioBlockVisitor::lookForCollapse(SPtr<Grid3D> grid, const int& ix1, const int& ix2, const int& ix3, const int& level)
{
   std::vector<SPtr<Block3D>> neighbors;
   grid->getAllNeighbors(ix1, ix2,ix3, level, this->levelDepth, neighbors);
   for(size_t i=0; i<neighbors.size(); i++)
   {     
      if(    ( neighbors[i]->isActive() || includeNotActiveBlocks )
         &&  neighbors[i]->getLevel() < level-this->maxLevelRatio) 
      {
         return true;
      }
   }

   return false;
}
//////////////////////////////////////////////////////////////////////////
void RatioBlockVisitor::setExpandByAdaptation(bool expandBlocks)
{
   if(this->expandBlocks != expandBlocks)
   {
      this->expandBlocks = expandBlocks;

      int l1 = Block3DVisitor::getStartLevel();
      int l2 = Block3DVisitor::getStopLevel();

      if(expandBlocks) { if(l1 < l2) { Block3DVisitor::setStartLevel(l2); Block3DVisitor::setStopLevel(l1); } }
      else             { if(l2 < l1) { Block3DVisitor::setStartLevel(l2); Block3DVisitor::setStopLevel(l1); } }
   }
}
//////////////////////////////////////////////////////////////////////////
void RatioBlockVisitor::setLevelRatio(int ratio)
{
   if(ratio < 1) throw UbException(UB_EXARGS,"illegal ratio specified");
   this->maxLevelRatio = ratio;
}
//////////////////////////////////////////////////////////////////////////
std::string RatioBlockVisitor::getSpecificDescription()
{
   std::string str("Ratio:");
   return str;
}
//////////////////////////////////////////////////////////////////////////
int RatioBlockVisitor::getStartLevel()
{
   int l1 = Block3DVisitor::getStartLevel();
   int l2 = Block3DVisitor::getStopLevel();

   if(this->expandBlocks) { if(l2 < l1) return(l1); else return(l2); }
   else                   { if(l2 < l1) return(l2); else return(l1); }
}
//////////////////////////////////////////////////////////////////////////
int RatioBlockVisitor::getStopLevel()
{
   int l1 = Block3DVisitor::getStartLevel();
   int l2 = Block3DVisitor::getStopLevel();

   if(this->expandBlocks) { if(l2 < l1) return(l2); else return(l1); }
   else                   { if(l2 < l1) return(l1); else return(l2); }
}
//////////////////////////////////////////////////////////////////////////
void RatioBlockVisitor::setStartLevel(int level)
{
   if(this->expandBlocks) { if(level >= Block3DVisitor::getStopLevel()) Block3DVisitor::setStartLevel(level); }
   else                   { if(level <= Block3DVisitor::getStopLevel()) Block3DVisitor::setStartLevel(level); }
}
//////////////////////////////////////////////////////////////////////////
void RatioBlockVisitor::setStopLevel(int level)
{
   if(this->expandBlocks) { if(level <= Block3DVisitor::getStartLevel()) Block3DVisitor::setStopLevel(level); }
   else                   { if(level >= Block3DVisitor::getStartLevel()) Block3DVisitor::setStopLevel(level); }
}
