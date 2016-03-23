#ifndef RatioSmoothBlockVisitor_H
#define RatioSmoothBlockVisitor_H

#include "Block3DVisitor.h"

class RatioSmoothBlockVisitor : public Block3DVisitor
{
public:
   RatioSmoothBlockVisitor(int levelDepth, bool includeNotActiveBlocks = false);

   virtual ~RatioSmoothBlockVisitor() {}

   bool expandsByAdaptation() { return this->expandBlocks; }

   void setExpandByAdaptation(bool expandBlocks);

   int  getLevelRatio() { return this->maxLevelRatio; }
   bool isIterative()   { return true;                }

   void setLevelRatio(int ratio);

   int  getStartLevel();
   int  getStopLevel();

   void setStartLevel(int level);
   void setStopLevel(int level);

   std::string getSpecificDescription();

   void visit(Grid3DPtr grid, Block3DPtr block);

protected:
   bool lookForExpand(Grid3DPtr grid, const int& ix1, const int& ix2, const int& ix3, const int& level);
   bool lookForCollapse(Grid3DPtr grid, const int& ix1, const int& ix2, const int& ix3, const int& level);

private:
   int  maxLevelRatio;
   bool expandBlocks;
   int  levelDepth;
   bool includeNotActiveBlocks;

};

#endif //RatioSmoothBlockVisitor_H
