#ifndef CheckRatioBlockVisitor_H
#define CheckRatioBlockVisitor_H

#include <string>
#include <memory>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;

class CheckRatioBlockVisitor : public Block3DVisitor
{
public:
   CheckRatioBlockVisitor(int levelDepth, bool includeNotActiveBlocks = true);

   virtual ~CheckRatioBlockVisitor() {}

   bool getState();
   void resetState();
   std::string getStateString();

      void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block) override;

private:
   int  levelDepth;
   bool includeNotActiveBlocks;
   bool state;
   std::shared_ptr<Block3D> falseBlock;
};

#endif //OverlapBlockVisitor_H

