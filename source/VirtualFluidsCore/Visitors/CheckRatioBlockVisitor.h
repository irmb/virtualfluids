#ifndef CheckRatioBlockVisitor_H
#define CheckRatioBlockVisitor_H

#include "Block3DVisitor.h"

class CheckRatioBlockVisitor : public Block3DVisitor
{
public:
   CheckRatioBlockVisitor(int levelDepth, bool includeNotActiveBlocks = true);

   virtual ~CheckRatioBlockVisitor() {}

   bool getState();
   void resetState();
   std::string getStateString();

   void visit(Grid3DPtr grid, Block3DPtr block);

protected:

private:
   int  levelDepth;
   bool includeNotActiveBlocks;
   bool state;
   Block3DPtr falseBlock;
};

#endif //OverlapBlockVisitor_H

