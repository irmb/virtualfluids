#ifndef OverlapBlockVisitor_H
#define OverlapBlockVisitor_H

#include "Block3DVisitor.h"

class OverlapBlockVisitor : public Block3DVisitor
{
public:
   OverlapBlockVisitor(int levelDepth, bool includeNotActiveBlocks = true);
   
   virtual ~OverlapBlockVisitor(){}

   bool isIterative()   { return false; }

   std::string getSpecificDescription();

   void visit(Grid3DPtr grid, Block3DPtr block);

protected:

private:
   int  levelDepth;
   bool includeNotActiveBlocks;
};

#endif //OverlapBlockVisitor_H
