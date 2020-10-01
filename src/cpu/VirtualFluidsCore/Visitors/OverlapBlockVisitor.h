#ifndef OverlapBlockVisitor_H
#define OverlapBlockVisitor_H

#include <string>
#include <PointerDefinitions.h>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;

class OverlapBlockVisitor : public Block3DVisitor
{
public:
   OverlapBlockVisitor(int levelDepth, bool includeNotActiveBlocks = true);
   
   virtual ~OverlapBlockVisitor(){}

   bool isIterative()   { return false; }

   std::string getSpecificDescription();

   void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;


private:
   int  levelDepth;
   bool includeNotActiveBlocks;
};

#endif //OverlapBlockVisitor_H
