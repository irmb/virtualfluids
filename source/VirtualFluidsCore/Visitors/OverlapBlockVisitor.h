#ifndef OverlapBlockVisitor_H
#define OverlapBlockVisitor_H

#include <string>
#include <memory>

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

   void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block) override;


private:
   int  levelDepth;
   bool includeNotActiveBlocks;
};

#endif //OverlapBlockVisitor_H
