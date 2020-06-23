#ifndef SetBcBlocksBlockVisitor_h__
#define SetBcBlocksBlockVisitor_h__

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;
class Interactor3D;

class SetBcBlocksBlockVisitor : public Block3DVisitor
{
public:
   SetBcBlocksBlockVisitor(SPtr<Interactor3D> interactor);
   virtual ~SetBcBlocksBlockVisitor() {}

   virtual void visit(SPtr<Grid3D> grid, SPtr<Block3D> block);

private:
   SPtr<Interactor3D> interactor;
};
#endif // SetBcBlocksBlockVisitor_h__



