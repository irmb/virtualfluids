#ifndef SetSolidBlocksBlockVisitor_h__
#define SetSolidBlocksBlockVisitor_h__

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;
class Interactor3D;

class SetSolidBlocksBlockVisitor : public Block3DVisitor
{
public:
   SetSolidBlocksBlockVisitor(SPtr<Interactor3D> interactor);
   virtual ~SetSolidBlocksBlockVisitor() {}

   virtual void visit(SPtr<Grid3D> grid, SPtr<Block3D> block);

private:
   SPtr<Interactor3D> interactor;
}
#endif // SetSolidBlocksBlockVisitor_h__
;
