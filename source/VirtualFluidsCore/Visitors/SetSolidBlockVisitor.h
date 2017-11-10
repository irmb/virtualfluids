#ifndef SetSolidOrTransBlockVisitor_h
#define SetSolidOrTransBlockVisitor_h

#include "Block3DVisitor.h"
#include "Interactor3D.h"


class SetSolidBlockVisitor : public Block3DVisitor
{
public:
   enum BlockType{SOLID, BC};
public:
   SetSolidBlockVisitor(Interactor3DPtr interactor, BlockType type);

   virtual ~SetSolidBlockVisitor() {}

   virtual void visit(Grid3DPtr grid, Block3DPtr block);

private:
   Interactor3DPtr interactor;
   BlockType type;
};

#endif

