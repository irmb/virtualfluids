#ifndef SetSolidOrTransBlockVisitor_h
#define SetSolidOrTransBlockVisitor_h

#include "Block3DVisitor.h"
#include "Interactor3D.h"


class SetSolidOrTransBlockVisitor : public Block3DVisitor
{
public:
   enum BlockType{SOLID, TRANS};
public:
   SetSolidOrTransBlockVisitor(Interactor3DPtr interactor, BlockType type);

   virtual ~SetSolidOrTransBlockVisitor() {}

   virtual void visit(Grid3DPtr grid, Block3DPtr block);

private:
   Interactor3DPtr interactor;
   BlockType type;
};

#endif

