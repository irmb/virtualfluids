#ifndef Grid3DVisitor_h
#define Grid3DVisitor_h

#include <PointerDefinitions.h>


class Grid3D;

class Grid3DVisitor
{
public:
   Grid3DVisitor() {}
   virtual ~Grid3DVisitor() {}

   virtual void visit(SPtr<Grid3D> grid) = 0;
};

#endif 
