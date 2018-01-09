#ifndef Grid3DVisitor_h
#define Grid3DVisitor_h

#include <memory>

class Grid3DVisitor;
typedef std::shared_ptr<Grid3DVisitor> Grid3DVisitorPtr;

class Grid3D;

class Grid3DVisitor
{
public:
   Grid3DVisitor() {}
   virtual ~Grid3DVisitor() {}

   virtual void visit(std::shared_ptr<Grid3D> grid) = 0;
};

#endif 
