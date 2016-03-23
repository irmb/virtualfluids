#ifndef Grid3DVisitor_h
#define Grid3DVisitor_h

#include <boost/shared_ptr.hpp>
class Grid3DVisitor;
typedef boost::shared_ptr<Grid3DVisitor> Grid3DVisitorPtr;

#include "Grid3D.h" 

class Grid3DVisitor
{
public:
   Grid3DVisitor() {}
   virtual ~Grid3DVisitor() {}

   virtual void visit(Grid3DPtr grid)=0;
};

#endif 
