#ifndef Block3DVisitor_h
#define Block3DVisitor_h

#include <boost/shared_ptr.hpp>
class Block3DVisitor;
typedef boost::shared_ptr<Block3DVisitor> Block3DVisitorPtr;

#include "Grid3D.h"
#include "Block3D.h"

class Block3DVisitor
{
public:
   Block3DVisitor() : startLevel(-1), stopLevel(-1)
   {
   }

   Block3DVisitor(int startLevel, int stopLevel) : startLevel(startLevel), stopLevel(stopLevel)
   {
   }

	virtual ~Block3DVisitor()
   {
   }
	
   virtual void visit(Grid3DPtr grid, Block3DPtr block)=0;
   
   int  getStartLevel() const; 
   int  getStopLevel() const;
   void setStartLevel(int level);
   void setStopLevel(int level);

private:
   int  startLevel;
   int  stopLevel;
};
//////////////////////////////////////////////////////////////////////////
inline int  Block3DVisitor::getStartLevel() const
{ 
   return this->startLevel;  
}
//////////////////////////////////////////////////////////////////////////
inline int  Block3DVisitor::getStopLevel() const
{ 
   return this->stopLevel;   
}
//////////////////////////////////////////////////////////////////////////
inline void Block3DVisitor::setStartLevel(int level)
{ 
   this->startLevel = level; 
}
//////////////////////////////////////////////////////////////////////////
inline void Block3DVisitor::setStopLevel(int level) 
{ 
   this->stopLevel = level;  
}

#endif 
