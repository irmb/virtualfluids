#ifndef Block3DVisitor_h
#define Block3DVisitor_h

#include <PointerDefinitions.h>

class Block3D;
class Grid3D;

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
	
   virtual void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) = 0;
   
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
