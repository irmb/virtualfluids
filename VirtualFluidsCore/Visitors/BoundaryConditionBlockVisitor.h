#ifndef BoundaryConditionBlockVisitor_h__
#define BoundaryConditionBlockVisitor_h__

#include "Block3DVisitor.h"
#include "BoundaryConditionProcessor.h"

class BoundaryConditionBlockVisitor : public Block3DVisitor
{
public:
   BoundaryConditionBlockVisitor();
   virtual ~BoundaryConditionBlockVisitor();
   void visit(Grid3DPtr grid, Block3DPtr block);
protected:
private:
   BoundaryConditionPtr velocity;
   BoundaryConditionPtr density;
   BoundaryConditionPtr noSlip;
   BoundaryConditionPtr slip;
};
#endif // BoundaryConditionBlockVisitor_h__
