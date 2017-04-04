#ifndef BoundaryConditionBlockVisitor_h__
#define BoundaryConditionBlockVisitor_h__

#include "Block3DVisitor.h"
#include "BCAdapter.h"

class BoundaryConditionsBlockVisitor : public Block3DVisitor
{
public:
   BoundaryConditionsBlockVisitor();
   virtual ~BoundaryConditionsBlockVisitor();
   
   void visit(Grid3DPtr grid, Block3DPtr block);
   void addBC(BCAdapterPtr bc);
protected:
private:
   std::map<char, BCAlgorithmPtr> bcMap;
};
#endif // BoundaryConditionBlockVisitor_h__
