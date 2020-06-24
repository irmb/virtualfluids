#ifndef BoundaryConditionBlockVisitor_h__
#define BoundaryConditionBlockVisitor_h__

#include <map>
#include <PointerDefinitions.h>

#include "Block3DVisitor.h"


class Grid3D;
class Block3D;
class BCAlgorithm;
class BCAdapter;

class BoundaryConditionsBlockVisitor : public Block3DVisitor
{
public:
   BoundaryConditionsBlockVisitor();
   virtual ~BoundaryConditionsBlockVisitor();
   
      void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;
   void addBC(SPtr<BCAdapter> bc);
protected:
private:
   std::map<char, SPtr<BCAlgorithm> > bcMap;
};
#endif // BoundaryConditionBlockVisitor_h__
