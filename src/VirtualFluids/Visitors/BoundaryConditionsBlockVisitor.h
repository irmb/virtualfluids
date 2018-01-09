#ifndef BoundaryConditionBlockVisitor_h__
#define BoundaryConditionBlockVisitor_h__

#include <map>
#include <memory>

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
   
      void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block) override;
   void addBC(std::shared_ptr<BCAdapter> bc);
protected:
private:
   std::map<char, std::shared_ptr<BCAlgorithm> > bcMap;
};
#endif // BoundaryConditionBlockVisitor_h__
