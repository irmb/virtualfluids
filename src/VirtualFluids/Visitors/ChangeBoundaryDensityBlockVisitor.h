#ifndef ChangeBoundaryDensityBlockVisitor_h__
#define ChangeBoundaryDensityBlockVisitor_h__

#include <memory>

#include "Block3DVisitor.h"

class Block3D;
class Grid3D;
class BoundaryConditions;

class ChangeBoundaryDensityBlockVisitor : public Block3DVisitor
{
public:
   ChangeBoundaryDensityBlockVisitor(float oldBoundaryDensity, float newBoundaryDensity);
   virtual ~ChangeBoundaryDensityBlockVisitor();

   void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block) override;
private:
   float oldBoundaryDensity; 
   float newBoundaryDensity;
   std::shared_ptr<BoundaryConditions> bcPtr;
};
#endif // ChangeBoundaryDensityBlockVisitor_h__
