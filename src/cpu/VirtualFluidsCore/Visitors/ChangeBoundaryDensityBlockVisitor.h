#ifndef ChangeBoundaryDensityBlockVisitor_h__
#define ChangeBoundaryDensityBlockVisitor_h__

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"

class Block3D;
class Grid3D;
class BoundaryConditions;

class ChangeBoundaryDensityBlockVisitor : public Block3DVisitor
{
public:
   ChangeBoundaryDensityBlockVisitor(float oldBoundaryDensity, float newBoundaryDensity);
   virtual ~ChangeBoundaryDensityBlockVisitor();

   void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;
private:
   float oldBoundaryDensity; 
   float newBoundaryDensity;
   SPtr<BoundaryConditions> bcPtr;
};
#endif // ChangeBoundaryDensityBlockVisitor_h__
