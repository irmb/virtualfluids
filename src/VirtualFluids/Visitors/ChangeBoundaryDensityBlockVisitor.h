#ifndef ChangeBoundaryDensityBlockVisitor_h__
#define ChangeBoundaryDensityBlockVisitor_h__

#include "Block3DVisitor.h"
#include "BoundaryConditions.h"

class ChangeBoundaryDensityBlockVisitor : public Block3DVisitor
{
public:
   ChangeBoundaryDensityBlockVisitor(float oldBoundaryDensity, float newBoundaryDensity);
   virtual ~ChangeBoundaryDensityBlockVisitor();
   virtual void visit(Grid3DPtr grid, Block3DPtr block);
protected:
private:
   float oldBoundaryDensity; 
   float newBoundaryDensity;
   BoundaryConditionsPtr bcPtr;
};
#endif // ChangeBoundaryDensityBlockVisitor_h__
