#ifndef ChangeBoundaryDensityBlockVisitor_h__
#define ChangeBoundaryDensityBlockVisitor_h__

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"
#include "lbm/constants/D3Q27.h"

class Block3D;
class Grid3D;
class BoundaryConditions;

class ChangeBoundaryDensityBlockVisitor : public Block3DVisitor
{
public:
    ChangeBoundaryDensityBlockVisitor(real oldBoundaryDensity, real newBoundaryDensity);
    ~ChangeBoundaryDensityBlockVisitor() override;

    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

private:
    real oldBoundaryDensity;
    real newBoundaryDensity;
    SPtr<BoundaryConditions> bcPtr;
};
#endif // ChangeBoundaryDensityBlockVisitor_h__
