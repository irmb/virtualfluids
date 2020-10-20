#include "ChangeBoundaryDensityBlockVisitor.h"
#include "BCArray3D.h"
#include "BCProcessor.h"
#include "Block3D.h"
#include "BoundaryConditions.h"
#include "Grid3D.h"
#include "Grid3DSystem.h"
#include "LBMKernel.h"

ChangeBoundaryDensityBlockVisitor::ChangeBoundaryDensityBlockVisitor(float oldBoundaryDensity, float newBoundaryDensity)
    : Block3DVisitor(0, Grid3DSystem::MAXLEVEL), oldBoundaryDensity(oldBoundaryDensity),
      newBoundaryDensity(newBoundaryDensity)
{
}
//////////////////////////////////////////////////////////////////////////
ChangeBoundaryDensityBlockVisitor::~ChangeBoundaryDensityBlockVisitor() = default;
//////////////////////////////////////////////////////////////////////////
void ChangeBoundaryDensityBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    if (block->getRank() == grid->getRank()) {
        SPtr<ILBMKernel> kernel = block->getKernel();
        SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();

        int minX1 = 0;
        int minX2 = 0;
        int minX3 = 0;
        int maxX1 = (int)bcArray->getNX1();
        int maxX2 = (int)bcArray->getNX2();
        int maxX3 = (int)bcArray->getNX3();

        for (int x3 = minX3; x3 < maxX3; x3++) {
            for (int x2 = minX2; x2 < maxX2; x2++) {
                for (int x1 = minX1; x1 < maxX1; x1++) {
                    if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {
                        bcPtr = bcArray->getBC(x1, x2, x3);
                        if (bcPtr) {
                            if (bcPtr->hasDensityBoundary()) {
                                float bcDensity = (float)bcPtr->getBoundaryDensity();
                                if (bcDensity == oldBoundaryDensity) {
                                    bcPtr->setBoundaryDensity(newBoundaryDensity);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
