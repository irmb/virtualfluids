#include "SpongeLayerBlockVisitor.h"

#include <vector>

#include "Grid3DSystem.h"
#include "LBMSystem.h"
#include <numerics/geometry3d/GbCuboid3D.h>
#include "Grid3D.h"
#include "Block3D.h"
#include "ILBMKernel.h"


SpongeLayerBlockVisitor::SpongeLayerBlockVisitor(GbCuboid3DPtr boundingBox) : Block3DVisitor(0, Grid3DSystem::MAXLEVEL), boundingBox(boundingBox)
{

}
//////////////////////////////////////////////////////////////////////////
SpongeLayerBlockVisitor::~SpongeLayerBlockVisitor()
{

}
//////////////////////////////////////////////////////////////////////////
void SpongeLayerBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
{
    if (block->getRank() == grid->getRank())
    {
        UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);
        UbTupleDouble3 blockLengths = grid->getBlockLengths(block);

        double minX1 = val<1>(org);
        double minX2 = val<2>(org);
        double minX3 = val<3>(org);
        double maxX1 = val<1>(org) + val<1>(blockLengths);
        double maxX2 = val<2>(org) + val<2>(blockLengths);
        double maxX3 = val<3>(org) + val<3>(blockLengths);

        if (boundingBox->isCellInsideGbObject3D(minX1, minX2, minX3, maxX1, maxX2, maxX3))
        {
            ILBMKernelPtr kernel = block->getKernel();
            kernel->setCollisionFactor(LBMSystem::calcCollisionFactor(0.01, block->getLevel()));
        }
    }
}

