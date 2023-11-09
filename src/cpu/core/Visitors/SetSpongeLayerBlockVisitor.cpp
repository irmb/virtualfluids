#include "SetSpongeLayerBlockVisitor.h"
#include "D3Q27System.h"
#include "LBMSystem.h"

#include "Block3D.h"
#include "Grid3D.h"
#include "LBMKernel.h"

SetSpongeLayerBlockVisitor::SetSpongeLayerBlockVisitor(const mu::Parser &spongeLayer)
    : Block3DVisitor(0, D3Q27System::MAXLEVEL), spongeLayer(spongeLayer)
{
}
//////////////////////////////////////////////////////////////////////////
SetSpongeLayerBlockVisitor::~SetSpongeLayerBlockVisitor() = default;
//////////////////////////////////////////////////////////////////////////
void SetSpongeLayerBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    if (block->getRank() == grid->getRank()) {
        SPtr<LBMKernel> kernel = dynamicPointerCast<LBMKernel>(block->getKernel());
        if (!kernel)
            throw std::runtime_error("SetSpongeLayerBlockVisitor: Kernel is not a LBMKernel");
        kernel->setWithSpongeLayer(true);
        kernel->setSpongeLayer(spongeLayer);
    }
}
