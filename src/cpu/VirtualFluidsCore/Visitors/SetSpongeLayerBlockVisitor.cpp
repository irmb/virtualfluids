#include "SetSpongeLayerBlockVisitor.h"
#include "Grid3DSystem.h"
#include "LBMSystem.h"

#include "LBMKernel.h"
#include "Grid3D.h"
#include "Block3D.h"

SetSpongeLayerBlockVisitor::SetSpongeLayerBlockVisitor( const mu::Parser& spongeLayer ) : Block3DVisitor(0, Grid3DSystem::MAXLEVEL), spongeLayer(spongeLayer)
{

}
//////////////////////////////////////////////////////////////////////////
SetSpongeLayerBlockVisitor::~SetSpongeLayerBlockVisitor()
{

}
//////////////////////////////////////////////////////////////////////////
void SetSpongeLayerBlockVisitor::visit( SPtr<Grid3D> grid, SPtr<Block3D> block )
{
   if(block->getRank() == grid->getRank())
   {
       SPtr<LBMKernel> kernel = dynamicPointerCast<LBMKernel>(block->getKernel());
       if (!kernel)
           throw std::runtime_error("SetSpongeLayerBlockVisitor: Kernel is not a LBMKernel");
      kernel->setWithSpongeLayer(true);
      kernel->setSpongeLayer(spongeLayer);
   }
}


