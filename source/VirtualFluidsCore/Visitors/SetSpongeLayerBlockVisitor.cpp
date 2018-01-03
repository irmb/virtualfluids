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
void SetSpongeLayerBlockVisitor::visit( Grid3DPtr grid, Block3DPtr block )
{
   if(block->getRank() == grid->getRank())
   {
       LBMKernelPtr kernel = std::dynamic_pointer_cast<LBMKernel>(block->getKernel());
       if (!kernel)
           throw std::runtime_error("SetSpongeLayerBlockVisitor: Kernel is not a LBMKernel");
      kernel->setWithSpongeLayer(true);
      kernel->setSpongeLayer(spongeLayer);
   }
}


