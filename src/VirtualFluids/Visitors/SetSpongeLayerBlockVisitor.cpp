#include "SetSpongeLayerBlockVisitor.h"
#include "Grid3DSystem.h"
#include "LBMSystem.h"

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
      block->getKernel()->setWithSpongeLayer(true);
      block->getKernel()->setSpongeLayer(spongeLayer);
   }
}


