#include "SpongeLayerBlockVisitor.h"
#include "Grid3DSystem.h"
#include "LBMSystem.h"

using namespace std;

SpongeLayerBlockVisitor::SpongeLayerBlockVisitor( GbCuboid3DPtr boundingBox ) : Block3DVisitor(0, Grid3DSystem::MAXLEVEL), boundingBox(boundingBox)
{

}
//////////////////////////////////////////////////////////////////////////
SpongeLayerBlockVisitor::~SpongeLayerBlockVisitor()
{

}
//////////////////////////////////////////////////////////////////////////
void SpongeLayerBlockVisitor::visit( Grid3DPtr grid, Block3DPtr block )
{
   double minX1,minX2,minX3,maxX1,maxX2,maxX3;
   int gridRank = grid->getRank();

   int minInitLevel = grid->getCoarsestInitializedLevel();
   int maxInitLevel = grid->getFinestInitializedLevel();

   double numSolids = 0.0;
   double numFluids = 0.0;
   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      vector<Block3DPtr> blockVector;
      grid->getBlocks(level, gridRank, blockVector);
      BOOST_FOREACH(Block3DPtr block, blockVector)
      {
         UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);
         UbTupleDouble3 blockLengths = grid->getBlockLengths(block);

         minX1 = val<1>(org);
         minX2 = val<2>(org);
         minX3 = val<3>(org);
         maxX1 = val<1>(org)+val<1>(blockLengths);
         maxX2 = val<2>(org)+val<2>(blockLengths);
         maxX3 = val<3>(org)+val<3>(blockLengths);

         if (boundingBox->isCellInsideGbObject3D(minX1,minX2,minX3,maxX1,maxX2,maxX3))
         {
            LBMKernelPtr kernel = block->getKernel();
            kernel->setCollisionFactor(kernel->getCollisionFactor()*0.5);
         }
      }
   }
}


