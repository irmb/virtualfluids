#ifndef VoxelMatrixUtil_h__
#define VoxelMatrixUtil_h__

#include "GbCuboid3D.h"
#include "D3Q27NoSlipBCAdapter.h"
#include "D3Q27Interactor.h"
#include "SetSolidOrTransBlockVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"


namespace Utilities
{
   void voxelMatrixDiscretisation(GbVoxelMatrix3DPtr matrix, std::string& pathname, int myid, int fileCounter, Grid3DPtr grid, int bounceBackOption, bool vmFile)
   {
      D3Q27BoundaryConditionAdapterPtr noSlipPM(new D3Q27NoSlipBCAdapter(bounceBackOption));
      D3Q27InteractorPtr vmInt = D3Q27InteractorPtr(new D3Q27Interactor(matrix, grid, noSlipPM, Interactor3D::SOLID));

      if (vmFile)
      {
         if (myid == 0) matrix->writeToVTKImageDataASCII(pathname + "/geo/vmatrix" + UbSystem::toString(fileCounter));
      } 

      GbCuboid3DPtr vmBox(new GbCuboid3D(matrix->getX1Minimum(), matrix->getX2Minimum(), matrix->getX3Minimum(), matrix->getX1Maximum(), matrix->getX2Maximum(), matrix->getX3Maximum()));
      if (myid == 0) GbSystem3D::writeGeoObject(vmBox.get(), pathname + "/geo/vmbox" + UbSystem::toString(fileCounter), WbWriterVtkXmlASCII::getInstance());
      D3Q27InteractorPtr vmBoxInt = D3Q27InteractorPtr(new D3Q27Interactor(vmBox, grid, noSlipPM, Interactor3D::SOLID));
      SetSolidOrTransBlockVisitor v1(vmBoxInt, SetSolidOrTransBlockVisitor::SOLID);
      grid->accept(v1);
      SetSolidOrTransBlockVisitor v2(vmBoxInt, SetSolidOrTransBlockVisitor::TRANS);
      grid->accept(v2);

      std::vector<Block3DPtr> blocks;
      std::vector<Block3DPtr>& sb = vmBoxInt->getSolidBlockSet();
      if (myid == 0) UBLOG(logINFO, "number of solid blocks = " << sb.size());
      blocks.insert(blocks.end(), sb.begin(), sb.end());
      std::vector<Block3DPtr>& tb = vmBoxInt->getTransBlockSet();
      if (myid == 0) UBLOG(logINFO, "number of trans blocks = " << tb.size());
      blocks.insert(blocks.end(), tb.begin(), tb.end());

      if (myid == 0) UBLOG(logINFO, "number of blocks = " << blocks.size());

      BOOST_FOREACH(Block3DPtr block, blocks)
      {
         block->setActive(true);
         vmInt->setDifferencesToGbObject3D(block);
      }
   }
}
#endif // VoxelMatrixUtil_h__

