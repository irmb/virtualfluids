#ifndef VoxelMatrixUtil_h__
#define VoxelMatrixUtil_h__

#include "GbCuboid3D.h"
#include "NoSlipBCAdapter.h"
#include "D3Q27Interactor.h"
#include "SetSolidBlockVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"


namespace Utilities
{
   void voxelMatrixDiscretisation(SPtr<GbVoxelMatrix3D> matrix, std::string& pathname, int myid, int fileCounter, SPtr<Grid3D> grid, int bounceBackOption, bool vmFile)
   {
      SPtr<BCAdapter> noSlipPM(new NoSlipBCAdapter(bounceBackOption));
      SPtr<D3Q27Interactor> vmInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(matrix, grid, noSlipPM, Interactor3D::SOLID));

      if (vmFile)
      {
         if (myid == 0) matrix->writeToVTKImageDataASCII(pathname + "/geo/vmatrix" + UbSystem::toString(fileCounter));
      } 

      GbCuboid3DPtr vmBox(new GbCuboid3D(matrix->getX1Minimum(), matrix->getX2Minimum(), matrix->getX3Minimum(), matrix->getX1Maximum(), matrix->getX2Maximum(), matrix->getX3Maximum()));
      if (myid == 0) GbSystem3D::writeGeoObject(vmBox.get(), pathname + "/geo/vmbox" + UbSystem::toString(fileCounter), WbWriterVtkXmlASCII::getInstance());
      SPtr<D3Q27Interactor> vmBoxInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(vmBox, grid, noSlipPM, Interactor3D::SOLID));
      SetSolidBlockVisitor v1(vmBoxInt, BlockType::SOLID);
      grid->accept(v1);
      SetSolidBlockVisitor v2(vmBoxInt, BlockType::BC);
      grid->accept(v2);

      std::vector<SPtr<Block3D>> blocks;
      std::vector<SPtr<Block3D>>& sb = vmBoxInt->getSolidBlockSet();
      if (myid == 0) UBLOG(logINFO, "number of solid blocks = " << sb.size());
      blocks.insert(blocks.end(), sb.begin(), sb.end());
      std::vector<SPtr<Block3D>>& tb = vmBoxInt->getBcBlocks();
      if (myid == 0) UBLOG(logINFO, "number of trans blocks = " << tb.size());
      blocks.insert(blocks.end(), tb.begin(), tb.end());

      if (myid == 0) UBLOG(logINFO, "number of blocks = " << blocks.size());

      for(SPtr<Block3D> block : blocks)
      {
         block->setActive(true);
         vmInt->setDifferencesToGbObject3D(block);
      }
   }
}
#endif // VoxelMatrixUtil_h__

