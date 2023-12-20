//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup cpu_Utilities Utilities
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef VoxelMatrixUtil_h__
#define VoxelMatrixUtil_h__

#include "GbCuboid3D.h"
#include "NoSlipBC.h"
#include "D3Q27Interactor.h"
#include "SetBcBlocksBlockVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"


namespace Utilities
{
   void voxelMatrixDiscretisation(SPtr<GbVoxelMatrix3D> matrix, std::string& pathname, int myid, int fileCounter, SPtr<Grid3D> grid, int bounceBackOption, bool vmFile)
   {
      SPtr<BC> noSlipPM(new NoSlipBC(bounceBackOption));
      SPtr<D3Q27Interactor> vmInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(matrix, grid, noSlipPM, Interactor3D::SOLID));

      if (vmFile)
      {
         if (myid == 0) matrix->writeToVTKImageDataASCII(pathname + "/geo/vmatrix" + UbSystem::toString(fileCounter));
      } 
      else
      {
         GbCuboid3DPtr vmBox(new GbCuboid3D(matrix->getX1Minimum(), matrix->getX2Minimum(), matrix->getX3Minimum(), matrix->getX1Maximum(), matrix->getX2Maximum(), matrix->getX3Maximum()));
         if (myid == 0) GbSystem3D::writeGeoObject(vmBox.get(), pathname + "/geo/vmbox" + UbSystem::toString(fileCounter), WbWriterVtkXmlASCII::getInstance());
      }

      //GbCuboid3DPtr vmBox(new GbCuboid3D(matrix->getX1Minimum(), matrix->getX2Minimum(), matrix->getX3Minimum(), matrix->getX1Maximum(), matrix->getX2Maximum(), matrix->getX3Maximum()));
      //if (myid == 0) GbSystem3D::writeGeoObject(vmBox.get(), pathname + "/geo/vmbox" + UbSystem::toString(fileCounter), WbWriterVtkXmlASCII::getInstance());
      //SPtr<D3Q27Interactor> vmBoxInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(vmBox, grid, noSlipPM, Interactor3D::SOLID));
      //SetSolidOrBoundaryBlockVisitor v1(vmBoxInt, SetSolidOrBoundaryBlockVisitor::SOLID);
      //grid->accept(v1);
      //SetSolidOrBoundaryBlockVisitor v2(vmBoxInt, SetSolidOrBoundaryBlockVisitor::BC);
      //grid->accept(v2);

      //std::vector<SPtr<Block3D>> blocks;
      //std::vector<SPtr<Block3D>>& sb = vmBoxInt->getSolidBlockSet();
      //if (myid == 0) UBLOG(logINFO, "number of solid blocks = " << sb.size());
      //blocks.insert(blocks.end(), sb.begin(), sb.end());
      //std::vector<SPtr<Block3D>>& tb = vmBoxInt->getBcBlocks();
      //if (myid == 0) UBLOG(logINFO, "number of trans blocks = " << tb.size());
      //blocks.insert(blocks.end(), tb.begin(), tb.end());

      //if (myid == 0) UBLOG(logINFO, "number of blocks = " << blocks.size());

      //for(SPtr<Block3D> block : blocks)
      //{
      //   block->setActive(true);
      //   vmInt->setDifferencesToGbObject3D(block);
      //}

      SetBcBlocksBlockVisitor v(vmInt);
      grid->accept(v);
      vmInt->initInteractor();
   }
}
#endif // VoxelMatrixUtil_h__


//! \}
