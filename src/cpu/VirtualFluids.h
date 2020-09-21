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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file VirtualFluids.h
//! \ingroup Applications
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef VirtualFluids_h__
#define VirtualFluids_h__

//VirtualFluids header files

#ifdef _OPENMP
#include <omp.h>
#endif

#include <basics/PointerDefinitions.h>

#include <muParser.h>

#include <basics/container/CbArray2D.h>
#include <basics/container/CbArray3D.h>
#include <basics/container/CbArray4D.h>
#include <basics/container/CbVector.h>

#include <basics/objects/ObObject.h>

#include <basics/utilities/UbComparators.h>
#include <basics/utilities/UbEqual.h>
#include <basics/utilities/UbException.h>
#include <basics/utilities/UbInfinity.h>
#include <basics/utilities/UbKeys.h>
#include <basics/utilities/UbLimits.h>
#include <basics/utilities/UbLogger.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbObservable.h>
#include <basics/utilities/UbObserver.h>
#include <basics/utilities/UbScheduler.h>
#include <basics/utilities/UbSystem.h>
#include <basics/utilities/UbTiming.h>
#include <basics/utilities/UbTuple.h>

#include <basics/writer/WbWriter.h>
#include <basics/writer/WbWriterVtkXmlASCII.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include <BoundaryConditions/BCArray3D.h>
#include <BoundaryConditions/BCProcessor.h>
#include <BoundaryConditions/BCAlgorithm.h>
#include <BoundaryConditions/BCFunction.h>
#include <BoundaryConditions/BoundaryConditions.h>
#include <BoundaryConditions/BCAdapter.h>
#include <BoundaryConditions/BCProcessor.h>
#include <BoundaryConditions/NoSlipBCAdapter.h>
#include <BoundaryConditions/VelocityBCAdapter.h>
#include <BoundaryConditions/BCAlgorithm.h>
#include <BoundaryConditions/VelocityBCAlgorithm.h>
#include <BoundaryConditions/NoSlipBCAlgorithm.h>

#include <Connectors/Block3DConnector.h>
#include <Connectors/D3Q27ETFullDirectConnector.h>
#include <Connectors/LocalBlock3DConnector.h>

#include <Data/D3Q27EsoTwist3DSplittedVector.h>
#include <Data/DataSet3D.h>
#include <Data/DistributionArray3D.h>
#include <Data/EsoTwist3D.h>
#include <Data/EsoTwistD3Q27System.h>

#include <Grid/Block3D.h>
#include <Grid/Calculator.h>
#include <Grid/BasicCalculator.h>
#include <Grid/Grid3D.h>
#include <Grid/Grid3DSystem.h>

#include <Interactors/D3Q27Interactor.h>
#include <Interactors/Interactor3D.h>
#include <Interactors/InteractorsHelper.h>

#include <CoProcessors/WriteBlocksCoProcessor.h>
#include <CoProcessors/WriteMacroscopicQuantitiesCoProcessor.h>
#include <CoProcessors/WriteBoundaryConditionsCoProcessor.h>
#include <CoProcessors/NUPSCounterCoProcessor.h>
#include <CoProcessors/CoProcessor.h>

#include <LBM/D3Q27System.h>
#include <LBM/LBMKernel.h>
#include <LBM/ILBMKernel.h>
#include <LBM/CumulantK17LBMKernel.h>
#include <LBM/LBMSystem.h>
#include <LBM/LBMUnitConverter.h>

#include <geometry3d/CoordinateTransformation3D.h>
#include <geometry3d/GbCuboid3D.h>
#include <geometry3d/GbLine3D.h>
#include <geometry3d/GbObject3D.h>
#include <geometry3d/GbPoint3D.h>
#include <geometry3d/GbPolygon3D.h>
#include <geometry3d/GbSystem3D.h>
#include <geometry3d/GbTriangle3D.h>
#include <geometry3d/GbVector3D.h>

#include <Parallel/Communicator.h>
#include <Parallel/NullCommunicator.h>

#include <Utilities/MemoryUtil.h>

#include <Visitors/Block3DVisitor.h>
#include <Visitors/InitDistributionsBlockVisitor.h>
#include <Visitors/SetConnectorsBlockVisitor.h>
#include <Visitors/GenBlocksGridVisitor.h>
#include <Visitors/Grid3DVisitor.h>
#include <Visitors/SetKernelBlockVisitor.h>
#include <Visitors/BoundaryConditionsBlockVisitor.h>

#endif // VirtualFluids_h__
