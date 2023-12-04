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

// VirtualFluids header files

#ifdef _OPENMP
#include <omp.h>
#endif

#include <parallel/Communicator.h>
#include <parallel/MPICommunicator.h>
#include <parallel/NullCommunicator.h>
#include <parallel/transmitter/TbTransmitter.h>
#include <parallel/transmitter/TbTransmitterLocal.h>
#include <parallel/transmitter/TbTransmitterMpiPool.h>

#include <basics/PointerDefinitions.h>

#include <basics/config/ConfigurationFile.h>
#include <logger/Logger.h>

#include <basics/container/CbArray2D.h>
#include <basics/container/CbArray3D.h>
#include <basics/container/CbArray4D.h>
#include <basics/container/CbVector.h>
#include <basics/container/CbVectorPool.h>
#include <basics/memory/MbSmartPtr.h>
#include <basics/memory/MbSmartPtrBase.h>
#include <basics/objects/ObObject.h>
#include <basics/utilities/UbComparators.h>
#include <basics/utilities/UbEqual.h>
#include <basics/utilities/UbException.h>
#include <basics/utilities/UbFileInput.h>
#include <basics/utilities/UbFileInputASCII.h>
#include <basics/utilities/UbFileInputBinary.h>
#include <basics/utilities/UbFileOutput.h>
#include <basics/utilities/UbFileOutputASCII.h>
#include <basics/utilities/UbFileOutputBinary.h>
#include <basics/utilities/UbInfinity.h>
#include <basics/utilities/UbKeys.h>
#include <basics/utilities/UbLimits.h>
#include <basics/utilities/UbLogger.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbObservable.h>
#include <basics/utilities/UbObserver.h>
#include <basics/utilities/UbRandom.h>
#include <basics/utilities/UbScheduler.h>
#include <basics/utilities/UbStringInputASCII.h>
#include <basics/utilities/UbSystem.h>
#include <basics/utilities/UbTuple.h>
#include <basics/writer/WbWriter.h>
#include <basics/writer/WbWriterAvsASCII.h>
#include <basics/writer/WbWriterAvsBinary.h>
#include <basics/writer/WbWriterBOBJ.h>
#include <basics/writer/WbWriterSunflow.h>
#include <basics/writer/WbWriterTecPlotASCII.h>
#include <basics/writer/WbWriterVtkASCII.h>
#include <basics/writer/WbWriterVtkBinary.h>
#include <basics/writer/WbWriterVtkXmlASCII.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>
#include <basics/writer/WbWriterX3D.h>
#include <muParser.h>
#include <muParserBase.h>
#include <muParserBytecode.h>
#include <muParserCallback.h>
#include <muParserDLL.h>
#include <muParserDef.h>
#include <muParserError.h>
#include <muParserFixes.h>
#include <muParserInt.h>
#include <muParserTemplateMagic.h>
#include <muParserTest.h>
#include <muParserToken.h>
#include <muParserTokenReader.h>

#include <cpu/core/BoundaryConditions/BC.h>
#include <cpu/core/BoundaryConditions/BCArray3D.h>
#include <cpu/core/BoundaryConditions/BCFunction.h>
#include <cpu/core/BoundaryConditions/BCSet.h>
#include <cpu/core/BoundaryConditions/BCStrategy.h>
#include <cpu/core/BoundaryConditions/BoundaryConditions.h>
#include <cpu/core/BoundaryConditions/NoSlipBC.h>
#include <cpu/core/BoundaryConditions/NoSlipInterpolated.h>
#include <cpu/core/BoundaryConditions/OutflowNonReflecting.h>
#include <cpu/core/BoundaryConditions/OutflowNonReflectingWithPressure.h>
#include <cpu/core/BoundaryConditions/PressureBC.h>
#include <cpu/core/BoundaryConditions/PressureNonEquilibrium.h>
#include <cpu/core/BoundaryConditions/SlipBC.h>
#include <cpu/core/BoundaryConditions/SlipBounceBack.h>
#include <cpu/core/BoundaryConditions/SlipInterpolated.h>
#include <cpu/core/BoundaryConditions/ThinWallBCSet.h>
#include <cpu/core/BoundaryConditions/ThinWallNoSlip.h>
#include <cpu/core/BoundaryConditions/VelocityBC.h>
#include <cpu/core/BoundaryConditions/VelocityBounceBack.h>
#include <cpu/core/BoundaryConditions/VelocityInterpolated.h>
#include <cpu/core/BoundaryConditions/VelocityNonReflecting.h>
#include <cpu/core/BoundaryConditions/VelocityWithPressureInterpolated.h>

#include <cpu/core/Connectors/Block3DConnector.h>
#include <cpu/core/Connectors/CoarseToFineVectorConnector.h>
#include <cpu/core/Connectors/FineToCoarseVectorConnector.h>
#include <cpu/core/Connectors/LocalBlock3DConnector.h>
#include <cpu/core/Connectors/OneDistributionFullDirectConnector.h>
#include <cpu/core/Connectors/OneDistributionFullVectorConnector.h>
#include <cpu/core/Connectors/RemoteBlock3DConnector.h>
#include <cpu/core/Connectors/TwoDistributionsFullDirectConnector.h>
#include <cpu/core/Connectors/TwoDistributionsFullVectorConnector.h>

#include <cpu/core/Data/DataSet3D.h>
#include <cpu/core/Data/DistributionArray3D.h>
#include <cpu/core/Data/EsoSplit.h>
#include <cpu/core/Data/EsoTwist3D.h>

#include <cpu/core/Simulation/Block3D.h>
#include <cpu/core/Simulation/Grid3D.h>
#include <cpu/core/Simulation/Simulation.h>

#include <cpu/core/Interactors/D3Q27Interactor.h>
#include <cpu/core/Interactors/D3Q27TriFaceMeshInteractor.h>
#include <cpu/core/Interactors/Interactor3D.h>
#include <cpu/core/Interactors/InteractorsHelper.h>

#include <cpu/core/SimulationObservers/AdjustForcingSimulationObserver.h>
#include <cpu/core/SimulationObservers/AverageValuesSimulationObserver.h>
#include <cpu/core/SimulationObservers/CalculateForcesSimulationObserver.h>
#include <cpu/core/SimulationObservers/DecreaseViscositySimulationObserver.h>
#include <cpu/core/SimulationObservers/EmergencyExitSimulationObserver.h>
#include <cpu/core/SimulationObservers/InSituCatalystSimulationObserver.h>
#include <cpu/core/SimulationObservers/InSituVTKSimulationObserver.h>
#include <cpu/core/SimulationObservers/IntegrateValuesHelper.h>
#include <cpu/core/SimulationObservers/LineTimeSeriesSimulationObserver.h>
#include <cpu/core/SimulationObservers/MPIIOMigrationBESimulationObserver.h>
#include <cpu/core/SimulationObservers/MPIIOMigrationSimulationObserver.h>
#include <cpu/core/SimulationObservers/MPIIORestartSimulationObserver.h>
#include <cpu/core/SimulationObservers/MicrophoneArraySimulationObserver.h>
#include <cpu/core/SimulationObservers/NUPSCounterSimulationObserver.h>
#include <cpu/core/SimulationObservers/PressureDifferenceSimulationObserver.h>
#include <cpu/core/SimulationObservers/QCriterionSimulationObserver.h>
#include <cpu/core/SimulationObservers/ShearStressSimulationObserver.h>
#include <cpu/core/SimulationObservers/SimulationObserver.h>
#include <cpu/core/SimulationObservers/TimeAveragedValuesSimulationObserver.h>
#include <cpu/core/SimulationObservers/TimeDependentBCSimulationObserver.h>
#include <cpu/core/SimulationObservers/TimeseriesSimulationObserver.h>
#include <cpu/core/SimulationObservers/TurbulenceIntensitySimulationObserver.h>
#include <cpu/core/SimulationObservers/WriteBlocksSimulationObserver.h>
#include <cpu/core/SimulationObservers/WriteBoundaryConditionsSimulationObserver.h>
#include <cpu/core/SimulationObservers/WriteMQFromSelectionSimulationObserver.h>
#include <cpu/core/SimulationObservers/WriteMacroscopicQuantitiesSimulationObserver.h>

#include <cpu/core/LBM/D3Q27System.h>
#include <cpu/core/LBM/Interpolation/CompressibleOffsetMomentsInterpolator.h>
#include <cpu/core/LBM/Interpolation/ICell.h>
#include <cpu/core/LBM/Interpolation/IncompressibleOffsetInterpolator.h>
#include <cpu/core/LBM/Interpolation/Interpolator.h>
#include <cpu/core/LBM/LBMKernel.h>
#include <cpu/core/LBM/LBMSystem.h>
#include <cpu/core/LBM/LBMUnitConverter.h>

#include <cpu/core/LBM/B92IncompressibleNavierStokes.h>
#include <cpu/core/LBM/K15CompressibleNavierStokes.h>
#include <cpu/core/LBM/K16IncompressibleNavierStokes.h>
#include <cpu/core/LBM/K17CompressibleNavierStokes.h>

#include <geometry3d/CoordinateTransformation3D.h>
#include <geometry3d/GbCuboid3D.h>
#include <geometry3d/GbCylinder3D.h>
#include <geometry3d/GbHalfSpace3D.h>
#include <geometry3d/GbHalfSpaceKrischan3D.h>
#include <geometry3d/GbGyroidThirdOrder.h>
#include <geometry3d/GbGyroidThirdOrderLong.h>
#include <geometry3d/GbLine3D.h>
#include <geometry3d/GbMeshTools3D.h>
#include <geometry3d/GbObject3D.h>
#include <geometry3d/GbObjectGroup3D.h>
#include <geometry3d/GbPoint3D.h>
#include <geometry3d/GbPolygon3D.h>
#include <geometry3d/GbQuadFaceMesh3D.h>
#include <geometry3d/GbSphere3D.h>
#include <geometry3d/GbSystem3D.h>
#include <geometry3d/GbTriFaceMesh3D.h>
#include <geometry3d/GbTriangle3D.h>
#include <geometry3d/GbTriangularMesh3D.h>
#include <geometry3d/GbVector3D.h>
#include <geometry3d/GbVoxelMatrix3D.h>
#include <geometry3d/KdTree/KdNode.h>
#include <geometry3d/KdTree/KdRay.h>
#include <geometry3d/KdTree/KdSplitCandidate.h>
#include <geometry3d/KdTree/KdSplitCandidateManager.h>
#include <geometry3d/KdTree/KdTree.h>
#include <geometry3d/KdTree/KdUtilities.h>
#include <geometry3d/KdTree/intersectionhandler/KdCountLineIntersectionHandler.h>
#include <geometry3d/KdTree/intersectionhandler/KdCountRayIntersectionHandler.h>
#include <geometry3d/KdTree/intersectionhandler/KdLineIntersectionHandler.h>
#include <geometry3d/KdTree/intersectionhandler/KdRayIntersectionHandler.h>
#include <geometry3d/KdTree/splitalgorithms/KdSAHSplit.h>
#include <geometry3d/KdTree/splitalgorithms/KdSpatiallMedianSplit.h>
#include <geometry3d/KdTree/splitalgorithms/KdSplitAlgorithm.h>

#include <cpu/core/Parallel/MetisPartitioner.h>

#include <cpu/core/Utilities/ChangeRandomQs.hpp>
#include <cpu/core/Utilities/CheckpointConverter.h>
#include <cpu/core/Utilities/MathUtil.hpp>
#include <cpu/core/Utilities/MemoryUtil.h>
#include <cpu/core/Utilities/VoxelMatrixUtil.hpp>

#include <cpu/core/Visitors/Block3DVisitor.h>
#include <cpu/core/Visitors/BoundaryConditionsBlockVisitor.h>
#include <cpu/core/Visitors/ChangeBoundaryDensityBlockVisitor.h>
#include <cpu/core/Visitors/CheckRatioBlockVisitor.h>
#include <cpu/core/Visitors/CoarsenCrossAndInsideGbObjectBlockVisitor.h>
#include <cpu/core/Visitors/CreateTransmittersHelper.h>
#include <cpu/core/Visitors/GenBlocksGridVisitor.h>
#include <cpu/core/Visitors/Grid3DVisitor.h>
#include <cpu/core/Visitors/InitDistributionsBlockVisitor.h>
#include <cpu/core/Visitors/InitDistributionsWithInterpolationGridVisitor.h>
#include <cpu/core/Visitors/MetisPartitioningGridVisitor.h>
#include <cpu/core/Visitors/OverlapBlockVisitor.h>
#include <cpu/core/Visitors/RatioBlockVisitor.h>
#include <cpu/core/Visitors/RatioSmoothBlockVisitor.h>
#include <cpu/core/Visitors/RefineAroundGbObjectHelper.h>
#include <cpu/core/Visitors/RefineCrossAndInsideGbObjectBlockVisitor.h>
#include <cpu/core/Visitors/RefineCrossAndInsideGbObjectHelper.h>
#include <cpu/core/Visitors/RefineInterGbObjectsVisitor.h>
#include <cpu/core/Visitors/RenumberBlockVisitor.h>
#include <cpu/core/Visitors/SetBcBlocksBlockVisitor.h>
#include <cpu/core/Visitors/SetConnectorsBlockVisitor.h>
#include <cpu/core/Visitors/SetForcingBlockVisitor.h>
#include <cpu/core/Visitors/SetInterpolationConnectorsBlockVisitor.h>
#include <cpu/core/Visitors/SetInterpolationDirsBlockVisitor.h>
#include <cpu/core/Visitors/SetKernelBlockVisitor.h>
#include <cpu/core/Visitors/SetSolidBlocksBlockVisitor.h>
#include <cpu/core/Visitors/SetUndefinedNodesBlockVisitor.h>
#include <cpu/core/Visitors/SpongeLayerBlockVisitor.h>
#include <cpu/core/Visitors/ViscosityBlockVisitor.h>

#endif // VirtualFluids_h__
