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

#include <BoundaryConditions/BC.h>
#include <BoundaryConditions/BCStrategy.h>
#include <BoundaryConditions/BCArray3D.h>
#include <BoundaryConditions/BCFunction.h>
#include <BoundaryConditions/BCSet.h>
#include <BoundaryConditions/BoundaryConditions.h>
#include <BoundaryConditions/PressureBC.h>
#include <BoundaryConditions/NoSlipBC.h>
#include <BoundaryConditions/NoSlipInterpolated.h>
#include <BoundaryConditions/PressureNonEquilibrium.h>
#include <BoundaryConditions/OutflowNonReflecting.h>
#include <BoundaryConditions/OutflowNonReflectingWithPressure.h>
#include <BoundaryConditions/VelocityNonReflecting.h>
#include <BoundaryConditions/SlipBC.h>
#include <BoundaryConditions/SlipInterpolated.h>
#include <BoundaryConditions/ThinWallBCSet.h>
#include <BoundaryConditions/ThinWallNoSlip.h>
#include <BoundaryConditions/VelocityBC.h>
#include <BoundaryConditions/VelocityInterpolated.h>
#include <BoundaryConditions/VelocityWithPressureInterpolated.h>
#include <BoundaryConditions/VelocityBounceBack.h>
#include <BoundaryConditions/SlipBounceBack.h>





#include <Connectors/Block3DConnector.h>
//#include <Connectors/Block3DConnectorFactory.h>
//#include <Connectors/CoarseToFineBlock3DConnector.h>
//#include <Connectors/CoarseToFineNodeSetBlock3DConnector.h>
//#include <Connectors/ConnectorFactory.h>
#include <Connectors/CoarseToFineVectorConnector.h>
#include <Connectors/FineToCoarseVectorConnector.h>
#include <Connectors/OneDistributionFullDirectConnector.h>
#include <Connectors/OneDistributionFullVectorConnector.h>
//#include <Connectors/FineToCoarseBlock3DConnector.h>
//#include <Connectors/FineToCoarseNodeSetBlock3DConnector.h>
#include <Connectors/LocalBlock3DConnector.h>
#include <Connectors/RemoteBlock3DConnector.h>
#include <Connectors/TwoDistributionsFullDirectConnector.h>
#include <Connectors/TwoDistributionsFullVectorConnector.h>


#include <Data/EsoSplit.h>

#include <Data/DataSet3D.h>
#include <Data/DistributionArray3D.h>
#include <Data/EsoTwist3D.h>


#include <Simulation/Block3D.h>
#include <Simulation/Simulation.h>
#include <Simulation/Grid3D.h>

#include <Interactors/D3Q27Interactor.h>
#include <Interactors/D3Q27TriFaceMeshInteractor.h>
#include <Interactors/Interactor3D.h>
#include <Interactors/InteractorsHelper.h>

#include <SimulationObservers/AdjustForcingSimulationObserver.h>
#include <SimulationObservers/CalculateForcesSimulationObserver.h>

#include <SimulationObservers/WriteMacroscopicQuantitiesSimulationObserver.h>
#include <SimulationObservers/WriteMQFromSelectionSimulationObserver.h>
#include <SimulationObservers/WriteBoundaryConditionsSimulationObserver.h>
#include <SimulationObservers/WriteMQFromSelectionSimulationObserver.h>
#include <SimulationObservers/WriteMacroscopicQuantitiesSimulationObserver.h>
#include <WriteBlocksSimulationObserver.h>
//#include <SimulationObservers/PathLineSimulationObserver.h>
//#include <SimulationObservers/PathLineSimulationObserverMcpart.h>
#include <SimulationObservers/EmergencyExitSimulationObserver.h>
#include <SimulationObservers/NUPSCounterSimulationObserver.h>
#include <SimulationObservers/PressureDifferenceSimulationObserver.h>
//#include <SimulationObservers/Particles.h>
#include <SimulationObservers/AverageValuesSimulationObserver.h>
#include <SimulationObservers/SimulationObserver.h>
#include <SimulationObservers/DecreaseViscositySimulationObserver.h>
#include <SimulationObservers/InSituVTKSimulationObserver.h>
#include <SimulationObservers/QCriterionSimulationObserver.h>
#include <SimulationObservers/ShearStressSimulationObserver.h>
#include <SimulationObservers/TimeseriesSimulationObserver.h>
#include <SimulationObservers/TurbulenceIntensitySimulationObserver.h>
#include <SimulationObservers/TimeAveragedValuesSimulationObserver.h>

//#include <SimulationObservers/MeanValuesSimulationObserver.h>
#include <SimulationObservers/InSituCatalystSimulationObserver.h>
#include <SimulationObservers/LineTimeSeriesSimulationObserver.h>
#include <SimulationObservers/MPIIOMigrationBESimulationObserver.h>
#include <SimulationObservers/MPIIOMigrationSimulationObserver.h>
#include <SimulationObservers/MPIIORestartSimulationObserver.h>
#include <SimulationObservers/MicrophoneArraySimulationObserver.h>


#include <TimeDependentBCSimulationObserver.h>

#include <IntegrateValuesHelper.h>
#include <LBM/Interpolation/CompressibleOffsetMomentsInterpolator.h>
#include <LBM/Interpolation/IncompressibleOffsetInterpolator.h>
#include <LBM/Interpolation/Interpolator.h>
#include <LBM/D3Q27System.h>
#include <LBM/Interpolation/ICell.h>
#include <LBM/LBMKernel.h>
#include <LBM/LBMSystem.h>
#include <LBM/LBMUnitConverter.h>


#include <LBM/B92IncompressibleNavierStokes.h>
#include <LBM/K15CompressibleNavierStokes.h>
#include <LBM/K16IncompressibleNavierStokes.h>
#include <LBM/K17CompressibleNavierStokes.h>


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

#include <Parallel/MetisPartitioner.h>

#include <Utilities/ChangeRandomQs.hpp>
#include <Utilities/CheckpointConverter.h>
#include <Utilities/MathUtil.hpp>
#include <Utilities/MemoryUtil.h>
#include <Utilities/VoxelMatrixUtil.hpp>

#include <CheckRatioBlockVisitor.h>
#include <InitDistributionsWithInterpolationGridVisitor.h>
#include <SpongeLayerBlockVisitor.h>
#include <Visitors/Block3DVisitor.h>
#include <Visitors/BoundaryConditionsBlockVisitor.h>
#include <Visitors/ChangeBoundaryDensityBlockVisitor.h>
#include <Visitors/CoarsenCrossAndInsideGbObjectBlockVisitor.h>
//#include <Visitors/ConnectorBlockVisitor.h>
#include <Visitors/CreateTransmittersHelper.h>
#include <Visitors/GenBlocksGridVisitor.h>
#include <Visitors/Grid3DVisitor.h>
#include <Visitors/InitDistributionsBlockVisitor.h>
#include <Visitors/MetisPartitioningGridVisitor.h>
#include <Visitors/OverlapBlockVisitor.h>
#include <Visitors/RatioBlockVisitor.h>
#include <Visitors/RatioSmoothBlockVisitor.h>
#include <Visitors/RefineCrossAndInsideGbObjectBlockVisitor.h>
#include <Visitors/RefineInterGbObjectsVisitor.h>
#include <Visitors/RenumberBlockVisitor.h>
#include <Visitors/SetBcBlocksBlockVisitor.h>
#include <Visitors/SetConnectorsBlockVisitor.h>
#include <Visitors/SetForcingBlockVisitor.h>
#include <Visitors/SetInterpolationDirsBlockVisitor.h>
#include <Visitors/SetKernelBlockVisitor.h>
#include <Visitors/SetSolidBlocksBlockVisitor.h>
#include <Visitors/SetUndefinedNodesBlockVisitor.h>
#include <Visitors/ViscosityBlockVisitor.h>
#include <Visitors/BoundaryConditionsBlockVisitor.h>
#include <Visitors/BoundaryConditionsBlockVisitor.h>
#include <Visitors/ChangeBoundaryDensityBlockVisitor.h>
#include <InitDistributionsWithInterpolationGridVisitor.h>

#include <CheckRatioBlockVisitor.h>
#include <SpongeLayerBlockVisitor.h>

#include <Visitors/SetInterpolationConnectorsBlockVisitor.h>

#include <RefineAroundGbObjectHelper.h>
#include <Visitors/RefineCrossAndInsideGbObjectHelper.h>

#endif // VirtualFluids_h__
