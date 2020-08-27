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

#include <basics/PointerDefinitions.h>

#include <basics/container/CbArray2D.h>
#include <basics/container/CbArray3D.h>
#include <basics/container/CbArray4D.h>
#include <basics/container/CbVector.h>
#include <basics/container/CbVectorPool.h>
#include <basics/memory/MbSmartPtr.h>
#include <basics/memory/MbSmartPtrBase.h>
#include <basics/objects/ObObject.h>
#include <basics/transmitter/TbTransmitter.h>
#include <basics/transmitter/TbTransmitterLocal.h>
#include <basics/transmitter/TbTransmitterMpiPool.h>
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
#include <basics/utilities/UbNupsTimer.h>
#include <basics/utilities/UbObservable.h>
#include <basics/utilities/UbObserver.h>
#include <basics/utilities/UbRandom.h>
#include <basics/utilities/UbScheduler.h>
#include <basics/utilities/UbStringInputASCII.h>
#include <basics/utilities/UbSystem.h>
#include <basics/utilities/UbTiming.h>
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

#include <BoundaryConditions/BCAdapter.h>
#include <BoundaryConditions/BCAlgorithm.h>
#include <BoundaryConditions/BCArray3D.h>
#include <BoundaryConditions/BCFunction.h>
#include <BoundaryConditions/BCProcessor.h>
#include <BoundaryConditions/BoundaryConditions.h>
#include <BoundaryConditions/DensityBCAdapter.h>
#include <BoundaryConditions/EqDensityBCAlgorithm.h>
#include <BoundaryConditions/HighViscosityNoSlipBCAlgorithm.h>
#include <BoundaryConditions/NoSlipBCAdapter.h>
#include <BoundaryConditions/NoSlipBCAlgorithm.h>
#include <BoundaryConditions/NonEqDensityBCAlgorithm.h>
#include <BoundaryConditions/NonReflectingOutflowBCAlgorithm.h>
#include <BoundaryConditions/SlipBCAdapter.h>
#include <BoundaryConditions/SlipBCAlgorithm.h>
#include <BoundaryConditions/ThinWallBCProcessor.h>
#include <BoundaryConditions/ThinWallNoSlipBCAlgorithm.h>
#include <BoundaryConditions/VelocityBCAdapter.h>
#include <BoundaryConditions/VelocityBCAlgorithm.h>
#include <BoundaryConditions/VelocityWithDensityBCAlgorithm.h>
#include <BoundaryConditions/DensityAndThixotropyBCAlgorithm.h>
#include <BoundaryConditions/NoSlipAndThixotropyBCAlgorithm.h>
#include <BoundaryConditions/VelocityAndThixotropyBCAlgorithm.h>
#include <BoundaryConditions/NonReflectingOutflowAndThixotropyBCAlgorithm.h>
#include <BoundaryConditions/VelocityWithDensityAndThixotropyBCAlgorithm.h>
#include <BoundaryConditions/SimpleVelocityBCAlgorithm.h>
#include <BoundaryConditions/ThixotropyNoSlipBCAlgorithm.h>
#include <BoundaryConditions/BinghamModelNoSlipBCAlgorithm.h>
#include <BoundaryConditions/HerschelBulkleyModelNoSlipBCAlgorithm.h>

#include <Connectors/Block3DConnector.h>
#include <Connectors/Block3DConnectorFactory.h>
#include <Connectors/CoarseToFineBlock3DConnector.h>
#include <Connectors/CoarseToFineNodeSetBlock3DConnector.h>
#include <Connectors/ConnectorFactory.h>
#include <Connectors/D3Q27ETCFOffVectorConnector.h>
#include <Connectors/D3Q27ETFCOffVectorConnector.h>
#include <Connectors/D3Q27ETFullDirectConnector.h>
#include <Connectors/D3Q27ETFullVectorConnector.h>
#include <Connectors/FineToCoarseBlock3DConnector.h>
#include <Connectors/FineToCoarseNodeSetBlock3DConnector.h>
#include <Connectors/LocalBlock3DConnector.h>
#include <Connectors/RemoteBlock3DConnector.h>

#include <Data/D3Q27EsoTwist3DSplittedVector.h>
#include <Data/D3Q27EsoTwist3DSplittedVectorEx.h>
#include <Data/DataSet3D.h>
#include <Data/DistributionArray3D.h>
#include <Data/EsoTwist3D.h>
#include <Data/EsoTwistD3Q27System.h>
#include <Data/VoidData3D.h>

#include <Grid/BasicCalculator.h>
#include <Grid/Block3D.h>
#include <Grid/Calculator.h>
#include <Grid/Grid3D.h>
#include <Grid/Grid3DSystem.h>

#include <Interactors/D3Q27Interactor.h>
#include <Interactors/D3Q27TriFaceMeshInteractor.h>
#include <Interactors/Interactor3D.h>
#include <Interactors/InteractorsHelper.h>

#include <CoProcessors/AdjustForcingCoProcessor.h>
#include <CoProcessors/CalculateForcesCoProcessor.h>
#include <CoProcessors/WriteBlocksCoProcessor.h>
#include <CoProcessors/WriteBoundaryConditionsCoProcessor.h>
#include <CoProcessors/WriteMQFromSelectionCoProcessor.h>
#include <CoProcessors/WriteMacroscopicQuantitiesCoProcessor.h>
//#include <CoProcessors/PathLineCoProcessor.h>
//#include <CoProcessors/PathLineCoProcessorMcpart.h>
#include <CoProcessors/EmergencyExitCoProcessor.h>
#include <CoProcessors/NUPSCounterCoProcessor.h>
#include <CoProcessors/PressureDifferenceCoProcessor.h>
//#include <CoProcessors/Particles.h>
#include <CoProcessors/AverageValuesCoProcessor.h>
#include <CoProcessors/CoProcessor.h>
#include <CoProcessors/DecreaseViscosityCoProcessor.h>
#include <CoProcessors/InSituVTKCoProcessor.h>
#include <CoProcessors/QCriterionCoProcessor.h>
#include <CoProcessors/ShearStressCoProcessor.h>
#include <CoProcessors/TimeseriesCoProcessor.h>
#include <CoProcessors/TurbulenceIntensityCoProcessor.h>
//#include <CoProcessors/MeanValuesCoProcessor.h>
#include <CoProcessors/InSituCatalystCoProcessor.h>
#include <CoProcessors/LineTimeSeriesCoProcessor.h>
#include <CoProcessors/MPIIOMigrationBECoProcessor.h>
#include <CoProcessors/MPIIOMigrationCoProcessor.h>
#include <CoProcessors/MPIIORestartCoProcessor.h>
#include <CoProcessors/MicrophoneArrayCoProcessor.h>
#include <WriteThixotropyQuantitiesCoProcessor.h>

#include <IntegrateValuesHelper.h>
//#include <LBM/D3Q27CompactInterpolationProcessor.h>
#include <LBM/CompressibleOffsetInterpolationProcessor.h>
#include <LBM/CompressibleOffsetMomentsInterpolationProcessor.h>
#include <LBM/CompressibleOffsetSquarePressureInterpolationProcessor.h>
#include <LBM/IncompressibleOffsetInterpolationProcessor.h>
#include <LBM/InterpolationHelper.h>
#include <LBM/InterpolationProcessor.h>
//#include <LBM/D3Q27OffsetInterpolationProcessor.h>
#include <IncompressibleCumulantWithSpongeLayerLBMKernel.h>
#include <LBM/CompressibleCumulant4thOrderViscosityLBMKernel.h>
#include <LBM/CompressibleCumulantLBMKernel.h>
#include <LBM/D3Q27System.h>
#include <LBM/ICell.h>
#include <LBM/IncompressibleCumulantLBMKernel.h>
#include <LBM/InitDensityLBMKernel.h>
#include <LBM/InterpolationProcessor.h>
#include <LBM/LBMKernel.h>
#include <LBM/LBMKernelETD3Q27BGK.h>
#include <LBM/LBMSystem.h>
#include <LBM/LBMUnitConverter.h>
//#include <LBM/BGKLBMKernel.h>
#include <LBM/ThixotropyLBMKernel.h>
#include <LBM/ThixotropyExpLBMKernel.h>
#include <LBM/CumulantLBMKernel.h>
#include <LBM/ThixotropyModelLBMKernel.h>
#include <LBM/BinghamModelLBMKernel.h>
#include <LBM/HerschelBulkleyModelLBMKernel.h>


#include <geometry3d/CoordinateTransformation3D.h>
#include <geometry3d/GbCuboid3D.h>
#include <geometry3d/GbCylinder3D.h>
#include <geometry3d/GbHalfSpace3D.h>
#include <geometry3d/GbHalfSpaceKrischan3D.h>
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

#include <Parallel/BlocksDistributor.h>
#include <Parallel/Communicator.h>
#include <Parallel/MPICommunicator.h>
#include <Parallel/MetisPartitioner.h>
#include <Parallel/NullCommunicator.h>
#include <Parallel/PriorityQueueDecompositor.h>
#include <Parallel/SimpleGeometricPartitioner.h>
#include <Parallel/ZoltanPartitioner.h>

#include <Utilities/ChangeRandomQs.hpp>
#include <Utilities/CheckpointConverter.h>
#include <Utilities/ConfigurationFile.hpp>
#include <Utilities/MathUtil.hpp>
#include <Utilities/MemoryUtil.h>
#include <Utilities/VoxelMatrixUtil.hpp>

#include <CheckRatioBlockVisitor.h>
#include <InitDistributionsFromFileBlockVisitor.h>
#include <InitDistributionsWithInterpolationGridVisitor.h>
#include <SpongeLayerBlockVisitor.h>
#include <Visitors/Block3DVisitor.h>
#include <Visitors/BoundaryConditionsBlockVisitor.h>
#include <Visitors/ChangeBoundaryDensityBlockVisitor.h>
#include <Visitors/CoarsenCrossAndInsideGbObjectBlockVisitor.h>
#include <Visitors/ConnectorBlockVisitor.h>
#include <Visitors/CreateTransmittersHelper.h>
#include <Visitors/GenBlocksGridVisitor.h>
#include <Visitors/Grid3DVisitor.h>
#include <Visitors/InitDistributionsBlockVisitor.h>
#include <Visitors/MetisPartitioningGridVisitor.h>
#include <Visitors/OverlapBlockVisitor.h>
#include <Visitors/PQueuePartitioningGridVisitor.h>
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
#include <Visitors/SetSpongeLayerBlockVisitor.h>
#include <Visitors/SetUndefinedNodesBlockVisitor.h>
#include <Visitors/ViscosityBlockVisitor.h>
#include <Visitors/BoundaryConditionsBlockVisitor.h>
#include <Visitors/BoundaryConditionsBlockVisitor.h>
#include <Visitors/ChangeBoundaryDensityBlockVisitor.h>
#include <InitDistributionsFromFileBlockVisitor.h>
#include <InitDistributionsWithInterpolationGridVisitor.h>
#include <InitThixotropyBlockVisitor.h>
#include <CheckRatioBlockVisitor.h>
#include <SpongeLayerBlockVisitor.h>
#include <ZoltanPartitioningGridVisitor.h>

#include <RefineAroundGbObjectHelper.h>
#include <Visitors/RefineCrossAndInsideGbObjectHelper.h>

#endif // VirtualFluids_h__
