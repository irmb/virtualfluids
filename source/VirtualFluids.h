//VirtualFluids header files
 
#if defined VF_FETOL
#define WIN32_LEAN_AND_MEAN
#include <JM.h>
#endif

 
#include <MuParser/include/muParser.h>
#include <MuParser/include/muParserBase.h>
#include <MuParser/include/muParserBytecode.h>
#include <MuParser/include/muParserCallback.h>
#include <MuParser/include/muParserDef.h>
#include <MuParser/include/muParserDLL.h>
#include <MuParser/include/muParserError.h>
#include <MuParser/include/muParserFixes.h>
#include <MuParser/include/muParserInt.h>
#include <MuParser/include/muParserStack.h>
#include <MuParser/include/muParserTemplateMagic.h>
#include <MuParser/include/muParserTest.h>
#include <MuParser/include/muParserToken.h>
#include <MuParser/include/muParserTokenReader.h>
#include <basics/container/CbArray2D.h>
#include <basics/container/CbArray3D.h>
#include <basics/container/CbArray4D.h>
#include <basics/container/CbVector.h>
#include <basics/container/CbVectorPool.h>
#include <basics/memory/MbMemPool.h>
#include <basics/memory/MbSharedPointerDefines.h>
#include <basics/memory/MbSmartPtr.h>
#include <basics/memory/MbSmartPtrBase.h>
#include <basics/objects/ObCreator.h>
#include <basics/objects/ObFactory.h>
#include <basics/objects/ObObject.h>
#include <basics/objects/ObObjectCreator.h>
#include <basics/objects/ObObjectFactory.h>
#include <basics/objects/ObObjectManager.h>
#include <basics/transmitter/TbTransmitter.h>
#include <basics/transmitter/TbTransmitterLocal.h>
#include <basics/transmitter/TbTransmitterMpiPool.h>
#include <basics/utilities/UbAutoRun.hpp>
#include <basics/utilities/UbComparators.h>
#include <basics/utilities/UbConverter.h>
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
#include <basics/utilities/UbPointerWrapper.h>
#include <basics/utilities/UbRandom.h>
#include <basics/utilities/UbScheduler.h>
#include <basics/utilities/UbStaticPathMap.h>
#include <basics/utilities/UbString.h>
#include <basics/utilities/UbStringInputASCII.h>
#include <basics/utilities/UbSystem.h>
#include <basics/utilities/UbTableModel.h>
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
// #include <BoundaryCondition/BCArray.h>
#include <BoundaryConditions/BCArray3D.h>
#include <BoundaryConditions/BCProcessor.h>
#include <BoundaryConditions/BCAlgorithm.h>
#include <BoundaryConditions/BCFunction.h>
#include <BoundaryConditions/BoundaryConditions.h>
#include <BoundaryConditions/BCAdapter.h>
#include <BoundaryConditions/DensityBCAdapter.h>
#include <BoundaryConditions/BCProcessor.h>
#include <BoundaryConditions/ThinWallBCProcessor.h>
#include <BoundaryConditions/NoSlipBCAdapter.h>
#include <BoundaryConditions/SlipBCAdapter.h>
#include <BoundaryConditions/VelocityBCAdapter.h>

//#include <BoundaryCondition/BoundaryConditionProcessor.h>
#include <BoundaryConditions/BCAlgorithm.h>
#include <BoundaryConditions/VelocityBCAlgorithm.h>
#include <BoundaryConditions/NonEqDensityBCAlgorithm.h>
#include <BoundaryConditions/EqDensityBCAlgorithm.h>
#include <BoundaryConditions/NoSlipBCAlgorithm.h>
#include <BoundaryConditions/ThinWallNoSlipBCAlgorithm.h>
#include <BoundaryConditions/HighViscosityNoSlipBCAlgorithm.h>
#include <BoundaryConditions/SlipBCAlgorithm.h>
#include <BoundaryConditions/NonReflectingDensityBCAlgorithm.h>
#include <BoundaryConditions/NonReflectingVelocityBCAlgorithm.h>
#include <BoundaryConditions/NonReflectingSlipBCAlgorithm.h>

#include <Connectors/Block3DConnector.h>
#include <Connectors/D3Q27ETCFOffVectorConnector.h>
#include <Connectors/D3Q27ETFCOffVectorConnector.h>
#include <Connectors/D3Q27ETFullDirectConnector.h>
#include <Connectors/D3Q27ETFullVectorConnector.h>
#include <Connectors/LocalBlock3DConnector.h>
#include <Connectors/RemoteBlock3DConnector.h>
#include <Connectors/CoarseToFineBlock3DConnector.h>
#include <Connectors/CoarseToFineNodeSetBlock3DConnector.h>
#include <Connectors/FineToCoarseBlock3DConnector.h>
#include <Connectors/FineToCoarseNodeSetBlock3DConnector.h>
#include <Connectors/ConnectorFactory.h>
#include <Connectors/Block3DConnectorFactory.h>
#include <Data/D3Q27EsoTwist3DSplittedVector.h>
#include <Data/D3Q27EsoTwist3DSplittedVectorEx.h>
#include <Data/DataSet3D.h>
#include <Data/DistributionArray3D.h>
#include <Data/EsoTwist3D.h>
#include <Data/EsoTwistD3Q27System.h>
#include <Data/EsoTwistD3Q27SparseData.h>
#include <Data/VoidData3D.h>
#include <Grid/Block3D.h>
//#include <Grid/BoostSerializationClassExportHelper.h>
#include <Grid/CalculationManager.h>
#include <Grid/Calculator.h>
//#include <Grid/Calculator2.h>
#include <Grid/Grid3D.h>
#include <Grid/Grid3DSystem.h>
#include <Interactors/D3Q27Interactor.h>
#include <Interactors/D3Q27TriFaceMeshInteractor.h>
#include <Interactors/Interactor3D.h>
#include <Interactors/InteractorsHelper.h>
#include <CoProcessors/WriteBlocksCoProcessor.h>
#include <CoProcessors/AdjustForcingCoProcessor.h>
#include <CoProcessors/CalculateForcesCoProcessor.h>
#include <CoProcessors/WriteMacroscopicQuantitiesCoProcessor.h>
#include <CoProcessors/WriteBoundaryConditionsCoProcessor.h>
//#include <CoProcessors/PathLineCoProcessor.h>
//#include <CoProcessors/PathLineCoProcessorMcpart.h>
#include <CoProcessors/PressureDifferenceCoProcessor.h>
#include <CoProcessors/EmergencyExitCoProcessor.h>
#include <CoProcessors/NUPSCounterCoProcessor.h>
//#include <CoProcessors/Particles.h>
#include <CoProcessors/CoProcessor.h>
#include <CoProcessors/RestartCoProcessor.h>
#include <CoProcessors/TurbulenceIntensityCoProcessor.h>
#include <CoProcessors/AverageValuesCoProcessor.h>
#include <CoProcessors/DecreaseViscosityCoProcessor.h>
#include <CoProcessors/TimeseriesCoProcessor.h>
#include <CoProcessors/ShearStressCoProcessor.h>
#include <CoProcessors/QCriterionCoProcessor.h>
#include <CoProcessors/InSituVTKCoProcessor.h>
//#include <CoProcessors/MeanValuesCoProcessor.h>
#include <CoProcessors/TimeAveragedValuesCoProcessor.h>
#include <CoProcessors/InSituCatalystCoProcessor.h>
#include <LineTimeSeriesCoProcessor.h>
#include <IntegrateValuesHelper.h>
//#include <LBM/D3Q27CompactInterpolationProcessor.h>
#include <LBM/IncompressibleOffsetInterpolationProcessor.h>
#include <LBM/CompressibleOffsetInterpolationProcessor.h>
#include <LBM/InterpolationHelper.h>
#include <LBM/InterpolationProcessor.h>
//#include <LBM/D3Q27OffsetInterpolationProcessor.h>
#include <LBM/D3Q27System.h>
#include <LBM/ICell.h>
#include <LBM/InterpolationProcessor.h>
#include <LBM/LBMKernel.h>
#include <IncompressibleCumulantWithSpongeLayerLBMKernel.h>
#include <LBM/LBMKernel.h>
#include <LBM/IncompressibleCumulantLBMKernel.h>
#include <LBM/CompressibleCumulantLBMKernel.h>
#include <LBM/InitDensityLBMKernel.h>
#include <LBM/VoidLBMKernel.h>
#include <LBM/LBMSystem.h>
#include <LBM/LBMSystems.h>
#include <LBM/LBMUnitConverter.h>
#include <numerics/geometry3d/CoordinateTransformation3D.h>
#include <numerics/geometry3d/GbCuboid3D.h>
#include <numerics/geometry3d/GbCylinder3D.h>
#include <numerics/geometry3d/GbHalfSpace3D.h>
#include <numerics/geometry3d/GbHalfSpaceKrischan3D.h>
#include <numerics/geometry3d/GbLine3D.h>
#include <numerics/geometry3d/GbMeshTools3D.h>
#include <numerics/geometry3d/GbObject3D.h>
#include <numerics/geometry3d/GbObject3DManager.h>
#include <numerics/geometry3d/GbObjectGroup3D.h>
#include <numerics/geometry3d/GbPoint3D.h>
#include <numerics/geometry3d/GbPolygon3D.h>
#include <numerics/geometry3d/GbQuadFaceMesh3D.h>
#include <numerics/geometry3d/GbSphere3D.h>
#include <numerics/geometry3d/GbSystem3D.h>
#include <numerics/geometry3d/GbTriangle3D.h>
#include <numerics/geometry3d/GbTriangularMesh3D.h>
#include <numerics/geometry3d/GbTriFaceMesh3D.h>
#include <numerics/geometry3d/GbVector3D.h>
#include <numerics/geometry3d/GbVoxelMatrix3D.h>
#include <numerics/geometry3d/creator/GbCuboid3DCreator.h>
#include <numerics/geometry3d/creator/GbCylinder3DCreator.h>
#include <numerics/geometry3d/creator/GbLine3DCreator.h>
#include <numerics/geometry3d/creator/GbObject3DCreator.h>
#include <numerics/geometry3d/creator/GbObject3DFactory.h>
#include <numerics/geometry3d/creator/GbPoint3DCreator.h>
#include <numerics/geometry3d/creator/GbPolygon3DCreator.h>
#include <numerics/geometry3d/creator/GbQuadFaceMesh3DCreator.h>
#include <numerics/geometry3d/creator/GbSphere3DCreator.h>
#include <numerics/geometry3d/creator/GbTriangle3DCreator.h>
#include <numerics/geometry3d/creator/GbTriangularMesh3DCreator.h>
#include <numerics/geometry3d/creator/GbTriFaceMesh3DCreator.h>
#include <numerics/geometry3d/creator/GbVoxelMatrix3DCreator.h>
// #include <numerics/geometry3d/examples/stl2inp/QDefineUniformMesh.h>
// #include <numerics/geometry3d/examples/stl2inp/stl2inp.h>
// #include <numerics/geometry3d/fem/FeAdhocTriFaceMesh3D.h>
// #include <numerics/geometry3d/fem/FeHalfDisc3D.h>
// #include <numerics/geometry3d/fem/FePlateTriangularMesh3D.h>
// #include <numerics/geometry3d/fem/FePoint3D.h>
// #include <numerics/geometry3d/fem/FeRing3D.h>
// #include <numerics/geometry3d/fem/FeTriFaceMesh3D.h>
// #include <numerics/geometry3d/fem/creator/FeTriFaceMesh3DCreator.h>
#include <numerics/geometry3d/KdTree/KdNode.h>
#include <numerics/geometry3d/KdTree/KdRay.h>
#include <numerics/geometry3d/KdTree/KdSplitCandidate.h>
#include <numerics/geometry3d/KdTree/KdSplitCandidateManager.h>
#include <numerics/geometry3d/KdTree/KdTree.h>
#include <numerics/geometry3d/KdTree/KdUtilities.h>
#include <numerics/geometry3d/KdTree/intersectionhandler/KdCountLineIntersectionHandler.h>
#include <numerics/geometry3d/KdTree/intersectionhandler/KdCountRayIntersectionHandler.h>
#include <numerics/geometry3d/KdTree/intersectionhandler/KdLineIntersectionHandler.h>
#include <numerics/geometry3d/KdTree/intersectionhandler/KdRayIntersectionHandler.h>
#include <numerics/geometry3d/KdTree/splitalgorithms/KdSAHSplit.h>
#include <numerics/geometry3d/KdTree/splitalgorithms/KdSpatiallMedianSplit.h>
#include <numerics/geometry3d/KdTree/splitalgorithms/KdSplitAlgorithm.h>
// #include <numerics/geometry3d/presentation/QGbCuboid3DInstrument.h>
// #include <numerics/geometry3d/presentation/QGbCylinder3DInstrument.h>
// #include <numerics/geometry3d/presentation/QGbObject3DInstrument.h>
// #include <numerics/geometry3d/presentation/QGbSphere3DInstrument.h>
// #include <numerics/geometry3d/presentation/QVTKGbObject3DViewer.h>
// #include <numerics/geometry3d/presentation/vtkGbCuboid3D.h>
// #include <numerics/geometry3d/presentation/vtkGbCylinder3D.h>
// #include <numerics/geometry3d/presentation/vtkGbSphere3D.h>
// #include <numerics/geometry3d/presentation/vtkGbTriangularMesh3D.h>

#include <Parallel/Communicator.h>
#include <Parallel/LoadBalancer.h>
#include <Parallel/MetisPartitioner.h>
#include <Parallel/MPICommunicator.h>
#include <Parallel/NullCommunicator.h>
#include <Parallel/PriorityQueueDecompositor.h>
#include <Parallel/SimpleGeometricPartitioner.h>
#include <Parallel/Synchronizer.h>
#include <Parallel/ZoltanPartitioner.h>
#include <Parallel/BlocksDistributor.h>
#include <ZoltanPartitioningGridVisitor.h>
#include <Utilities/MathUtil.hpp>
#include <Utilities/MemoryUtil.h>
#include <Utilities/StringUtil.hpp>
#include <Utilities/ConfigurationFile.hpp>
#include <Utilities/VoxelMatrixUtil.hpp>
#include <Utilities/ChangeRandomQs.hpp>
#include <Utilities/ConfigFileReader.h>
#include <Visitors/Block3DVisitor.h>
#include <Visitors/CreateTransmittersHelper.h>
#include <Visitors/InitDistributionsBlockVisitor.h>
#include <Visitors/SetConnectorsBlockVisitor.h>
#include <Visitors/SetUndefinedNodesBlockVisitor.h>
#include <Visitors/GenBlocksGridVisitor.h>
#include <Visitors/Grid3DVisitor.h>
#include <Visitors/MetisPartitioningGridVisitor.h>
#include <Visitors/OverlapBlockVisitor.h>
#include <Visitors/PQueuePartitioningGridVisitor.h>
#include <Visitors/RatioBlockVisitor.h>
#include <Visitors/RatioSmoothBlockVisitor.h>
#include <Visitors/RefineCrossAndInsideGbObjectBlockVisitor.h>
#include <Visitors/RefineInterGbObjectsVisitor.h>
#include <Visitors/CoarsenCrossAndInsideGbObjectBlockVisitor.h>
#include <Visitors/SetInterpolationDirsBlockVisitor.h>
#include <Visitors/SetKernelBlockVisitor.h>
#include <Visitors/SetForcingBlockVisitor.h>
#include <Visitors/SetSpongeLayerBlockVisitor.h>
#include <Visitors/SetSolidOrTransBlockVisitor.h>
#include <Visitors/RenumberBlockVisitor.h>
#include <Visitors/ConnectorBlockVisitor.h>
#include <Visitors/ViscosityBlockVisitor.h>
#include <Visitors/BoundaryConditionsBlockVisitor.h>
#include <Visitors/BoundaryConditionsBlockVisitor.h>
#include <Visitors/ChangeBoundaryDensityBlockVisitor.h>
#include <InitDistributionsFromFileBlockVisitor.h>
#include <CheckRatioBlockVisitor.h>
#include <SpongeLayerBlockVisitor.h>

#include <Visitors/RefineCrossAndInsideGbObjectHelper.h>
#include <RefineAroundGbObjectHelper.h>

#if defined VF_FETOL
   #include <FETOL/FETOLCalculator.h>
   #include <FETOL/FETOLCommunicator.h>
   #include <FETOL/FETOLSetConnectorsBlockVisitor.h>
   #include <FETOL/FETOLTransmitterBondPool.h>   
#endif

