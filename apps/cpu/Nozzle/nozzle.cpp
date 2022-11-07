#include <iostream>
#include <string>
#include <memory>

#include "VirtualFluids.h"

#include "LiggghtsCouplingCoProcessor.h"
#include "LiggghtsCouplingWrapper.h"
#include "IBcumulantK17LBMKernel.h"

using namespace std;


int main(int argc, char *argv[])
{
    //Sleep(30000);

    std::shared_ptr<vf::mpi::Communicator> comm = vf::mpi::MPICommunicator::getInstance();
    int myid = comm->getProcessID();


    // bounding box
    //double g_minX1 = -1341.81e-3;
    //double g_minX2 =  348.087e-3;
    //double g_minX3 = -210e-3;

    //double g_maxX1 = -1260.81e-3;
    //double g_maxX2 =  429.087e-3;
    //double g_maxX3 =  214.5e-3;

    double g_minX1 = -1341.81e-3 + 10e-3;
    double g_minX2 =  0.360872;
    double g_minX3 = -210e-3;

    double g_maxX1 = -1260.81e-3 - 10e-3;
    double g_maxX2 =  0.416302;
    double g_maxX3 =  210e-3;

    int blockNX[3] = { 10, 10, 10 };

    double dx = 1e-3;

    double nuLB = 1e-3;
    double uLB  = -0.01;
    double rhoLB = 0.0;

    SPtr<LBMKernel> kernel   = make_shared<IBcumulantK17LBMKernel>();
    //SPtr<LBMKernel> kernel   = make_shared<CumulantK17LBMKernel>();
    SPtr<BCProcessor> bcProc = make_shared<BCProcessor>();
    kernel->setBCProcessor(bcProc);

    SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
    noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

    mu::Parser fct;
    fct.SetExpr("U");
    fct.DefineConst("U", uLB);
    SPtr<BCAdapter> inflowBCAdapter(new VelocityBCAdapter(false, false, true, fct, 0, BCFunction::INFCONST));
    inflowBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityBCAlgorithm()));

    SPtr<BCAdapter> outflowBCAdapter(new DensityBCAdapter(rhoLB));
    outflowBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));
    //////////////////////////////////////////////////////////////////////////////////
    // BC visitor
    BoundaryConditionsBlockVisitor bcVisitor;
    bcVisitor.addBC(noSlipBCAdapter);
    bcVisitor.addBC(inflowBCAdapter);
    bcVisitor.addBC(outflowBCAdapter);

    SPtr<Grid3D> grid = make_shared<Grid3D>(comm);
    grid->setPeriodicX1(true);
    grid->setPeriodicX2(true);
    grid->setPeriodicX3(false);
    grid->setDeltaX(dx);
    grid->setBlockNX(blockNX[0], blockNX[1], blockNX[2]);

    string geoPath = "d:/Projects/TRR277/Project/WP4/NozzleGeo";

    string outputPath = "f:/temp/NozzleFlowTestSerial";
    UbSystem::makeDirectory(outputPath);
    UbSystem::makeDirectory(outputPath + "/liggghts");

    if (myid == 0) {
        stringstream logFilename;
        logFilename << outputPath + "/logfile" + UbSystem::toString(UbSystem::getTimeStamp()) + ".txt";
        UbLog::output_policy::setStream(logFilename.str());
    }

    SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::DIR_MMM, MetisPartitioner::RECURSIVE));
    
    SPtr<GbObject3D> gridCube = make_shared <GbCuboid3D>(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3);
    if (myid == 0)
        GbSystem3D::writeGeoObject(gridCube.get(), outputPath + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

    GenBlocksGridVisitor genBlocks(gridCube);
    grid->accept(genBlocks);

    //geo
    //////////////////////////////////////////////////////////
    int accuracy = Interactor3D::EDGES;
    ///////////////////////////////////
    SPtr<GbTriFaceMesh3D> meshNozzleAirDistributor = std::make_shared<GbTriFaceMesh3D>();
    if (myid == 0) UBLOG(logINFO, "Read meshNozzleAirDistributor:start");
    meshNozzleAirDistributor->readMeshFromSTLFileASCII(geoPath + "/01_Nozzle_Air_Distributor.stl", false);
    if (myid == 0) UBLOG(logINFO, "Read meshNozzleAirDistributor:end");
    if (myid == 0) GbSystem3D::writeGeoObject(meshNozzleAirDistributor.get(), outputPath + "/geo/meshNozzleAirDistributor", WbWriterVtkXmlBinary::getInstance());
    SPtr<Interactor3D> intrNozzleAirDistributor = std::make_shared<D3Q27TriFaceMeshInteractor>(meshNozzleAirDistributor, grid, noSlipBCAdapter, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy);
    ///////////////////////////////////////////////////////////
    SPtr<GbTriFaceMesh3D> meshNozzleAirInlet = std::make_shared<GbTriFaceMesh3D>();
    if (myid == 0) UBLOG(logINFO, "Read meshNozzleAirInlet:start");
    meshNozzleAirInlet->readMeshFromSTLFileASCII(geoPath + "/02_Nozzle_Air_Inlet.stl", false);
    if (myid == 0) UBLOG(logINFO, "Read meshNozzleAirInlet:end");
    if (myid == 0) GbSystem3D::writeGeoObject(meshNozzleAirInlet.get(), outputPath + "/geo/meshNozzleAirInlet", WbWriterVtkXmlBinary::getInstance());
    SPtr<Interactor3D> intrNozzleAirInlet = std::make_shared<D3Q27TriFaceMeshInteractor>(meshNozzleAirInlet, grid, noSlipBCAdapter, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy);
    ///////////////////////////////////////////////////////////
    SPtr<GbTriFaceMesh3D> meshNozzleSpacer = std::make_shared<GbTriFaceMesh3D>();
    if (myid == 0) UBLOG(logINFO, "Read meshNozzleSpacer:start");
    meshNozzleSpacer->readMeshFromSTLFileASCII(geoPath + "/03_Nozzle_Spacer.stl", true);
    if (myid == 0) UBLOG(logINFO, "Read meshNozzleSpacer:end");
    if (myid == 0) GbSystem3D::writeGeoObject(meshNozzleSpacer.get(), outputPath + "/geo/meshNozzleSpacer", WbWriterVtkXmlBinary::getInstance());
    SPtr<Interactor3D> intrNozzleSpacer = std::make_shared<D3Q27TriFaceMeshInteractor>(meshNozzleSpacer, grid, noSlipBCAdapter, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy);
    ///////////////////////////////////////////////////////////
    SPtr<GbTriFaceMesh3D> meshNozzleAccDistributor = std::make_shared<GbTriFaceMesh3D>();
    if (myid == 0) UBLOG(logINFO, "Read meshNozzleAccDistributor:start");
    meshNozzleAccDistributor->readMeshFromSTLFileASCII(geoPath + "/04_Nozzle_Acc_Distributor.stl", false);
    if (myid == 0) UBLOG(logINFO, "Read meshNozzleAccDistributor:end");
    if (myid == 0) GbSystem3D::writeGeoObject(meshNozzleAccDistributor.get(), outputPath + "/geo/meshNozzleAccDistributor", WbWriterVtkXmlBinary::getInstance());
    SPtr<Interactor3D> intrNozzleAccDistributor = std::make_shared<D3Q27TriFaceMeshInteractor>(meshNozzleAccDistributor, grid, noSlipBCAdapter, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy);
    ///////////////////////////////////////////////////////////
    SPtr<GbTriFaceMesh3D> meshNozzleAccInlet = std::make_shared<GbTriFaceMesh3D>();
    if (myid == 0) UBLOG(logINFO, "Read meshNozzleAccInlet:start");
    meshNozzleAccInlet->readMeshFromSTLFileASCII(geoPath + "/05_Nozzle_Acc_Inlet.stl", false);
    if (myid == 0) UBLOG(logINFO, "Read meshNozzleAccInlet:end");
    if (myid == 0) GbSystem3D::writeGeoObject(meshNozzleAccInlet.get(), outputPath + "/geo/meshNozzleAccInlet", WbWriterVtkXmlBinary::getInstance());
    SPtr<Interactor3D> intrNozzleAccInlet = std::make_shared<D3Q27TriFaceMeshInteractor>(meshNozzleAccInlet, grid, noSlipBCAdapter, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy);
    ///////////////////////////////////////////////////////////
    SPtr<GbTriFaceMesh3D> meshNozzleVolcanNozzle1 = std::make_shared<GbTriFaceMesh3D>();
    if (myid == 0) UBLOG(logINFO, "Read meshNozzleVolcanNozzle1:start");
    meshNozzleVolcanNozzle1->readMeshFromSTLFileBinary(geoPath + "/06_1_Nozzle_Volcan_Nozzle.stl", true);
    if (myid == 0) UBLOG(logINFO, "Read meshNozzleVolcanNozzle1:end");
    if (myid == 0) GbSystem3D::writeGeoObject(meshNozzleVolcanNozzle1.get(), outputPath + "/geo/meshNozzleVolcanNozzle1", WbWriterVtkXmlBinary::getInstance());
    SPtr<Interactor3D> intrNozzleVolcanNozzle1 = std::make_shared<D3Q27TriFaceMeshInteractor>(meshNozzleVolcanNozzle1, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES);
    ///////////////////////////////////////////////////////////
    SPtr<GbTriFaceMesh3D> meshNozzleVolcanNozzle2 = std::make_shared<GbTriFaceMesh3D>();
    if (myid == 0) UBLOG(logINFO, "Read meshNozzleVolcanNozzle2:start");
    meshNozzleVolcanNozzle2->readMeshFromSTLFileBinary(geoPath + "/06_2_Nozzle_Volcan_Nozzle.stl", true);
    if (myid == 0) UBLOG(logINFO, "Read meshNozzleVolcanNozzle2:end");
    if (myid == 0) GbSystem3D::writeGeoObject(meshNozzleVolcanNozzle2.get(), outputPath + "/geo/meshNozzleVolcanNozzle2", WbWriterVtkXmlBinary::getInstance());
    SPtr<Interactor3D> intrNozzleVolcanNozzle2 = std::make_shared<D3Q27TriFaceMeshInteractor>(meshNozzleVolcanNozzle2, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES);
    ///////////////////////////////////////////////////////////
    //box
    SPtr<D3Q27Interactor> intrBox = SPtr<D3Q27Interactor>(new D3Q27Interactor(gridCube, grid, noSlipBCAdapter, Interactor3D::INVERSESOLID));
    ///////////////////////////////////////////////////////////
    //inflow
    GbCylinder3DPtr geoInflow(new GbCylinder3D(-1.30181+0.0005, 0.390872-0.00229, 0.20105, -1.30181+0.0005, 0.390872-0.00229, 0.23, 0.013));
    if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), outputPath + "/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());
    SPtr<D3Q27Interactor> intrInflow = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, inflowBCAdapter, Interactor3D::SOLID));
    ///////////////////////////////////////////////////////////
    //outflow
    GbCylinder3DPtr geoOutflow(new GbCylinder3D(-1.30181+0.0005, 0.390872-0.00229, -0.22, -1.30181+0.0005, 0.390872-0.00229, -0.21, 0.013));
    if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), outputPath + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());
    SPtr<D3Q27Interactor> intrOutflow = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, inflowBCAdapter, Interactor3D::SOLID));
    ///////////////////////////////////////////////////////////

    InteractorsHelper intHelper(grid, metisVisitor, true);
    intHelper.addInteractor(intrBox);
    intHelper.addInteractor(intrInflow);
    intHelper.addInteractor(intrOutflow);
    intHelper.addInteractor(intrNozzleAirDistributor);
    intHelper.addInteractor(intrNozzleAirInlet);
    intHelper.addInteractor(intrNozzleSpacer);
    intHelper.addInteractor(intrNozzleAccDistributor);
    intHelper.addInteractor(intrNozzleAccInlet);
    intHelper.addInteractor(intrNozzleVolcanNozzle1);
    intHelper.addInteractor(intrNozzleVolcanNozzle2);


    intHelper.selectBlocks();

    SPtr<CoProcessor> ppblocks = make_shared<WriteBlocksCoProcessor>(
         grid, SPtr<UbScheduler>(new UbScheduler(1)), outputPath, WbWriterVtkXmlBinary::getInstance(), comm);
     ppblocks->process(0);
     ppblocks.reset();

     if (myid == 0) UBLOG(logINFO, Utilities::toString(grid, comm->getNumberOfProcesses()));


    SetKernelBlockVisitor kernelVisitor(kernel, nuLB, comm->getNumberOfProcesses());
    grid->accept(kernelVisitor);

    intHelper.setBC();

    InitDistributionsBlockVisitor initVisitor;
    grid->accept(initVisitor);

  
    string inFile1 = "d:/Projects/VirtualFluids_Develop/apps/cpu/Nozzle/in.nozzle";
    //string inFile2 = "d:/Projects/VirtualFluids_LIGGGHTS_coupling/apps/cpu/LiggghtsApp/in2.lbdem";
    MPI_Comm mpi_comm = *(MPI_Comm*)(comm->getNativeCommunicator());
    LiggghtsCouplingWrapper wrapper(argv, mpi_comm);


    double d_part = 1e-3;
 
    // SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(r_p, 1.480, 2060, r_p/dx);
    //SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(r_p, LBMUnitConverter::AIR_20C, r_p / dx);
    SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(d_part, 1., 1000, d_part / dx, 0.01);
    if (myid == 0) std::cout << units->toString() << std::endl;

    //return 0;

    double v_frac = 0.1;
    double dt_phys   = units->getFactorTimeLbToW();
    int demSubsteps = 10;
    double dt_dem   = dt_phys / (double)demSubsteps;
    int vtkSteps    = 100;
    string demOutDir = outputPath + "/liggghts";

    //wrapper.execCommand("echo none");

    //wrapper.execFile((char*)inFile1.c_str());

    //// set timestep and output directory
    wrapper.setVariable("t_step", dt_dem);
    wrapper.setVariable("dmp_stp", vtkSteps * demSubsteps);
    wrapper.setVariable("dmp_dir", demOutDir);

    wrapper.execFile((char *)inFile1.c_str());
    wrapper.runUpto(demSubsteps - 1);
    //wrapper.runUpto(1000);

    SPtr<UbScheduler> lScheduler = make_shared<UbScheduler>(1); 
    SPtr<LiggghtsCouplingCoProcessor> lcCoProcessor =
        make_shared<LiggghtsCouplingCoProcessor>(grid, lScheduler, comm, wrapper, demSubsteps, units);

    // boundary conditions grid
    {
        SPtr<UbScheduler> geoSch(new UbScheduler(1));
        SPtr<WriteBoundaryConditionsCoProcessor> ppgeo(new WriteBoundaryConditionsCoProcessor(
            grid, geoSch, outputPath, WbWriterVtkXmlBinary::getInstance(), comm));
        ppgeo->process(0);
        ppgeo.reset();
    }

    grid->accept(bcVisitor);

    OneDistributionSetConnectorsBlockVisitor setConnsVisitor(comm);
    grid->accept(setConnsVisitor);

    int numOfThreads          = 18;
    omp_set_num_threads(numOfThreads);

    SPtr<UbScheduler> nupsSch = std::make_shared<UbScheduler>(10, 10, 100);
    SPtr<NUPSCounterCoProcessor> nupsCoProcessor = make_shared<NUPSCounterCoProcessor>(grid, nupsSch, numOfThreads, comm);

    //// write data for visualization of macroscopic quantities
    SPtr < UbScheduler> visSch(new UbScheduler(vtkSteps));
    //SPtr<UbScheduler> visSch(new UbScheduler(1, 8700, 8800));
    visSch->addSchedule(1, 8700, 8800);
    SPtr<WriteMacroscopicQuantitiesCoProcessor> writeMQCoProcessor(
        new WriteMacroscopicQuantitiesCoProcessor(grid, visSch, outputPath, WbWriterVtkXmlBinary::getInstance(),
                                                  SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));
    writeMQCoProcessor->process(0);

    int endTime = 10000000; //20;
    SPtr<Calculator> calculator(new BasicCalculator(grid, lScheduler, endTime));
    calculator->addCoProcessor(nupsCoProcessor);
    calculator->addCoProcessor(lcCoProcessor);
    calculator->addCoProcessor(writeMQCoProcessor);

    if (myid == 0) UBLOG(logINFO, "Simulation-start");
    calculator->calculate();
    if (myid == 0) UBLOG(logINFO, "Simulation-end");


    return 0;
}
