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
    std::shared_ptr<vf::mpi::Communicator> comm = vf::mpi::MPICommunicator::getInstance();
    int myid                                        = comm->getProcessID();


    // bounding box
    double g_minX1 = 0;
    double g_minX2 = 0;
    double g_minX3 = 0;

    double g_maxX1 = 1;
    double g_maxX2 = 1;
    double g_maxX3 = 10;

    int blockNX[3] = { 16, 16, 16 };
    double dx = 1./32.;

    double d_part = 0.25;
    double r_p    = d_part / 2.0;

    //int blockNX[3] = { 10, 10, 10 };
    //double dx      = 0.05;


    double nuLB = 1e-2;

    SPtr<LBMKernel> kernel   = make_shared<IBcumulantK17LBMKernel>();
    SPtr<BCProcessor> bcProc = make_shared<BCProcessor>();
    kernel->setBCProcessor(bcProc);

    SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
    noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));
    //////////////////////////////////////////////////////////////////////////////////
    // BC visitor
    BoundaryConditionsBlockVisitor bcVisitor;
    bcVisitor.addBC(noSlipBCAdapter);

    SPtr<Grid3D> grid = make_shared<Grid3D>(comm);
    grid->setPeriodicX1(true);
    grid->setPeriodicX2(true);
    grid->setPeriodicX3(false);
    grid->setDeltaX(dx);
    grid->setBlockNX(blockNX[0], blockNX[1], blockNX[2]);

    string outputPath = "f:/temp/FallingSpheresTest";

    UbSystem::makeDirectory(outputPath);
    UbSystem::makeDirectory(outputPath + "/liggghts");

    SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::DIR_MMM, MetisPartitioner::RECURSIVE));
    
    SPtr<GbObject3D> gridCube = make_shared <GbCuboid3D>(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3);
    if (myid == 0)
        GbSystem3D::writeGeoObject(gridCube.get(), outputPath + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

    GenBlocksGridVisitor genBlocks(gridCube);
    grid->accept(genBlocks);

    SPtr<CoProcessor> ppblocks =
        make_shared <WriteBlocksCoProcessor>(grid, SPtr<UbScheduler>(new UbScheduler(1)), outputPath,
                                                          WbWriterVtkXmlBinary::getInstance(), comm);
    ppblocks->process(0);
    ppblocks.reset();

    double dx2 = 2.0 * dx;
    GbCuboid3DPtr wallZmin(
        new GbCuboid3D(g_minX1 - dx2, g_minX2 - dx2, g_minX3 - dx2, g_maxX1 + dx2, g_maxX2 + dx2, g_minX3));
    GbSystem3D::writeGeoObject(wallZmin.get(), outputPath + "/geo/wallZmin", WbWriterVtkXmlASCII::getInstance());
    GbCuboid3DPtr wallZmax(
        new GbCuboid3D(g_minX1 - dx2, g_minX2 - dx2, g_maxX3, g_maxX1 + dx2, g_maxX2 + dx2, g_maxX3 + dx2));
    GbSystem3D::writeGeoObject(wallZmax.get(), outputPath + "/geo/wallZmax", WbWriterVtkXmlASCII::getInstance());

    SPtr<D3Q27Interactor> wallZminInt(new D3Q27Interactor(wallZmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
    SPtr<D3Q27Interactor> wallZmaxInt(new D3Q27Interactor(wallZmax, grid, noSlipBCAdapter, Interactor3D::SOLID));

    InteractorsHelper intHelper(grid, metisVisitor, true);
    intHelper.addInteractor(wallZminInt);
    intHelper.addInteractor(wallZmaxInt);
    intHelper.selectBlocks();

    SetKernelBlockVisitor kernelVisitor(kernel, nuLB, 1e9, 1e9);
    grid->accept(kernelVisitor);

    intHelper.setBC();

    InitDistributionsBlockVisitor initVisitor;
    grid->accept(initVisitor);

    SPtr<UbScheduler> lScheduler = make_shared<UbScheduler>(1);
    string inFile1 = "d:/Projects/VirtualFluids_Develop/apps/cpu/FallingSphere/in.lbdem";
    string inFile2 = "d:/Projects/VirtualFluids_Develop/apps/cpu/FallingSphere/in2.lbdem";
    MPI_Comm mpi_comm = *(MPI_Comm*)(comm->getNativeCommunicator());
    LiggghtsCouplingWrapper wrapper(argv, mpi_comm);


 
    // SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(r_p, 1.480, 2060, r_p/dx);
    //SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(r_p, LBMUnitConverter::AIR_20C, r_p / dx);
    SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(r_p, 0.1, 1000, r_p / dx, 0.01);
    std::cout << units->toString() << std::endl;

    double v_frac = 0.1;
    double dt_phys   = units->getFactorTimeLbToW();
    int demSubsteps = 10;
    double dt_dem   = dt_phys / (double)demSubsteps;
    int vtkSteps    = 100;
    string demOutDir = outputPath; 

    wrapper.execCommand("echo none");

    wrapper.setVariable("d_part", d_part);
    //wrapper.setVariable("r_part", d_part/2.);
    //wrapper.setVariable("v_frac", v_frac);

    wrapper.execFile((char*)inFile1.c_str());
 
    //// set timestep and output directory
    wrapper.setVariable("t_step", dt_dem);
    wrapper.setVariable("dmp_stp", vtkSteps * demSubsteps);
    wrapper.setVariable("dmp_dir", demOutDir);

    wrapper.execFile((char *)inFile2.c_str());
    wrapper.runUpto(demSubsteps - 1);
  
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


    // write data for visualization of macroscopic quantities
    SPtr<UbScheduler> visSch(new UbScheduler(vtkSteps));
    SPtr<WriteMacroscopicQuantitiesCoProcessor> writeMQCoProcessor(
        new WriteMacroscopicQuantitiesCoProcessor(grid, visSch, outputPath, WbWriterVtkXmlBinary::getInstance(),
                                                  SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));

    int endTime = 3000; //20;
    SPtr<Calculator> calculator(new BasicCalculator(grid, lScheduler, endTime));
    calculator->addCoProcessor(lcCoProcessor);
    calculator->addCoProcessor(writeMQCoProcessor);

    if (myid == 0) UBLOG(logINFO, "Simulation-start");
    calculator->calculate();
    if (myid == 0) UBLOG(logINFO, "Simulation-end");


    return 0;
}
