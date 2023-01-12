#include <iostream>
#include <string>
#include <memory>

#include "VirtualFluids.h"

using namespace std;

void run(string configname)
{
    using namespace vf::lbm::dir;

    try {

        //Sleep(30000);

        vf::basics::ConfigurationFile config;
        config.load(configname);

        string pathname            = config.getValue<string>("pathname");
        string pathGeo             = config.getValue<string>("pathGeo");
        string geoFile             = config.getValue<string>("geoFile");
        int numOfThreads           = config.getValue<int>("numOfThreads");
        vector<int> blocknx        = config.getVector<int>("blocknx");
        vector<double> boundingBox = config.getVector<double>("boundingBox");
        // vector<double>  length = config.getVector<double>("length");
        double uLB = config.getValue<double>("uLB");
        // double uF2                         = config.getValue<double>("uF2");
        double nuL             = config.getValue<double>("nuL");
        double nuG             = config.getValue<double>("nuG");
        double densityRatio    = config.getValue<double>("densityRatio");
        double sigma           = config.getValue<double>("sigma");
        int interfaceWidth = config.getValue<int>("interfaceWidth");
        //double radius          = config.getValue<double>("radius");
        double theta           = config.getValue<double>("contactAngle");
        double gr              = config.getValue<double>("gravity");
        double phiL            = config.getValue<double>("phi_L");
        double phiH            = config.getValue<double>("phi_H");
        double tauH            = config.getValue<double>("Phase-field Relaxation");
        double mob             = config.getValue<double>("Mobility");

        double endTime     = config.getValue<double>("endTime");
        double outTime     = config.getValue<double>("outTime");
        double availMem    = config.getValue<double>("availMem");
        int refineLevel    = config.getValue<int>("refineLevel");
        double Re          = config.getValue<double>("Re");
        double dx          = config.getValue<double>("dx");
        bool logToFile     = config.getValue<bool>("logToFile");
        double restartStep = config.getValue<double>("restartStep");
        double cpStart     = config.getValue<double>("cpStart");
        double cpStep      = config.getValue<double>("cpStep");
        bool newStart      = config.getValue<bool>("newStart");

        double beta = 12 * sigma / interfaceWidth;
        double kappa = 1.5 * interfaceWidth * sigma;

        SPtr<vf::mpi::Communicator> comm = vf::mpi::MPICommunicator::getInstance();
        int myid                = comm->getProcessID();

        if (myid == 0)
            UBLOG(logINFO, "Jet Breakup: Start!");

        if (logToFile) {
#if defined(__unix__)
            if (myid == 0) {
                const char *str = pathname.c_str();
                mkdir(str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            }
#endif

            if (myid == 0) {
                stringstream logFilename;
                logFilename << pathname + "/logfile" + UbSystem::toString(UbSystem::getTimeStamp()) + ".txt";
                UbLog::output_policy::setStream(logFilename.str());
            }
        }

        // Sleep(30000);

        // LBMReal dLB = 0; // = length[1] / dx;
        LBMReal rhoLB = 0.0;
        LBMReal nuLB  = nuL; //(uLB*dLB) / Re;

        SPtr<LBMUnitConverter> conv(new LBMUnitConverter());

        //const int baseLevel = 0;

        SPtr<LBMKernel> kernel;

        //kernel = SPtr<LBMKernel>(new MultiphaseScratchCumulantLBMKernel());
        //kernel = SPtr<LBMKernel>(new MultiphaseCumulantLBMKernel());
        //kernel = SPtr<LBMKernel>(new MultiphaseTwoPhaseFieldsCumulantLBMKernel());
        //kernel = SPtr<LBMKernel>(new MultiphaseTwoPhaseFieldsVelocityCumulantLBMKernel());
       // kernel = SPtr<LBMKernel>(new MultiphaseTwoPhaseFieldsPressureFilterLBMKernel());
        kernel = SPtr<LBMKernel>(new MultiphasePressureFilterLBMKernel());

        kernel->setWithForcing(true);
        kernel->setForcingX1(0.0);
        kernel->setForcingX2(gr);
        kernel->setForcingX3(0.0);

        kernel->setPhiL(phiL);
        kernel->setPhiH(phiH);
        kernel->setPhaseFieldRelaxation(tauH);
        kernel->setMobility(mob);

        //nuL, nuG, densityRatio, beta, kappa, theta,

        kernel->setCollisionFactorMultiphase(nuL, nuG);
        kernel->setDensityRatio(densityRatio);
        kernel->setMultiphaseModelParameters(beta, kappa);
        kernel->setContactAngle(theta);
        kernel->setInterfaceWidth(interfaceWidth);

        SPtr<BCProcessor> bcProc(new BCProcessor());
        // BCProcessorPtr bcProc(new ThinWallBCProcessor());

        kernel->setBCProcessor(bcProc);

        SPtr<Grid3D> grid(new Grid3D(comm));
         //grid->setPeriodicX1(true);
         //grid->setPeriodicX2(true);
         //grid->setPeriodicX3(true);
        grid->setGhostLayerWidth(2);

       
        SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, DIR_MMM, MetisPartitioner::RECURSIVE));

        //////////////////////////////////////////////////////////////////////////
        // restart
        SPtr<UbScheduler> rSch(new UbScheduler(cpStep, cpStart));
        //SPtr<MPIIORestartCoProcessor> rcp(new MPIIORestartCoProcessor(grid, rSch, pathname, comm));
        SPtr<MPIIOMigrationCoProcessor> rcp(new MPIIOMigrationCoProcessor(grid, rSch, metisVisitor, pathname, comm));
        //SPtr<MPIIOMigrationBECoProcessor> rcp(new MPIIOMigrationBECoProcessor(grid, rSch, pathname, comm));
        //rcp->setNu(nuLB);
        //rcp->setNuLG(nuL, nuG);
        //rcp->setDensityRatio(densityRatio);

        rcp->setLBMKernel(kernel);
        rcp->setBCProcessor(bcProc);
        //////////////////////////////////////////////////////////////////////////
        // BC Adapter
        //////////////////////////////////////////////////////////////////////////////
        mu::Parser fctF1;
        // fctF1.SetExpr("vy1*(1-((x1-x0)^2+(x3-z0)^2)/(R^2))");
        // fctF1.SetExpr("vy1*(1-(sqrt((x1-x0)^2+(x3-z0)^2)/R))^0.1");
        fctF1.SetExpr("vy1");
        fctF1.DefineConst("vy1", 0.0);
        fctF1.DefineConst("R", 8.0);
        fctF1.DefineConst("x0", 0.0);
        fctF1.DefineConst("z0", 0.0);
        //SPtr<BCAdapter> velBCAdapterF1(
        //    new MultiphaseVelocityBCAdapter(false, true, false, fctF1, phiH, 0.0, BCFunction::INFCONST));

        mu::Parser fctF2;
        fctF2.SetExpr("vy1");
        fctF2.DefineConst("vy1", uLB);

        double startTime = 30;
        SPtr<BCAdapter> velBCAdapterF1(new MultiphaseVelocityBCAdapter(true, false, false, fctF1, phiH, 0.0, startTime));
        SPtr<BCAdapter> velBCAdapterF2(new MultiphaseVelocityBCAdapter(true, false, false, fctF2, phiH, startTime, endTime));

        SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
        noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new MultiphaseNoSlipBCAlgorithm()));

        SPtr<BCAdapter> denBCAdapter(new DensityBCAdapter(rhoLB));
        denBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new MultiphaseNonReflectingOutflowBCAlgorithm()));

        mu::Parser fctPhi_F1;
        fctPhi_F1.SetExpr("phiH");
        fctPhi_F1.DefineConst("phiH", phiH);

        mu::Parser fctPhi_F2;
        fctPhi_F2.SetExpr("phiL");
        fctPhi_F2.DefineConst("phiL", phiL);

        mu::Parser fctvel_F2_init;
        fctvel_F2_init.SetExpr("U");
        fctvel_F2_init.DefineConst("U", 0);

        velBCAdapterF1->setBcAlgorithm(SPtr<BCAlgorithm>(new MultiphaseVelocityBCAlgorithm()));
        //////////////////////////////////////////////////////////////////////////////////
        // BC visitor
        MultiphaseBoundaryConditionsBlockVisitor bcVisitor;
        bcVisitor.addBC(noSlipBCAdapter);
        bcVisitor.addBC(denBCAdapter); //Ohne das BB?
        bcVisitor.addBC(velBCAdapterF1);

        SPtr<D3Q27Interactor> inflowF1Int;
        SPtr<D3Q27Interactor> cylInt;
        if (newStart) {

            //  if (newStart) {

            // bounding box
            /*double g_minX1 = 0.0;
            double g_minX2 = -length[1] / 2.0;
            double g_minX3 = -length[2] / 2.0;

            double g_maxX1 = length[0];
            double g_maxX2 = length[1] / 2.0;
            double g_maxX3 = length[2] / 2.0;*/

            double g_minX1 = boundingBox[0];
            double g_minX2 = boundingBox[2];
            double g_minX3 = boundingBox[4];

            double g_maxX1 = boundingBox[1];
            double g_maxX2 = boundingBox[3];
            double g_maxX3 = boundingBox[5];

            // geometry
            SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
            if (myid == 0)
                GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube",
                    WbWriterVtkXmlBinary::getInstance());

            if (myid == 0) UBLOG(logINFO, "Read geoFile:start");
            SPtr<GbTriFaceMesh3D> cylinder = make_shared<GbTriFaceMesh3D>();
            cylinder->readMeshFromSTLFileBinary(pathGeo + "/" + geoFile, false);
            GbSystem3D::writeGeoObject(cylinder.get(), pathname + "/geo/Stlgeo", WbWriterVtkXmlBinary::getInstance());
            if (myid == 0) UBLOG(logINFO, "Read geoFile:stop");
            // inflow
            //GbCuboid3DPtr geoInflowF1(new GbCuboid3D(g_minX1, g_minX2 - 0.5 * dx, g_minX3, g_maxX1, g_minX2 - 1.0 * dx, g_maxX3));
            GbCuboid3DPtr geoInflowF1(
                new GbCuboid3D(g_minX1*0.5 - dx, g_minX2 - dx, g_minX3*0.5 - dx, g_maxX1*0.5 + dx, g_minX2, g_maxX3*0.5 + dx));
            if (myid == 0)  GbSystem3D::writeGeoObject(geoInflowF1.get(), pathname + "/geo/geoInflowF1",                                           WbWriterVtkXmlASCII::getInstance());

            GbCylinder3DPtr cylinder1(new GbCylinder3D(g_minX1-dx, 0.0, 0.0, g_minX1+dx, 0.0, 0.0, 3e-3));
            if (myid == 0)
                GbSystem3D::writeGeoObject(cylinder1.get(), pathname + "/geo/cylinder1",
                    WbWriterVtkXmlASCII::getInstance());

            //GbCylinder3DPtr cylinder2(
            //    new GbCylinder3D(0.0, g_minX2 - 2.0 * dx / 2.0, 0.0, 0.0, g_minX2 + 4.0 * dx, 0.0, 8.0+2.0*dx));
            //if (myid == 0)
            //    GbSystem3D::writeGeoObject(cylinder2.get(), pathname + "/geo/cylinder2",
            //                               WbWriterVtkXmlASCII::getInstance());
            // outflow
            // GbCuboid3DPtr geoOutflow(new GbCuboid3D(-1.0, -1, -1.0, 121.0, 1.0, 121.0)); // For JetBreakup (Original)
            //GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_minX1, g_maxX2 - 40 * dx, g_minX3, g_maxX1, g_maxX2, g_maxX3));
            GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_minX1, g_maxX2, g_minX3, g_maxX1, g_maxX2 + dx, g_maxX3));
            if (myid == 0)
                GbSystem3D::writeGeoObject(geoOutflow.get(), pathname + "/geo/geoOutflow",
                    WbWriterVtkXmlASCII::getInstance());

            // double blockLength = blocknx[0] * dx;

            if (myid == 0) {
                UBLOG(logINFO, "uLb = " << uLB);
                UBLOG(logINFO, "rho = " << rhoLB);
                UBLOG(logINFO, "nuLb = " << nuLB);
                UBLOG(logINFO, "Re = " << Re);
                UBLOG(logINFO, "dx = " << dx);
                UBLOG(logINFO, "Preprocess - start");
            }

            grid->setDeltaX(dx);
            grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

            grid->setPeriodicX1(false);
            grid->setPeriodicX2(false);
            grid->setPeriodicX3(false);

            GenBlocksGridVisitor genBlocks(gridCube);
            grid->accept(genBlocks);

            SPtr<WriteBlocksCoProcessor> ppblocks(new WriteBlocksCoProcessor(
                grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

            SPtr<Interactor3D> tubes(new D3Q27TriFaceMeshInteractor(cylinder, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::POINTS));

            //inflowF1Int =
            //    SPtr<D3Q27Interactor>(new D3Q27Interactor(cylinder1, grid, noSlipBCAdapter, Interactor3D::SOLID));
            //inflowF1Int->addBCAdapter(velBCAdapterF2);

            SPtr<D3Q27Interactor> outflowInt(new D3Q27Interactor(geoOutflow, grid, denBCAdapter, Interactor3D::SOLID));

            // Create boundary conditions geometry
            GbCuboid3DPtr wallXmin(new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_minX3 - dx, g_minX1, g_maxX2 + dx, g_maxX3));
            GbSystem3D::writeGeoObject(wallXmin.get(), pathname + "/geo/wallXmin", WbWriterVtkXmlASCII::getInstance());
            GbCuboid3DPtr wallXmax(new GbCuboid3D(g_maxX1, g_minX2 - dx, g_minX3 - dx, g_maxX1 + dx, g_maxX2 + dx, g_maxX3));
            GbSystem3D::writeGeoObject(wallXmax.get(), pathname + "/geo/wallXmax", WbWriterVtkXmlASCII::getInstance());
            GbCuboid3DPtr wallZmin(new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_minX3 - dx, g_maxX1 + dx, g_maxX2 + dx, g_minX3));
            GbSystem3D::writeGeoObject(wallZmin.get(), pathname + "/geo/wallZmin", WbWriterVtkXmlASCII::getInstance());
            GbCuboid3DPtr wallZmax(new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_maxX3, g_maxX1 + dx, g_maxX2 + dx, g_maxX3 + dx));
            GbSystem3D::writeGeoObject(wallZmax.get(), pathname + "/geo/wallZmax", WbWriterVtkXmlASCII::getInstance());
            GbCuboid3DPtr wallYmin(new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_minX3 - dx, g_maxX1 + dx, g_minX2, g_maxX3));
            GbSystem3D::writeGeoObject(wallYmin.get(), pathname + "/geo/wallYmin", WbWriterVtkXmlASCII::getInstance());
            GbCuboid3DPtr wallYmax(new GbCuboid3D(g_minX1 - dx, g_maxX2, g_minX3 - dx, g_maxX1 + dx, g_maxX2 + dx, g_maxX3));
            GbSystem3D::writeGeoObject(wallYmax.get(), pathname + "/geo/wallYmax", WbWriterVtkXmlASCII::getInstance());

            // Add boundary conditions to grid generator
            SPtr<D3Q27Interactor> wallXminInt(new D3Q27Interactor(wallXmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
            SPtr<D3Q27Interactor> wallXmaxInt(new D3Q27Interactor(wallXmax, grid, noSlipBCAdapter, Interactor3D::SOLID));
            SPtr<D3Q27Interactor> wallZminInt(new D3Q27Interactor(wallZmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
            SPtr<D3Q27Interactor> wallZmaxInt(new D3Q27Interactor(wallZmax, grid, noSlipBCAdapter, Interactor3D::SOLID));
            SPtr<D3Q27Interactor> wallYminInt(new D3Q27Interactor(wallYmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
            SPtr<D3Q27Interactor> wallYmaxInt(new D3Q27Interactor(wallYmax, grid, noSlipBCAdapter, Interactor3D::SOLID));


            cylInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(cylinder1, grid, velBCAdapterF1, Interactor3D::SOLID));
            cylInt->addBCAdapter(velBCAdapterF2);
            //SPtr<D3Q27Interactor> cyl2Int(new D3Q27Interactor(cylinder2, grid, noSlipBCAdapter, Interactor3D::SOLID));


            InteractorsHelper intHelper(grid, metisVisitor, true);
            intHelper.addInteractor(cylInt);
            intHelper.addInteractor(tubes);
            //intHelper.addInteractor(outflowInt);
            //intHelper.addInteractor(cyl2Int);


            intHelper.addInteractor(wallXminInt);
            intHelper.addInteractor(wallXmaxInt);
            intHelper.addInteractor(wallZminInt);
            intHelper.addInteractor(wallZmaxInt);
            intHelper.addInteractor(wallYminInt);
            intHelper.addInteractor(wallYmaxInt);
            //intHelper.addInteractor(inflowF1Int);


            intHelper.selectBlocks();

            ppblocks->process(0);
            ppblocks.reset();

            unsigned long long numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
            int ghostLayer                    = 3;
            unsigned long long numberOfNodesPerBlock =
                (unsigned long long)(blocknx[0]) * (unsigned long long)(blocknx[1]) * (unsigned long long)(blocknx[2]);
            unsigned long long numberOfNodes = numberOfBlocks * numberOfNodesPerBlock;
            unsigned long long numberOfNodesPerBlockWithGhostLayer =
                numberOfBlocks * (blocknx[0] + ghostLayer) * (blocknx[1] + ghostLayer) * (blocknx[2] + ghostLayer);
            double needMemAll =
                double(numberOfNodesPerBlockWithGhostLayer * (27 * sizeof(double) + sizeof(int) + sizeof(float) * 4));
            double needMem = needMemAll / double(comm->getNumberOfProcesses());

            if (myid == 0) {
                UBLOG(logINFO, "Number of blocks = " << numberOfBlocks);
                UBLOG(logINFO, "Number of nodes  = " << numberOfNodes);
                int minInitLevel = grid->getCoarsestInitializedLevel();
                int maxInitLevel = grid->getFinestInitializedLevel();
                for (int level = minInitLevel; level <= maxInitLevel; level++) {
                    int nobl = grid->getNumberOfBlocks(level);
                    UBLOG(logINFO, "Number of blocks for level " << level << " = " << nobl);
                    UBLOG(logINFO, "Number of nodes for level " << level << " = " << nobl * numberOfNodesPerBlock);
                }
                UBLOG(logINFO, "Necessary memory  = " << needMemAll << " bytes");
                UBLOG(logINFO, "Necessary memory per process = " << needMem << " bytes");
                UBLOG(logINFO, "Available memory per process = " << availMem << " bytes");
            }

            MultiphaseSetKernelBlockVisitor kernelVisitor(kernel, nuL, nuG, availMem, needMem);

            grid->accept(kernelVisitor);

            if (refineLevel > 0) {
                SetUndefinedNodesBlockVisitor undefNodesVisitor;
                grid->accept(undefNodesVisitor);
            }

            intHelper.setBC();

            // initialization of distributions
            mu::Parser fct1;
            fct1.SetExpr("phiL");
            fct1.DefineConst("phiL", phiL);
            //MultiphaseInitDistributionsBlockVisitor initVisitor(interfaceThickness);
            MultiphaseVelocityFormInitDistributionsBlockVisitor initVisitor;
            initVisitor.setPhi(fct1);
            grid->accept(initVisitor);
///////////////////////////////////////////////////////////////////////////////////////////
            //{
                // std::vector<std::vector<SPtr<Block3D>>> blockVector;
                // int gridRank = comm->getProcessID();
                // int minInitLevel = grid->getCoarsestInitializedLevel();
                // int maxInitLevel = grid->getFinestInitializedLevel();
                // blockVector.resize(maxInitLevel + 1);
                // for (int level = minInitLevel; level <= maxInitLevel; level++) {
                //    grid->getBlocks(level, gridRank, true, blockVector[level]);
                //}
                //    for (int level = minInitLevel; level <= maxInitLevel; level++) {
                //    for (SPtr<Block3D> block : blockVector[level]) {
                //        if (block) {
                //            int ix1 = block->getX1();
                //            int ix2 = block->getX2();
                //            int ix3 = block->getX3();
                //            int level = block->getLevel();

                //            for (int dir = 0; dir < D3Q27System::ENDDIR; dir++) {
                //                SPtr<Block3D> neighBlock = grid->getNeighborBlock(dir, ix1, ix2, ix3, level);

                //                if (!neighBlock) {

                //                }
                //            }
                //        }
                //    }
                //}
            //    SPtr<Block3D> block = grid->getBlock(0, 0, 0, 0);
            //    SPtr<LBMKernel> kernel = dynamicPointerCast<LBMKernel>(block->getKernel());
            //    SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();

            //    for (int ix3 = 0; ix3 <= 13; ix3++) {
            //        for (int ix2 = 0; ix2 <= 13; ix2++) {
            //            for (int ix1 = 0; ix1 <= 13; ix1++) {
            //                if (ix1 == 0 || ix2 == 0 || ix3 == 0 || ix1 == 13 || ix2 == 13 || ix3 == 13)
            //                    bcArray->setUndefined(ix1, ix2, ix3);
            //            }
            //        }
            //    }
            //}
            ////////////////////////////////////////////////////////////////////////////////////////////
            // boundary conditions grid
            {
                SPtr<UbScheduler> geoSch(new UbScheduler(1));
                SPtr<WriteBoundaryConditionsCoProcessor> ppgeo(new WriteBoundaryConditionsCoProcessor(
                    grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), comm));
                ppgeo->process(0);
                ppgeo.reset();
            }

            if (myid == 0)
                UBLOG(logINFO, "Preprocess - end");
        } else {
            if (myid == 0) {
                UBLOG(logINFO, "Parameters:");
                UBLOG(logINFO, "uLb = " << uLB);
                UBLOG(logINFO, "rho = " << rhoLB);
                UBLOG(logINFO, "nuLb = " << nuLB);
                UBLOG(logINFO, "Re = " << Re);
                UBLOG(logINFO, "dx = " << dx);
                UBLOG(logINFO, "number of levels = " << refineLevel + 1);
                UBLOG(logINFO, "numOfThreads = " << numOfThreads);
                UBLOG(logINFO, "path = " << pathname);
            }

            rcp->restart((int)restartStep);
            grid->setTimeStep(restartStep);

            if (myid == 0)
                UBLOG(logINFO, "Restart - end");
        }

      //  TwoDistributionsSetConnectorsBlockVisitor setConnsVisitor(comm);
      //  grid->accept(setConnsVisitor);

       //ThreeDistributionsSetConnectorsBlockVisitor setConnsVisitor(comm);

        grid->accept(bcVisitor);

        //ThreeDistributionsDoubleGhostLayerSetConnectorsBlockVisitor setConnsVisitor(comm);
        TwoDistributionsDoubleGhostLayerSetConnectorsBlockVisitor setConnsVisitor(comm);
        grid->accept(setConnsVisitor);

        SPtr<UbScheduler> visSch(new UbScheduler(outTime));
        SPtr<WriteMultiphaseQuantitiesCoProcessor> pp(new WriteMultiphaseQuantitiesCoProcessor(
            //SPtr<WriteMacroscopicQuantitiesCoProcessor> pp(new WriteMacroscopicQuantitiesCoProcessor(
            grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
        pp->process(0);

        SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
        SPtr<NUPSCounterCoProcessor> npr(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

        SPtr<UbScheduler> timeBCSch(new UbScheduler(1, startTime, startTime));
        auto timeDepBC = make_shared<TimeDependentBCCoProcessor>(TimeDependentBCCoProcessor(grid, timeBCSch));
        timeDepBC->addInteractor(cylInt);

#ifdef _OPENMP
        omp_set_num_threads(numOfThreads);
#endif

        SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
        SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, endTime));
        calculator->addCoProcessor(npr);
        calculator->addCoProcessor(pp);
        calculator->addCoProcessor(timeDepBC);
        calculator->addCoProcessor(rcp);




        if (myid == 0)
            UBLOG(logINFO, "Simulation-start");
        calculator->calculate();
        if (myid == 0)
            UBLOG(logINFO, "Simulation-end");
    } catch (std::exception &e) {
        cerr << e.what() << endl << flush;
    } catch (std::string &s) {
        cerr << s << endl;
    } catch (...) {
        cerr << "unknown exception" << endl;
    }
}
int main(int argc, char *argv[])
{
    // Sleep(30000);
    if (argv != NULL) {
        if (argv[1] != NULL) {
            run(string(argv[1]));
        } else {
            cout << "Configuration file is missing!" << endl;
        }
    }
}
