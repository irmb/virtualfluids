#include <iostream>
#include <string>

#include "VirtualFluids.h"

using namespace std;

void run(string configname)
{
    try {
        ConfigurationFile config;
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
        int interfaceThickness = config.getValue<int>("interfaceThickness");
        double radius          = config.getValue<double>("radius");
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

        double beta  = 12 * sigma / interfaceThickness;
        double kappa = 1.5 * interfaceThickness * sigma;

        SPtr<Communicator> comm = MPICommunicator::getInstance();
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

        kernel = SPtr<LBMKernel>(new MultiphaseScratchCumulantLBMKernel());
  //      kernel = SPtr<LBMKernel>(new MultiphaseCumulantLBMKernel());

        kernel->setWithForcing(true);
        kernel->setForcingX1(0.0);
        kernel->setForcingX2(gr);
        kernel->setForcingX3(0.0);

        kernel->setPhiL(phiL);
        kernel->setPhiH(phiH);
        kernel->setPhaseFieldRelaxation(tauH);
        kernel->setMobility(mob);

        SPtr<BCProcessor> bcProc(new BCProcessor());
        // BCProcessorPtr bcProc(new ThinWallBCProcessor());

        kernel->setBCProcessor(bcProc);

        SPtr<Grid3D> grid(new Grid3D(comm));
        // grid->setPeriodicX1(true);
        // grid->setPeriodicX2(true);
        // grid->setPeriodicX3(true);
        //////////////////////////////////////////////////////////////////////////
        // restart
        SPtr<UbScheduler> rSch(new UbScheduler(cpStep, cpStart));
        // RestartCoProcessor rp(grid, rSch, comm, pathname, RestartCoProcessor::TXT);
        MPIIORestartCoProcessor rcp(grid, rSch, pathname, comm);
        rcp.setLBMKernel(kernel);
        rcp.setBCProcessor(bcProc);
        //////////////////////////////////////////////////////////////////////////

        mu::Parser fctF1;
        // fctF1.SetExpr("vy1*(1-((x1-x0)^2+(x3-z0)^2)/(R^2))");
        // fctF1.SetExpr("vy1*(1-(sqrt((x1-x0)^2+(x3-z0)^2)/R))^0.1");
        fctF1.SetExpr("vy1");
        fctF1.DefineConst("vy1", -uLB*0);
        //fctF1.DefineConst("vy1", -uLB);
        fctF1.DefineConst("R", 8.0);
        fctF1.DefineConst("x0", 0.0);
        fctF1.DefineConst("z0", 0.0);

      //  mu::Parser fctFStart;
      //  fctFStart.SetExpr("0");

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

            // GbObject3DPtr innerCube(new GbCuboid3D(g_minX1+2, g_minX2+2, g_minX3+2, g_maxX1-2, g_maxX2-2,
            // g_maxX3-2));

            // GbObject3DPtr cylinder1(new GbCylinder3D(g_minX1 - 2.0*dx, g_maxX2/2, g_maxX3/2, g_minX1 + 12.0*dx,
            // g_maxX2/2, g_maxX3/2, radius)); GbObject3DPtr cylinder2(new GbCylinder3D(g_minX1 + 12.0*dx, g_maxX2/2,
            // g_maxX3/2, g_maxX1 + 2.0*dx, g_maxX2/2, g_maxX3/2, dLB / 2.0));

            // GbObject3DPtr cylinder(new GbCylinder3D(g_minX1 - 2.0*dx, g_maxX2/2, g_maxX3/2, g_maxX1 + 2.0*dx,
            // g_maxX2/2, g_maxX3/2, dLB / 2.0)); GbObject3DPtr cylinders(new GbObject3DManager()); GbObject3DPtr
            // cylinders1(new GbObjectGroup3D());

            SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
            if (myid == 0)
                GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube",
                                           WbWriterVtkXmlBinary::getInstance());


            if (myid == 0) UBLOG(logINFO, "Read geoFile:start");
            SPtr<GbTriFaceMesh3D> cylinder = make_shared<GbTriFaceMesh3D>();
            cylinder->readMeshFromSTLFileASCII(pathGeo + "/" + geoFile, false);
            GbSystem3D::writeGeoObject(cylinder.get(), pathname + "/geo/Stlgeo", WbWriterVtkXmlBinary::getInstance());
            if (myid == 0) UBLOG(logINFO, "Read geoFile:stop");
            // inflow
            // GbCuboid3DPtr geoInflowF1(new GbCuboid3D(40.0, 628.0, 40.0, 80, 631.0, 80.0));  // For JetBreakup
            // (Original) GbCuboid3DPtr geoInflowF1(new GbCuboid3D(g_minX1-2.0*dx, g_minX2-2.0*dx, g_minX3-2.0*dx,
            // g_maxX1+2.0*dx, g_minX2+2.0*dx, g_maxX3+2.0*dx)); if (myid == 0)
            // GbSystem3D::writeGeoObject(geoInflowF1.get(), pathname + "/geo/geoInflowF1",
            // WbWriterVtkXmlASCII::getInstance());

            ////outflow
            ////GbCuboid3DPtr geoOutflow(new GbCuboid3D(-1.0, -1, -1.0, 121.0, 1.0, 121.0)); // For JetBreakup
            ///(Original)
            // GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_minX1-2.0*dx, g_maxX2, g_minX3-2.0*dx, g_maxX1+2.0*dx,
            // g_maxX2+2.0*dx, g_maxX3+2.0*dx)); if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname +
            // "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

            GbCuboid3DPtr geoInflowF1(
                new GbCuboid3D(g_minX1, g_minX2 - 0.5 * dx, g_minX3, g_maxX1, g_minX2 - 1.0 * dx, g_maxX3));
            if (myid == 0)
                GbSystem3D::writeGeoObject(geoInflowF1.get(), pathname + "/geo/geoInflowF1",
                                           WbWriterVtkXmlASCII::getInstance());

            // outflow
            // GbCuboid3DPtr geoOutflow(new GbCuboid3D(-1.0, -1, -1.0, 121.0, 1.0, 121.0)); // For JetBreakup (Original)
            GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_minX1, g_maxX2 - 1 * dx, g_minX3, g_maxX1, g_maxX2, g_maxX3));
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

            // BC Adapter
            //////////////////////////////////////////////////////////////////////////////
            SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
            noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new MultiphaseNoSlipBCAlgorithm()));

            SPtr<BCAdapter> denBCAdapter(new DensityBCAdapter(rhoLB));
            denBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new MultiphaseNonReflectingOutflowBCAlgorithm()));

            // double r = 5.0; //boost::dynamic_pointer_cast<GbCylinder3D>(cylinder)->getRadius();
            // double cx1 = g_minX1;
            // double cx2 = 0.0; //cylinder->getX2Centroid();
            // double cx3 = 0.0; //cylinder->getX3Centroid();

            mu::Parser fctPhi_F1;
            fctPhi_F1.SetExpr("phiH");
            fctPhi_F1.DefineConst("phiH", phiH);

            mu::Parser fctPhi_F2;
            fctPhi_F2.SetExpr("phiL");
            fctPhi_F2.DefineConst("phiL", phiL);

            mu::Parser fctvel_F2_init;
            fctvel_F2_init.SetExpr("U");
            fctvel_F2_init.DefineConst("U", 0);

            // fct.SetExpr("U");
            // fct.DefineConst("U", uLB);
            // BCAdapterPtr velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));

            //Am Anfang keine Geschwindigkeit
     //       double startGeschwindigkeit = 200;
     //       SPtr<BCAdapter> velBCAdapterFStart(
     //           new MultiphaseVelocityBCAdapter(false, true, false, fctFStart, phiH, 0.0, startGeschwindigkeit));


            SPtr<BCAdapter> velBCAdapterF1(
                new MultiphaseVelocityBCAdapter(false, true, false, fctF1, phiH, 0, endTime));

            // BCAdapterPtr velBCAdapterF2_1_init(new MultiphaseVelocityBCAdapter(false, false, true, fctF2_1, phiH,
            // 0.0, endTime)); BCAdapterPtr velBCAdapterF2_2_init(new MultiphaseVelocityBCAdapter(false, false, true,
            // fctF2_2, phiH, 0.0, endTime));

            // BCAdapterPtr velBCAdapterF2_1_init(new MultiphaseVelocityBCAdapter(false, false, true, fctvel_F2_init,
            // phiL, 0.0, endTime)); BCAdapterPtr velBCAdapterF2_2_init(new MultiphaseVelocityBCAdapter(false, false,
            // true, fctvel_F2_init, phiL, 0.0, endTime));

            velBCAdapterF1->setBcAlgorithm(SPtr<BCAlgorithm>(new MultiphaseVelocityBCAlgorithm()));
           // velBCAdapterFStart->setBcAlgorithm(SPtr<BCAlgorithm>(new MultiphaseVelocityBCAlgorithm()));
            // velBCAdapterF2_1_init->setBcAlgorithm(BCAlgorithmPtr(new MultiphaseVelocityBCAlgorithm()));
            // velBCAdapterF2_2_init->setBcAlgorithm(BCAlgorithmPtr(new MultiphaseVelocityBCAlgorithm()));

            // velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new VelocityWithDensityBCAlgorithm()));
            // mu::Parser fct;
            // fct.SetExpr("U");
            // fct.DefineConst("U", uLB);
            // BCAdapterPtr velBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
            // velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new NonReflectingVelocityBCAlgorithm()));

            //////////////////////////////////////////////////////////////////////////////////
            // BC visitor
            MultiphaseBoundaryConditionsBlockVisitor bcVisitor;
            bcVisitor.addBC(noSlipBCAdapter);
           // bcVisitor.addBC(denBCAdapter); //Ohne das BB?
            bcVisitor.addBC(velBCAdapterF1);
            // bcVisitor.addBC(velBCAdapterF2_1_init);
            // bcVisitor.addBC(velBCAdapterF2_2_init);

            SPtr<WriteBlocksCoProcessor> ppblocks(new WriteBlocksCoProcessor(
                grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

            //ppblocks->process(0);

            SPtr<Interactor3D> tubes(
                new D3Q27TriFaceMeshInteractor(cylinder, grid, noSlipBCAdapter, Interactor3D::SOLID));

            SPtr<D3Q27Interactor> inflowF1Int(
                new D3Q27Interactor(geoInflowF1, grid, velBCAdapterF1, Interactor3D::SOLID));

           // inflowF1Int->addBCAdapter(velBCAdapterFStart);

            // D3Q27InteractorPtr inflowF2_1Int_init = D3Q27InteractorPtr(new D3Q27Interactor(geoInflowF2_1, grid,
            // velBCAdapterF2_1_init, Interactor3D::SOLID));

            // D3Q27InteractorPtr inflowF2_2Int_init = D3Q27InteractorPtr(new D3Q27Interactor(geoInflowF2_2, grid,
            // velBCAdapterF2_2_init, Interactor3D::SOLID));

            SPtr<D3Q27Interactor> outflowInt(new D3Q27Interactor(geoOutflow, grid, denBCAdapter, Interactor3D::SOLID));

            // SetSolidBlockVisitor visitor1(inflowF2_1Int, SetSolidBlockVisitor::BC);
            // grid->accept(visitor1);
            // SetSolidBlockVisitor visitor2(inflowF2_2Int, SetSolidBlockVisitor::BC);
            // grid->accept(visitor2);

            SPtr<Grid3DVisitor> metisVisitor(
                new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW));
            InteractorsHelper intHelper(grid, metisVisitor);
            intHelper.addInteractor(tubes);
            intHelper.addInteractor(inflowF1Int);
            intHelper.addInteractor(outflowInt);
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

            // LBMKernelPtr kernel;

            // kernel = LBMKernelPtr(new MultiphaseCumulantLBMKernel(blocknx[0], blocknx[1], blocknx[2],
            // MultiphaseCumulantLBMKernel::NORMAL));

            // kernel->setWithForcing(true);
            // kernel->setForcingX1(0.0);
            // kernel->setForcingX2(gr);
            // kernel->setForcingX3(0.0);

            // kernel->setPhiL(phiL);
            // kernel->setPhiH(phiH);
            // kernel->setPhaseFieldRelaxation(tauH);
            // kernel->setMobility(mob);

            // BCProcessorPtr bcProc(new BCProcessor());
            // //BCProcessorPtr bcProc(new ThinWallBCProcessor());

            // kernel->setBCProcessor(bcProc);

            MultiphaseSetKernelBlockVisitor kernelVisitor(kernel, nuL, nuG, densityRatio, beta, kappa, theta, availMem,
                                                          needMem);

            grid->accept(kernelVisitor);

            if (refineLevel > 0) {
                SetUndefinedNodesBlockVisitor undefNodesVisitor;
                grid->accept(undefNodesVisitor);
            }

            // inflowF2_1Int->initInteractor();
            // inflowF2_2Int->initInteractor();

            intHelper.setBC();

            grid->accept(bcVisitor);

            // initialization of distributions
            LBMReal x1c = radius;                  // g_minX1; //radius; //19; //(g_maxX1+g_minX1)/2;
            LBMReal x2c = (g_maxX2 + g_minX2) / 2; // g_minX2 + 2;
            LBMReal x3c = (g_maxX3 + g_minX3) / 2;
            mu::Parser fct1;

            // fct1.SetExpr("0.5-0.5*tanh(2*(sqrt((x1-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)");
            // fct1.SetExpr("phiM-phiM*tanh((sqrt((x1-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/(interfaceThickness*phiM))");

            // fct1.SetExpr("0.5*(phiH + phiL)-0.5*(phiH -
            // phiL)*tanh(2*(sqrt((x1-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)");

            // fct1.SetExpr("0.5*(phiH + phiL) + 0.5*(phiH - phiL)*tanh(2*((x2-radius))/interfaceThickness)");
            fct1.SetExpr("phiL");
            fct1.DefineConst("x1c", x1c);
            fct1.DefineConst("x2c", x2c);
            fct1.DefineConst("x3c", x3c);
            fct1.DefineConst("phiL", phiL);
            fct1.DefineConst("phiH", phiH);
            fct1.DefineConst("radius", radius);
            fct1.DefineConst("interfaceThickness", interfaceThickness);

            mu::Parser fct2;
            // fct2.SetExpr("vx1*(1-((x2-y0)^2+(x3-z0)^2)/(R^2))");
            fct2.SetExpr("vx1");
            fct2.DefineConst("R", 10.0);
            fct2.DefineConst("vx1", uLB);
            fct2.DefineConst("y0", 1.0);
            fct2.DefineConst("z0", 31.0);
            /*fct2.SetExpr("0.5*uLB-uLB*0.5*tanh(2*(sqrt((x1-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)");
            fct2.DefineConst("uLB", uLB);
            fct2.DefineConst("x1c", x1c);
            fct2.DefineConst("x2c", x2c);
            fct2.DefineConst("x3c", x3c);
            fct2.DefineConst("radius", radius);
            fct2.DefineConst("interfaceThickness", interfaceThickness);*/

            MultiphaseInitDistributionsBlockVisitor initVisitor(densityRatio, interfaceThickness, radius);
            initVisitor.setPhi(fct1);
            // initVisitor.setVx1(fct2);
            grid->accept(initVisitor);

            // set connectors
            //InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
            // InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
            //SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
            // ConnectorFactoryPtr factory(new Block3DConnectorFactory());
            // ConnectorBlockVisitor setConnsVisitor(comm, nuLB, iProcessor, factory);
           // TwoDistributionsSetConnectorsBlockVisitor setConnsVisitor(comm);
           // grid->accept(setConnsVisitor);

            // domain decomposition for threads
            // PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
            // grid->accept(pqPartVisitor);

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
        //} else {
        //    if (myid == 0) {
        //        UBLOG(logINFO, "Parameters:");
        //        UBLOG(logINFO, "uLb = " << uLB);
        //        UBLOG(logINFO, "rho = " << rhoLB);
        //        UBLOG(logINFO, "nuLb = " << nuLB);
        //        UBLOG(logINFO, "Re = " << Re);
        //        UBLOG(logINFO, "dx = " << dx);
        //        UBLOG(logINFO, "number of levels = " << refineLevel + 1);
        //        UBLOG(logINFO, "numOfThreads = " << numOfThreads);
        //        UBLOG(logINFO, "path = " << pathname);
        //    }

        //    rcp.restart((int)restartStep);
        //    grid->setTimeStep(restartStep);

        //    // BCAdapterPtr velBCAdapter(new VelocityBCAdapter());
        //    // velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new VelocityBCAlgorithm()));
        //    // velBCAdapter->setBcAlgorithm(BCAlgorithmPtr(new VelocityWithDensityBCAlgorithm()));
        //    // bcVisitor.addBC(velBCAdapter);
        //    // grid->accept(bcVisitor);

        //    // set connectors
        //    // InterpolationProcessorPtr iProcessor(new IncompressibleOffsetInterpolationProcessor());
        //    //InterpolationProcessorPtr iProcessor(new CompressibleOffsetInterpolationProcessor());
        //    //SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nuLB, iProcessor);
        //    //grid->accept(setConnsVisitor);

        //    if (myid == 0)
        //        UBLOG(logINFO, "Restart - end");
        //}
        
        TwoDistributionsSetConnectorsBlockVisitor setConnsVisitor(comm);
        grid->accept(setConnsVisitor);

        SPtr<UbScheduler> visSch(new UbScheduler(outTime));
        SPtr<WriteMultiphaseQuantitiesCoProcessor> pp(new WriteMultiphaseQuantitiesCoProcessor(
            grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));

        SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
        SPtr<NUPSCounterCoProcessor> npr(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

//////
        double fluxStart = 500.0;
        SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
        SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, fluxStart));
        calculator->addCoProcessor(npr);
        calculator->addCoProcessor(pp);

        calculator->calculate();
        dynamicPointerCast<MultiphaseVelocityBCAdapter>(velBCAdapterF1)->setNewVelocities(0.0,0.0,endTime, -uLB,0.0,endTime,0.0,0.0,endTime);
        inflowF1Int->initInteractor();

       // SPtr<D3Q27Interactor> inflowF1Int(
       //     new D3Q27Interactor(geoInflowF1, grid, velBCAdapterF1, Interactor3D::SOLID));

//////



       // SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
        grid->setTimeStep(fluxStart);
        SPtr<Calculator> calculator2(new BasicCalculator(grid, stepGhostLayer, endTime));
        calculator2->addCoProcessor(npr);
        calculator2->addCoProcessor(pp);


      


        if (myid == 0)
            UBLOG(logINFO, "Simulation-start");
        calculator2->calculate();
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
