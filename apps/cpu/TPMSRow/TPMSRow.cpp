#include <iostream>
#include <string>

//#include <boost/pointer_cast.hpp>

#include "VirtualFluids.h"

using namespace std;
using namespace vf::lbm::dir;
using namespace vf::basics::constant;

void run(string configname)
{
    try {
        vf::basics::ConfigurationFile config;
        config.load(configname);

        string pathname             = config.getValue<string>("pathname");
        int numOfThreads            = config.getValue<int>("numOfThreads");
        vector<int> blocknx         = config.getVector<int>("blocknx");
        double beginTime            = config.getValue<double>("beginTime");
        double endTime              = config.getValue<double>("endTime");
        double outTime              = config.getValue<double>("outTime");
        double availMem             = config.getValue<double>("availMem");
        double nu                   = config.getValue<double>("nu");
        double dx                   = config.getValue<double>("dx");
        double UnitEdgeLength       = config.getValue<double>("UnitEdgeLength");
        double Re                   = config.getValue<double>("Re");
        double Re0                  = config.getValue<double>("Re0");
        //double rhoIn                = config.getValue<double>("rhoIn");
        //string geometry             = config.getValue<string>("geometry");
        vector<double> length       = config.getVector<double>("length");
        //vector<double> FunnelL      = config.getVector<double>("FunnelL");
        //vector<double> FunnelOrigin = config.getVector<double>("FunnelOrigin");
        
        double          timeAvStart       = config.getValue<double>("timeAvStart");
        double          timeAvStop        = config.getValue<double>("timeAvStop");

        vector<double> TPMSL        = config.getVector<double>("TPMSL");
        vector<double> TPMSOrigin   = config.getVector<double>("TPMSOrigin");
        vector<double> gridCubeOrigin = config.getVector<double>("gridCubeOrigin");
        int refineLevel             = config.getValue<int>("refineLevel");
        bool logToFile              = config.getValue<bool>("logToFile");
        double restartStep          = config.getValue<double>("restartStep");
        double cpStart              = config.getValue<double>("cpStart");
        double cpStep               = config.getValue<double>("cpStep");
        bool newStart               = config.getValue<bool>("newStart");

        //SPtr<Communicator> comm = MPICommunicator::getInstance();
        SPtr<vf::mpi::Communicator> comm = vf::mpi::MPICommunicator::getInstance();
        int myid                = comm->getProcessID();
        //int numOfProcesses      = comm->getNumberOfProcesses();

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
        //dx = 1. / 100. / 112.;
        double vx = Re * nu / (UnitEdgeLength / dx);

        SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

        //UbSystem::makeDirectory(pathname);
         //UbSystem::makeDirectory(pathname+ "/mig");
         //UbSystem::makeDirectory(pathname+ "/geo");
         //UbSystem::makeDirectory(pathname+ "/blocks/blocks_");
      

        ////////////////////////////////////////////////////////////////////////
        // BC Adapter
        // BCAdapterPtr gradientAdapter(new VelocityBCAdapter(true, true, true, pdxC, pdyC, pdzC, 0.0,
        // BCFunction::INFCONST));
        // gradientAdapter->setBcAlgorithm(BCAlgorithmPtr(new FluxBCAlgorithm()));
        // BCAdapterPtr cubeNoslipAdapter(new NoSlipBCAdapter(1));
        SPtr<BCAdapter> tpmsNoslipAdapter(new NoSlipBCAdapter());
        //SPtr<BCAdapter> funnelNoslipAdapter(new NoSlipBCAdapter(1));

           // SPtr<BCAdapter> xMinApr(new DensityBCAdapter(0.0000001));
         SPtr<BCAdapter> xMinApr(new DensityBCAdapter());
        //  SPtr<BCAdapter> xMinApr(new VelocityBCAdapter(vx, 0., BCFunction::INFCONST, 0., 0., BCFunction::INFCONST,
         //  0.,0., BCFunction::INFCONST));

        SPtr<BCAdapter> xMaxApr(new DensityBCAdapter(0.));
        //SPtr<BCAdapter> yMinApr(new NoSlipBCAdapter(1));
        //SPtr<BCAdapter> yMaxApr(new NoSlipBCAdapter(1));
        SPtr<BCAdapter> zMinApr(new NoSlipBCAdapter());
        SPtr<BCAdapter> zMaxApr(new NoSlipBCAdapter());

        //SPtr<BCAdapter> zMinFunnelApr(new NoSlipBCAdapter(1));
        //SPtr<BCAdapter> zMaxFunnelApr(new NoSlipBCAdapter(1));

         //tpmsNoslipAdapter->setBcAlgorithm(BCAlgorithmPtr(new NoSlipBCAlgorithm()));
         //tpmsNoslipAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new ThinWallNoSlipBCAlgorithm()));

        tpmsNoslipAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));
        //funnelNoslipAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

         //xMinApr->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));
        // xMinApr->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityBCAlgorithm()));
         xMinApr->setBcAlgorithm(SPtr<BCAlgorithm>(new NonReflectingInflowBCAlgorithm())); 
        // xMinApr->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));
         //xMaxApr->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));
         xMaxApr->setBcAlgorithm(SPtr<BCAlgorithm>(new NonReflectingOutflowBCAlgorithmWithRelaxation()));
        //yMinApr->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));
        //yMaxApr->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));
        zMinApr->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));
        zMaxApr->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

        //zMinFunnelApr->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));
        //zMaxFunnelApr->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

        ////////////////////////////////////////////////////////////////////////
        // BC visitor
        BoundaryConditionsBlockVisitor bcVisitor;
        // bcVisitor.addBC(cubeNoslipAdapter);
        bcVisitor.addBC(tpmsNoslipAdapter);
        //bcVisitor.addBC(funnelNoslipAdapter);
        bcVisitor.addBC(xMinApr);
        bcVisitor.addBC(xMaxApr);
        //bcVisitor.addBC(yMinApr);
        //bcVisitor.addBC(yMaxApr);
        bcVisitor.addBC(zMinApr);
        bcVisitor.addBC(zMaxApr);
        //bcVisitor.addBC(zMinFunnelApr);
        //bcVisitor.addBC(zMaxFunnelApr);

        ////////////////////////////////////////////////////////////////////////    
        //spnonge layer
        //mu::Parser spongeLayer;
        //spongeLayer.SetExpr("x1>=(sizeX-sizeSP)/dx ? (sizeX/dx-(x1+1))/sizeSP/dx/2.0 + 0.5 : 1.0");
        //spongeLayer.DefineConst("sizeX", length[0]);
        //spongeLayer.DefineConst("sizeSP", 0.005);
        //spongeLayer.DefineConst("dx", dx);

        ////////////////////////////////////////////////////////////////////////
        // grid, kernel and BCProcessor
        SPtr<Grid3D> grid(new Grid3D(comm));
        SPtr<LBMKernel> kernel;
        //kernel = SPtr<LBMKernel>(new IncompressibleCumulantLBMKernel());
         kernel = SPtr<LBMKernel>(new CompressibleCumulantLBMKernel());
        //kernel = SPtr<LBMKernel>(new IncompressibleCumulantWithSpongeLayerLBMKernel());       
        //kernel->setWithSpongeLayer(true);
        //kernel->setSpongeLayer(spongeLayer);
        // kernel = ;
         // kernel = SPtr<LBMKernel>(new CumulantK17LBMKernel());
        // 		 mu::Parser fctForcingX1;
        // 		 fctForcingX1.SetExpr("Fx2");
        // 		 fctForcingX1.DefineConst("Fx2", 5e-4);
        // 		 kernel->setForcingX1(fctForcingX1);
        // 		 kernel->setWithForcing(true);
        //
        // SPtr<ThinWallBCProcessor> bcProc(new ThinWallBCProcessor());
        SPtr<BCProcessor> bcProc(new BCProcessor());
        kernel->setBCProcessor(bcProc);

       

        

            SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(
                comm, MetisPartitioningGridVisitor::LevelIntersected, DIR_00M, MetisPartitioner::RECURSIVE));

        //////////////////////////////////////////////////////////////////////////
        // restart
        SPtr<UbScheduler> mSch(new UbScheduler(cpStep, cpStart));
        SPtr<MPIIOMigrationCoProcessor> migCoProcessor(
            new MPIIOMigrationCoProcessor(grid, mSch,metisVisitor, pathname + "/mig", comm));
        migCoProcessor->setLBMKernel(kernel);
        migCoProcessor->setBCProcessor(bcProc);
        //////////////////////////////////////////////////////////////////////////

        if (newStart) {

        
            //GbGyroidThirdOrderPtr tpms;
            // // tpms = GbGyroidThirdOrderPtr(new GbGyroidThirdOrder(0, 0, 0, TPMSL[0], TPMSL[1], TPMSL[2], UnitEdgeLength,
            // // dx));
            // tpms = GbGyroidThirdOrderPtr(new GbGyroidThirdOrder(TPMSOrigin[0], TPMSOrigin[1], TPMSOrigin[2],
            //                                                   TPMSOrigin[0] + TPMSL[0],
            //                                                   TPMSOrigin[1] + TPMSL[1],
            //                                                   TPMSOrigin[2] + TPMSL[2],
            //                                                   UnitEdgeLength, dx, 2.5e-4));

            
            GbGeneralThirdOrderImplicitSurfacePtr tpms;
            tpms = GbGeneralThirdOrderImplicitSurfacePtr(new GbGeneralThirdOrderImplicitSurface(TPMSOrigin[0], TPMSOrigin[1], TPMSOrigin[2],
                                                              TPMSOrigin[0] + TPMSL[0],
                                                              TPMSOrigin[1] + TPMSL[1],
                                                              TPMSOrigin[2] + TPMSL[2],
                                                              UnitEdgeLength, dx, 2.5e-4));

            // 	for (int i = 0; i < 12; i++)
            // 	{
            // 	  cout << tpms->evaluateImplicitFunction(0.002, 0.002, i/1000., 1.)<<endl;
            // 	}

            if (myid == 0)
                GbSystem3D::writeGeoObject(tpms.get(), pathname + "/geo/tpms", WbWriterVtkXmlBinary::getInstance());


            //SPtr<GbTriFaceMesh3D> funnel;
            //SPtr<GbTriFaceMesh3D> funnel(new GbTriFaceMesh3D());
            //funnel->readMeshFromSTLFileBinary(geometry, true);

          

            //funnel = SPtr<GbTriFaceMesh3D>(GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile2(geometry, "tpmsMeshBody", GbTriFaceMesh3D::KDTREE_SAHPLIT, false));
            // funnel->rotate(0.,180,0.);

            //funnel->translate(-funnel->getX1Minimum() - funnel->getLengthX1(),
                              //tpms->getX2Centroid() - funnel->getX2Centroid(),
                              //tpms->getX3Centroid() - funnel->getX3Centroid());
            //if (myid == 0)
                //GbSystem3D::writeGeoObject(funnel.get(), pathname + "/geo/funnel", WbWriterVtkXmlBinary::getInstance());

            double g_minX1 = gridCubeOrigin[0];
            double g_minX2 = gridCubeOrigin[1];
            double g_minX3 = gridCubeOrigin[2];

            double g_maxX1 = gridCubeOrigin[0] + length[0];
            double g_maxX2 = gridCubeOrigin[1] + length[1];
            double g_maxX3 = gridCubeOrigin[2] + length[2];

            SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
            if (myid == 0)
                GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube",
                                           WbWriterVtkXmlBinary::getInstance());

            
            SPtr<GbCuboid3D> spongecube(new GbCuboid3D(TPMSOrigin[0] + TPMSL[0], g_minX2 - dx, g_minX3 - dx,
                                                       g_maxX1 + dx, g_maxX2 + dx, g_maxX3 + dx));
            if (myid == 0)
                GbSystem3D::writeGeoObject(spongecube.get(), pathname + "/geo/spongecube",
                                           WbWriterVtkXmlBinary::getInstance());
            if (myid == 0) {
                // UBLOG(logINFO,"rho = " << rhoLB );
                UBLOG(logINFO, "nu = " << nu);
                UBLOG(logINFO, "Re = " << Re);
                UBLOG(logINFO, "vx = " << vx);
                UBLOG(logINFO, "dx = " << dx);
                UBLOG(logINFO, "Preprocess - start");
            }

            grid->setDeltaX(dx);
            grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
            grid->setPeriodicX1(false);
            grid->setPeriodicX2(true);
            grid->setPeriodicX3(false);

            GenBlocksGridVisitor genBlocks(gridCube);
            grid->accept(genBlocks);

        if (refineLevel > 0)
         {
             //refinement area
            //SPtr<GbObject3D> refCube(new  GbCuboid3D(TPMSOrigin[0],TPMSOrigin[1],TPMSOrigin[2],TPMSL[0],TPMSL[1],TPMSL[2]));
            //if (myid == 0) GbSystem3D::writeGeoObject(refCube.get(), pathname + "/geo/refCube", WbWriterVtkXmlASCII::getInstance());
            
            if (myid == 0) UBLOG(logINFO, "Refinement - start");
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel, comm);
            //refineHelper.addGbObject(sphere, refineLevel);
            refineHelper.addGbObject(tpms, refineLevel);
            refineHelper.refine();
            if (myid == 0) UBLOG(logINFO, "Refinement - end");
         }

            SPtr<CoProcessor> ppblocks(new WriteBlocksCoProcessor(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname,
                                                                  WbWriterVtkXmlBinary::getInstance(), comm));

            ppblocks->process(0);
            
            // GbObject3DPtr solidcube(new GbCuboid3D(0, g_minX2, g_minX3, TPMSL[0], g_maxX2, g_maxX3));
            // if (myid == 0) GbSystem3D::writeGeoObject(solidcube.get(), pathname + "/geo/solidcube",
            // WbWriterVtkXmlBinary::getInstance());

            GbCuboid3DPtr xMin(
                new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_minX3 - dx, g_minX1, g_maxX2 + dx, g_maxX3 + dx));

            /*GbCuboid3DPtr yMin(
                new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_minX3 - dx, g_maxX1, g_minX2, g_maxX3 + dx));
            GbCuboid3DPtr yMax(
                new GbCuboid3D(g_minX1 - dx, g_maxX2, g_minX3 - dx, g_maxX1 + dx, g_maxX2 + dx, g_maxX3 + dx));*/

           /* GbCuboid3DPtr zMinFunnel(
                new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_minX3 - dx, g_maxX1, g_maxX2 + dx, g_minX3));
            GbCuboid3DPtr zMaxFunnel(
                new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_maxX3, g_maxX1 + dx, g_maxX2 + dx, g_maxX3 + dx));*/

            //g_minX1 = 0.;
            // g_minX2 = -length[1] / 2.0;
            // g_minX3 = -length[2] / 2.0;

            //g_maxX1 = TPMSL[0];
            // g_maxX2 = length[1] / 2.0;
            // g_maxX3 -= TPMSL[2] / 2.0;

            GbCuboid3DPtr xMax(new GbCuboid3D(g_maxX1 , g_minX2 - dx, g_minX3 - dx, g_maxX1 + dx, g_maxX2 + dx,
                                              g_maxX3 + dx));

            //GbCuboid3DPtr zMin(new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_minX3 - dx, 1.1 * g_maxX1, g_maxX2 + dx,
            //                                  g_minX3 + 0.5 * (length[2] - TPMSL[2])));
            //GbCuboid3DPtr zMax(new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_maxX3 - 0.5 * (length[2] - TPMSL[2]),
            //                                  1.1 * g_maxX1, g_maxX2 + dx, g_maxX3));

            GbCuboid3DPtr zMin(new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_minX3 - dx, g_maxX1 + dx, g_maxX2 + dx, g_minX3));
            GbCuboid3DPtr zMax(new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_maxX3, g_maxX1 + dx, g_maxX2 + dx, g_maxX3 + dx));

            if (myid == 0)
                GbSystem3D::writeGeoObject(xMin.get(), pathname + "/geo/xMin", WbWriterVtkXmlBinary::getInstance());
            if (myid == 0)
                GbSystem3D::writeGeoObject(xMax.get(), pathname + "/geo/xMax", WbWriterVtkXmlBinary::getInstance());
           /* if (myid == 0)
                GbSystem3D::writeGeoObject(yMin.get(), pathname + "/geo/yMin", WbWriterVtkXmlBinary::getInstance());
            if (myid == 0)
                GbSystem3D::writeGeoObject(yMax.get(), pathname + "/geo/yMax", WbWriterVtkXmlBinary::getInstance());*/
            if (myid == 0)
                GbSystem3D::writeGeoObject(zMin.get(), pathname + "/geo/zMin", WbWriterVtkXmlBinary::getInstance());
            if (myid == 0)
                GbSystem3D::writeGeoObject(zMax.get(), pathname + "/geo/zMax", WbWriterVtkXmlBinary::getInstance());

 /*           if (myid == 0)
                GbSystem3D::writeGeoObject(zMinFunnel.get(), pathname + "/geo/zMinFunnel",
                                           WbWriterVtkXmlBinary::getInstance());
            if (myid == 0)
                GbSystem3D::writeGeoObject(zMaxFunnel.get(), pathname + "/geo/zMaxFunnel",
                                           WbWriterVtkXmlBinary::getInstance());*/

            // D3Q27InteractorPtr cubeInt = D3Q27InteractorPtr(new D3Q27Interactor(solidcube, grid, cubeNoslipAdapter,
            // Interactor3D::SOLID));
            SPtr<D3Q27Interactor> tpmsInt = SPtr<D3Q27Interactor>(
                new D3Q27Interactor(tpms, grid, tpmsNoslipAdapter, Interactor3D::SOLID, Interactor3D::POINTS));
            //SPtr<Interactor3D> funnelInt = SPtr<D3Q27TriFaceMeshInteractor>(
                //new D3Q27TriFaceMeshInteractor(funnel, grid, funnelNoslipAdapter, Interactor3D::SOLID));
            // D3Q27TriFaceMeshInteractorPtr tpmsInt = D3Q27TriFaceMeshInteractorPtr(new
            // D3Q27TriFaceMeshInteractor(tpms, grid, tpmsNoslipAdapter, Interactor3D::SOLID));
            //  tpmsInt->setQs2(0);

            SPtr<D3Q27Interactor> xMinInt = SPtr<D3Q27Interactor>(
                new D3Q27Interactor(xMin, grid, xMinApr, Interactor3D::SOLID, Interactor3D::POINTS));
            SPtr<D3Q27Interactor> xMaxInt = SPtr<D3Q27Interactor>(
                new D3Q27Interactor(xMax, grid, xMaxApr, Interactor3D::SOLID, Interactor3D::POINTS));
          /*  SPtr<D3Q27Interactor> yMinInt =
                SPtr<D3Q27Interactor>(new D3Q27Interactor(yMin, grid, yMinApr, Interactor3D::SOLID));
            SPtr<D3Q27Interactor> yMaxInt =
                SPtr<D3Q27Interactor>(new D3Q27Interactor(yMax, grid, yMaxApr, Interactor3D::SOLID));*/
            SPtr<D3Q27Interactor> zMinInt = SPtr<D3Q27Interactor>(
                new D3Q27Interactor(zMin, grid, zMinApr, Interactor3D::SOLID, Interactor3D::POINTS));
            SPtr<D3Q27Interactor> zMaxInt = SPtr<D3Q27Interactor>(
                new D3Q27Interactor(zMax, grid, zMaxApr, Interactor3D::SOLID, Interactor3D::POINTS));

            /*SPtr<D3Q27Interactor> zMinFunnelInt =
                SPtr<D3Q27Interactor>(new D3Q27Interactor(zMinFunnel, grid, zMinFunnelApr, Interactor3D::SOLID));
            SPtr<D3Q27Interactor> zMaxFunnelInt =
                SPtr<D3Q27Interactor>(new D3Q27Interactor(zMaxFunnel, grid, zMaxFunnelApr, Interactor3D::SOLID));*/

            // return;

            InteractorsHelper intHelper(grid, metisVisitor,false);

            //intHelper.addInteractor(cubeInt);
            //intHelper.addInteractor(zMinFunnelInt);
            //intHelper.addInteractor(zMaxFunnelInt);
            //intHelper.addInteractor(funnelInt);

            intHelper.addInteractor(tpmsInt);
            intHelper.addInteractor(zMinInt);
            intHelper.addInteractor(zMaxInt);

            intHelper.addInteractor(xMinInt);
            intHelper.addInteractor(xMaxInt);
            //intHelper.addInteractor(yMinInt);
            //intHelper.addInteractor(yMaxInt);


            intHelper.selectBlocks();
            // intHelper.selectBlocks2();

            
            // domain decomposition for threads
            PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
            grid->accept(pqPartVisitor);

            ppblocks->process(0);
            ppblocks.reset();

            //////////////////////////////////////////////////////////////////////////
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
            //////////////////////////////////////////////////////////////////////////

            SetKernelBlockVisitor kernelVisitor(kernel, nu, availMem, needMem);
            grid->accept(kernelVisitor);

            //          if (refineLevel > 0)
            //          {
            // 			 SetUndefinedNodesBlockVisitor undefNodesVisitor;
            //             grid->accept(undefNodesVisitor);
            //          }

            intHelper.setBC();

            SpongeLayerBlockVisitor spongeLayerVisitor(spongecube, kernel, nu, DIR_P00);
            grid->accept(spongeLayerVisitor);

            grid->accept(bcVisitor);

            // initialization of distributions
            InitDistributionsBlockVisitor initVisitor;
             //initVisitor.setVx1(0.001);
            // initVisitor.setVx1(uLB);
            grid->accept(initVisitor);

            // boundary conditions grid
            {
                SPtr<UbScheduler> geoSch(new UbScheduler(1));
                SPtr<CoProcessor> ppgeo(new WriteBoundaryConditionsCoProcessor(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), comm));
                ppgeo->process(0);
                ppgeo.reset();
            }
            if (myid == 0)
                UBLOG(logINFO, "Preprocess - end");
        } 
        else 
        {
            if (myid == 0) {
                UBLOG(logINFO, "Parameters:");
                //UBLOG(logINFO, "uLb = " << uLB);
                //UBLOG(logINFO, "rho = " << rhoLB);
                //UBLOG(logINFO, "nuLb = " << nuLB);
                UBLOG(logINFO, "Re = " << Re);
                UBLOG(logINFO, "dx = " << dx);
                UBLOG(logINFO, "number of levels = " << refineLevel + 1);
                UBLOG(logINFO, "numOfThreads = " << numOfThreads);
                UBLOG(logINFO, "path = " << pathname);
            }

            migCoProcessor->restart((int)restartStep);
            grid->setTimeStep(restartStep);

            if (myid == 0)
                UBLOG(logINFO, "Restart - end");
        }
        // set connectors
        SPtr<InterpolationProcessor> iProcessor(new CompressibleOffsetInterpolationProcessor());
        //SetConnectorsBlockVisitor setConnsVisitor(comm, true, D3Q27System::ENDDIR, nu, iProcessor);
        OneDistributionSetConnectorsBlockVisitor setConnsVisitor(comm);
        grid->accept(setConnsVisitor);


        

        SPtr<UbScheduler> visSch(new UbScheduler(outTime/*,beginTime,endTime*/));
        SPtr<CoProcessor> pp(new WriteMacroscopicQuantitiesCoProcessor(grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
        
        SPtr<UbScheduler> tavSch(new UbScheduler(100, timeAvStart, timeAvStop));
        SPtr<TimeAveragedValuesCoProcessor> tav(new TimeAveragedValuesCoProcessor(grid, pathname, WbWriterVtkXmlBinary::getInstance(), tavSch, comm,
        TimeAveragedValuesCoProcessor::Density | TimeAveragedValuesCoProcessor::Velocity | TimeAveragedValuesCoProcessor::Fluctuations));
        tav->setWithGhostLayer(true);        
        
        SPtr<UbScheduler> nuSch(new UbScheduler(100, 0, endTime / 2));
        mu::Parser fnu;
        fnu.SetExpr("(L*u/T)*(((T-2*t)/Re0)+(2*t/Re))");
        fnu.DefineConst("Re0", Re0);
        fnu.DefineConst("Re", Re);
        fnu.DefineConst("T", endTime);
        fnu.DefineConst("L", (UnitEdgeLength / dx));
        fnu.DefineConst("u", vx);
        SPtr<CoProcessor> nupr(new DecreaseViscosityCoProcessor(grid, nuSch, &fnu, comm));

        SPtr<UbScheduler> nupsSch(new UbScheduler(10, 10, 100000000));
        SPtr<CoProcessor> npr(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));

        //omp_set_num_threads(numOfThreads);
        numOfThreads = 1;
        SPtr<UbScheduler> stepGhostLayer(visSch);
        SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, int(endTime)));

        //calculator->addCoProcessor(nupr);
        calculator->addCoProcessor(npr);
        calculator->addCoProcessor(pp);
        calculator->addCoProcessor(migCoProcessor);
        calculator->addCoProcessor(tav);

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
     //Sleep(25000);
    if (argv != NULL) {
        if (argv[1] != NULL) {
            run(string(argv[1]));
        } else {
            cout << "Configuration file is missing!" << endl;
        }
    }
}
