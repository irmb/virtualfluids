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
        double endTime              = config.getValue<double>("endTime");
        double outTime              = config.getValue<double>("outTime");
        double availMem             = config.getValue<double>("availMem");
        double nu                   = config.getValue<double>("nu");
        double dx                   = config.getValue<double>("dx");
        double UnitEdgeLength       = config.getValue<double>("UnitEdgeLength");
        double Re                   = config.getValue<double>("Re");
        double vx                   = config.getValue<double>("vx");
        vector<double> length       = config.getVector<double>("length"); 
        //double          timeAvStart       = config.getValue<double>("timeAvStart");
        //double          timeAvStop        = config.getValue<double>("timeAvStop");
        vector<double> TPMSL        = config.getVector<double>("TPMSL");
        vector<double> TPMSOrigin   = config.getVector<double>("TPMSOrigin");
        vector<double> gridCubeOrigin = config.getVector<double>("gridCubeOrigin");
        int refineLevel             = config.getValue<int>("refineLevel");
        bool logToFile              = config.getValue<bool>("logToFile");
        double restartStep          = config.getValue<double>("restartStep");
        double cpStart              = config.getValue<double>("cpStart");
        double cpStep               = config.getValue<double>("cpStep");
        bool newStart               = config.getValue<bool>("newStart");

        SPtr<vf::parallel::Communicator> comm = vf::parallel::MPICommunicator::getInstance();
        int myid                = comm->getProcessID();

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


        SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());



        ////////////////////////////////////////////////////////////////////////
        // BC Adapter
        SPtr<BC> tpmsNoslipAdapter(new NoSlipBC());
        SPtr<BC> xMinApr(new VelocityBC(vx, 0., BCFunction::INFCONST, 0., 0., BCFunction::INFCONST, 0.,0., BCFunction::INFCONST));
        SPtr<BC> xMaxApr(new PressureBC(0.));
        SPtr<BC> zMinApr(new NoSlipBC());
        SPtr<BC> zMaxApr(new NoSlipBC());

        tpmsNoslipAdapter->setBCStrategy(SPtr<BCStrategy>(new NoSlipInterpolated()));
        xMinApr->setBCStrategy(SPtr<BCStrategy>(new VelocityNonReflecting(c1o2))); 
        xMaxApr->setBCStrategy(SPtr<BCStrategy>(new OutflowNonReflectingWithPressure(c1o100)));
        zMinApr->setBCStrategy(SPtr<BCStrategy>(new NoSlipInterpolated()));
        zMaxApr->setBCStrategy(SPtr<BCStrategy>(new NoSlipInterpolated()));

        ////////////////////////////////////////////////////////////////////////234
        // BC visitor
        BoundaryConditionsBlockVisitor bcVisitor;
        ////////////////////////////////////////////////////////////////////////
        // grid, kernel and BCProcessor
        SPtr<Grid3D> grid(new Grid3D(comm));
        SPtr<LBMKernel> kernel;
         kernel = SPtr<LBMKernel>(new K15CompressibleNavierStokes());;
         SPtr<BCSet> bcProc(new BCSet());
         kernel->setBCSet(bcProc);
        SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelIntersected, d00M, MetisPartitioner::RECURSIVE));

        //////////////////////////////////////////////////////////////////////////
        // restart
        SPtr<UbScheduler> mSch(new UbScheduler(cpStep, cpStart));
        SPtr<MPIIOMigrationSimulationObserver> migSimulationObserver(
            new MPIIOMigrationSimulationObserver(grid, mSch,metisVisitor, pathname + "/mig", comm));
        migSimulationObserver->setLBMKernel(kernel);
        migSimulationObserver->setBCSet(bcProc);
        //////////////////////////////////////////////////////////////////////////

        if (newStart) {
            //GbGyroidThirdOrderPtr tpms;            
            // tpms = GbGyroidThirdOrderPtr(new GbGyroidThirdOrder(TPMSOrigin[0], TPMSOrigin[1], TPMSOrigin[2],
            //                                                   TPMSOrigin[0] + TPMSL[0],
            //                                                   TPMSOrigin[1] + TPMSL[1],
            //                                                   TPMSOrigin[2] + TPMSL[2],
            //                                                   UnitEdgeLength, dx, 2.5e-4));
            GbGyroidThirdOrderLongPtr tpms;
            tpms = GbGyroidThirdOrderLongPtr(new GbGyroidThirdOrderLong(TPMSOrigin[0], TPMSOrigin[1], TPMSOrigin[2],
                                                              TPMSOrigin[0] + TPMSL[0],
                                                              TPMSOrigin[1] + TPMSL[1],
                                                              TPMSOrigin[2] + TPMSL[2],
                                                              UnitEdgeLength, dx, 2.5e-4));

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

            SPtr<SimulationObserver> ppblocks(new WriteBlocksSimulationObserver(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname,
                                                                  WbWriterVtkXmlBinary::getInstance(), comm));

            ppblocks->update(0);


            GbCuboid3DPtr xMin( new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_minX3 - dx, g_minX1, g_maxX2 + dx, g_maxX3 + dx));
            GbCuboid3DPtr xMax(new GbCuboid3D(g_maxX1 , g_minX2 - dx, g_minX3 - dx, g_maxX1 + dx, g_maxX2 + dx, g_maxX3 + dx));
            GbCuboid3DPtr zMin(new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_minX3 - dx, g_maxX1 + dx, g_maxX2 + dx, g_minX3));
            GbCuboid3DPtr zMax(new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_maxX3, g_maxX1 + dx, g_maxX2 + dx, g_maxX3 + dx));

            if (myid == 0)
                GbSystem3D::writeGeoObject(xMin.get(), pathname + "/geo/xMin", WbWriterVtkXmlBinary::getInstance());
            if (myid == 0)
                GbSystem3D::writeGeoObject(xMax.get(), pathname + "/geo/xMax", WbWriterVtkXmlBinary::getInstance());
            if (myid == 0)
                GbSystem3D::writeGeoObject(zMin.get(), pathname + "/geo/zMin", WbWriterVtkXmlBinary::getInstance());
            if (myid == 0)
                GbSystem3D::writeGeoObject(zMax.get(), pathname + "/geo/zMax", WbWriterVtkXmlBinary::getInstance());

            SPtr<D3Q27Interactor> tpmsInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(tpms, grid, tpmsNoslipAdapter, Interactor3D::SOLID, Interactor3D::POINTS));
            SPtr<D3Q27Interactor> xMinInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(xMin, grid, xMinApr, Interactor3D::SOLID, Interactor3D::POINTS));
            SPtr<D3Q27Interactor> xMaxInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(xMax, grid, xMaxApr, Interactor3D::SOLID, Interactor3D::POINTS));
            SPtr<D3Q27Interactor> zMinInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(zMin, grid, zMinApr, Interactor3D::SOLID, Interactor3D::POINTS));
            SPtr<D3Q27Interactor> zMaxInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(zMax, grid, zMaxApr, Interactor3D::SOLID, Interactor3D::POINTS));

            InteractorsHelper intHelper(grid, metisVisitor,false);

            intHelper.addInteractor(tpmsInt);
            intHelper.addInteractor(zMinInt);
            intHelper.addInteractor(zMaxInt);

            intHelper.addInteractor(xMinInt);
            intHelper.addInteractor(xMaxInt);



            intHelper.selectBlocks();
            
            // domain decomposition for threads
            //PQueuePartitioningGridVisitor pqPartVisitor(numOfThreads);
            //grid->accept(pqPartVisitor);

            ppblocks->update(0);
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

            SetKernelBlockVisitor kernelVisitor(kernel, nu);
            grid->accept(kernelVisitor);

            intHelper.setBC();

            SpongeLayerBlockVisitor spongeLayerVisitor(spongecube, kernel, nu, dP00);
            grid->accept(spongeLayerVisitor);

            grid->accept(bcVisitor);

            // initialization of distributions
            InitDistributionsBlockVisitor initVisitor;
            grid->accept(initVisitor);

            // boundary conditions grid
            {
                SPtr<UbScheduler> geoSch(new UbScheduler(1));
                SPtr<SimulationObserver> ppgeo(new WriteBoundaryConditionsSimulationObserver(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), comm));
                ppgeo->update(0);
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

            migSimulationObserver->restart((int)restartStep);
            grid->setTimeStep(restartStep);

            if (myid == 0)
                UBLOG(logINFO, "Restart - end");
        }
        // set connectors
        SPtr<Interpolator> iProcessor(new CompressibleOffsetMomentsInterpolator());
        OneDistributionSetConnectorsBlockVisitor setConnsVisitor(comm);
        grid->accept(setConnsVisitor);

       
        SPtr<UbScheduler> visSch(new UbScheduler(outTime/*,beginTime,endTime*/));
        SPtr<SimulationObserver> pp(new WriteMacroscopicQuantitiesSimulationObserver(grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
        
        //SPtr<UbScheduler> tavSch(new UbScheduler(100, timeAvStart, timeAvStop));
        //SPtr<TimeAveragedValuesSimulationObserver> tav(new TimeAveragedValuesSimulationObserver(grid, pathname, WbWriterVtkXmlBinary::getInstance(), tavSch, comm,
        //TimeAveragedValuesSimulationObserver::Density | TimeAveragedValuesSimulationObserver::Velocity | TimeAveragedValuesSimulationObserver::Fluctuations));
        //tav->setWithGhostLayer(true);        
        

        SPtr<UbScheduler> nupsSch(new UbScheduler(500, 1000, 3000));
        //OpenMP threads control
#ifdef _OPENMP
      omp_set_num_threads(numOfThreads);
#endif

        SPtr<SimulationObserver> npr(new NUPSCounterSimulationObserver(grid, nupsSch, numOfThreads, comm));


        SPtr<UbScheduler> stepGhostLayer(visSch);
        SPtr<Simulation> calculator(new Simulation(grid, stepGhostLayer, int(endTime)));

        //calculator->addSimulationObserver(nupr);
        calculator->addSimulationObserver(npr);
        calculator->addSimulationObserver(pp);
        calculator->addSimulationObserver(migSimulationObserver);
        //calculator->addSimulationObserver(tav);

        if (myid == 0)
            UBLOG(logINFO, "Simulation-start");
        calculator->run();
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
    if (argv != NULL) {
        if (argv[1] != NULL) {
            run(string(argv[1]));
        } else {
            cout << "Configuration file is missing!" << endl;
        }
    }
}
