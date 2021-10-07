#include <iostream>
#include <string>
#include <memory>

#include "VirtualFluids.h"

#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "modify.h"
#include "fix_lb_coupling_onetoone.h"

using namespace std;

//void run(string configname)
//{
//    try {
//        ConfigurationFile config;
//        config.load(configname);
//
//        string pathname            = config.getValue<string>("pathname");
//        int numOfThreads           = config.getValue<int>("numOfThreads");
//        vector<int> blocknx        = config.getVector<int>("blocknx");
//        vector<double> boundingBox = config.getVector<double>("boundingBox");
//        double uLB             = config.getValue<double>("uLB");
//        double nuL             = config.getValue<double>("nuL");
//        double nuG             = config.getValue<double>("nuG");
//        double densityRatio    = config.getValue<double>("densityRatio");
//        double sigma           = config.getValue<double>("sigma");
//        int interfaceThickness = config.getValue<int>("interfaceThickness");
//        double radius          = config.getValue<double>("radius");
//        double theta           = config.getValue<double>("contactAngle");
//        double gr              = config.getValue<double>("gravity");
//        double phiL            = config.getValue<double>("phi_L");
//        double phiH            = config.getValue<double>("phi_H");
//        double tauH            = config.getValue<double>("Phase-field Relaxation");
//        double mob             = config.getValue<double>("Mobility");
//
//        double endTime     = config.getValue<double>("endTime");
//        double outTime     = config.getValue<double>("outTime");
//        double availMem    = config.getValue<double>("availMem");
//        int refineLevel    = config.getValue<int>("refineLevel");
//        double Re          = config.getValue<double>("Re");
//        double dx          = config.getValue<double>("dx");
//        bool logToFile     = config.getValue<bool>("logToFile");
//        //double restartStep = config.getValue<double>("restartStep");
//        //double cpStart     = config.getValue<double>("cpStart");
//        //double cpStep      = config.getValue<double>("cpStep");
//        bool newStart      = config.getValue<bool>("newStart");
//
//        double beta  = 12 * sigma / interfaceThickness;
//        double kappa = 1.5 * interfaceThickness * sigma;
//
//        SPtr<Communicator> comm = MPICommunicator::getInstance();
//        int myid                = comm->getProcessID();
//
//        if (myid == 0)
//            UBLOG(logINFO, "Droplet Test: Start!");
//
//        if (logToFile) {
//#if defined(__unix__)
//            if (myid == 0) {
//                const char *str = pathname.c_str();
//                mkdir(str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
//            }
//#endif
//
//            if (myid == 0) {
//                stringstream logFilename;
//                logFilename << pathname + "/logfile" + UbSystem::toString(UbSystem::getTimeStamp()) + ".txt";
//                UbLog::output_policy::setStream(logFilename.str());
//            }
//        }
//
//        //Sleep(30000);
//
//        // LBMReal dLB = 0; // = length[1] / dx;
//        LBMReal rhoLB = 0.0;
//        LBMReal nuLB  = nuL; //(uLB*dLB) / Re;
//
//        LBMReal D  = 2.0*radius;
//
//        SPtr<LBMUnitConverter> conv(new LBMUnitConverter());
//
//        //const int baseLevel = 0;
//
//        SPtr<LBMKernel> kernel;
//
//        //kernel = SPtr<LBMKernel>(new MultiphaseScratchCumulantLBMKernel());
//        //kernel = SPtr<LBMKernel>(new MultiphaseCumulantLBMKernel());
//        kernel = SPtr<LBMKernel>(new MultiphaseTwoPhaseFieldsPressureFilterLBMKernel());
//
//
//        kernel->setWithForcing(true);
//        kernel->setForcingX1(0.0);
//        kernel->setForcingX2(gr);
//        kernel->setForcingX3(0.0);
//
//        kernel->setPhiL(phiL);
//        kernel->setPhiH(phiH);
//        kernel->setPhaseFieldRelaxation(tauH);
//        kernel->setMobility(mob);
//
//        kernel->setCollisionFactorMultiphase(nuL, nuG);
//        kernel->setDensityRatio(densityRatio);
//        kernel->setMultiphaseModelParameters(beta, kappa);
//        kernel->setContactAngle(theta);
//
//        SPtr<BCProcessor> bcProc(new BCProcessor());
//        // BCProcessorPtr bcProc(new ThinWallBCProcessor());
//
//        kernel->setBCProcessor(bcProc);
//
//        SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
//        noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new MultiphaseNoSlipBCAlgorithm()));
//        //////////////////////////////////////////////////////////////////////////////////
//        // BC visitor
//        MultiphaseBoundaryConditionsBlockVisitor bcVisitor;
//        bcVisitor.addBC(noSlipBCAdapter);
//
//        SPtr<Grid3D> grid(new Grid3D(comm));
//        grid->setDeltaX(dx);
//        grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
//        grid->setPeriodicX1(true);
//        grid->setPeriodicX2(false);
//        grid->setPeriodicX3(true);
//        grid->setGhostLayerWidth(2);
//
//        //////////////////////////////////////////////////////////////////////////
//        // restart
//        //SPtr<UbScheduler> rSch(new UbScheduler(cpStep, cpStart));
//        ////SPtr<MPIIORestartCoProcessor> rcp(new MPIIORestartCoProcessor(grid, rSch, pathname, comm));
//        ////SPtr<MPIIOMigrationCoProcessor> rcp(new MPIIOMigrationCoProcessor(grid, rSch, pathname, comm));
//        //SPtr<MPIIOMigrationBECoProcessor> rcp(new MPIIOMigrationBECoProcessor(grid, rSch, pathname, comm));
//        //rcp->setNu(nuLB);
//        //rcp->setNuLG(nuL, nuG);
//        //rcp->setDensityRatio(densityRatio);
//
//        //rcp->setLBMKernel(kernel);
//        //rcp->setBCProcessor(bcProc);
//        //////////////////////////////////////////////////////////////////////////
//
//        if (newStart) {
//
//            // bounding box
//            double g_minX1 = boundingBox[0];
//            double g_minX2 = boundingBox[2];
//            double g_minX3 = boundingBox[4];
//
//            double g_maxX1 = boundingBox[1];
//            double g_maxX2 = boundingBox[3];
//            double g_maxX3 = boundingBox[5];
//
//            // geometry
//            SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
//            if (myid == 0)
//                GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube",
//                    WbWriterVtkXmlBinary::getInstance());
//
//            if (myid == 0) {
//                UBLOG(logINFO, "uLb = " << uLB);
//                UBLOG(logINFO, "rho = " << rhoLB);
//                UBLOG(logINFO, "nuLb = " << nuLB);
//                UBLOG(logINFO, "Re = " << Re);
//                UBLOG(logINFO, "dx = " << dx);
//                UBLOG(logINFO, "Preprocess - start");
//            }
//
//            GenBlocksGridVisitor genBlocks(gridCube);
//            grid->accept(genBlocks);
//
//            double dx2 = 2.0 * dx;
//            GbCuboid3DPtr wallYmin(new GbCuboid3D(g_minX1 - dx2, g_minX2 - dx2, g_minX3 - dx2, g_maxX1 + dx2, g_minX2, g_maxX3 + dx2));
//            GbSystem3D::writeGeoObject(wallYmin.get(), pathname + "/geo/wallYmin", WbWriterVtkXmlASCII::getInstance());
//            GbCuboid3DPtr wallYmax(new GbCuboid3D(g_minX1 - dx2, g_maxX2, g_minX3 - dx2, g_maxX1 + dx2, g_maxX2 + dx2, g_maxX3 + dx2));
//            GbSystem3D::writeGeoObject(wallYmax.get(), pathname + "/geo/wallYmax", WbWriterVtkXmlASCII::getInstance());
//
//            SPtr<D3Q27Interactor> wallYminInt(new D3Q27Interactor(wallYmin, grid, noSlipBCAdapter, Interactor3D::SOLID));
//            SPtr<D3Q27Interactor> wallYmaxInt(new D3Q27Interactor(wallYmax, grid, noSlipBCAdapter, Interactor3D::SOLID));
// 
//            SPtr<WriteBlocksCoProcessor> ppblocks(new WriteBlocksCoProcessor(
//                grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));
//
//            SPtr<Grid3DVisitor> metisVisitor(
//                new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::RECURSIVE));
//            InteractorsHelper intHelper(grid, metisVisitor, true);
//            intHelper.addInteractor(wallYminInt);
//            intHelper.addInteractor(wallYmaxInt);
//            intHelper.selectBlocks();
//
//            ppblocks->process(0);
//            ppblocks.reset();
//
//            unsigned long long numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
//            int ghostLayer                    = 5;
//            unsigned long long numberOfNodesPerBlock =
//                (unsigned long long)(blocknx[0]) * (unsigned long long)(blocknx[1]) * (unsigned long long)(blocknx[2]);
//            unsigned long long numberOfNodes = numberOfBlocks * numberOfNodesPerBlock;
//            unsigned long long numberOfNodesPerBlockWithGhostLayer =
//                numberOfBlocks * (blocknx[0] + ghostLayer) * (blocknx[1] + ghostLayer) * (blocknx[2] + ghostLayer);
//            double needMemAll =
//                double(numberOfNodesPerBlockWithGhostLayer * (27 * sizeof(double) + sizeof(int) + sizeof(float) * 4));
//            double needMem = needMemAll / double(comm->getNumberOfProcesses());
//
//            if (myid == 0) {
//                UBLOG(logINFO, "Number of blocks = " << numberOfBlocks);
//                UBLOG(logINFO, "Number of nodes  = " << numberOfNodes);
//                int minInitLevel = grid->getCoarsestInitializedLevel();
//                int maxInitLevel = grid->getFinestInitializedLevel();
//                for (int level = minInitLevel; level <= maxInitLevel; level++) {
//                    int nobl = grid->getNumberOfBlocks(level);
//                    UBLOG(logINFO, "Number of blocks for level " << level << " = " << nobl);
//                    UBLOG(logINFO, "Number of nodes for level " << level << " = " << nobl * numberOfNodesPerBlock);
//                }
//                UBLOG(logINFO, "Necessary memory  = " << needMemAll << " bytes");
//                UBLOG(logINFO, "Necessary memory per process = " << needMem << " bytes");
//                UBLOG(logINFO, "Available memory per process = " << availMem << " bytes");
//            }
//
//            MultiphaseSetKernelBlockVisitor kernelVisitor(kernel, nuL, nuG, densityRatio, beta, kappa, theta, availMem,
//                needMem);
//
//            grid->accept(kernelVisitor);
//
//            if (refineLevel > 0) {
//                SetUndefinedNodesBlockVisitor undefNodesVisitor;
//                grid->accept(undefNodesVisitor);
//            }
//
//
//            intHelper.setBC();
//
//            // initialization of distributions
//            LBMReal x1c = 2.5 * D; // (g_maxX1 - g_minX1-1)/2; //
//            LBMReal x2c = 12.5 * D; //(g_maxX2 - g_minX2-1)/2;
//            LBMReal x3c = 2.5 * D; //(g_maxX3 - g_minX3-1)/2;
//            mu::Parser fct1;
//            fct1.SetExpr("0.5-0.5*tanh(2*(sqrt((x1-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)");
//            fct1.DefineConst("x1c", x1c);
//            fct1.DefineConst("x2c", x2c);
//            fct1.DefineConst("x3c", x3c);
//            fct1.DefineConst("radius", radius);
//            fct1.DefineConst("interfaceThickness", interfaceThickness);
//
//            mu::Parser fct2;
//            fct2.SetExpr("0.5*uLB-uLB*0.5*tanh(2*(sqrt((x1-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)");
//            fct2.DefineConst("uLB", uLB);
//            fct2.DefineConst("x1c", x1c);
//            fct2.DefineConst("x2c", x2c);
//            fct2.DefineConst("x3c", x3c);
//            fct2.DefineConst("radius", radius);
//            fct2.DefineConst("interfaceThickness", interfaceThickness);
//
//            MultiphaseInitDistributionsBlockVisitor initVisitor(densityRatio);
//            initVisitor.setPhi(fct1);
//            initVisitor.setVx1(fct2);
//            grid->accept(initVisitor);
//
//            // boundary conditions grid
//            {
//                SPtr<UbScheduler> geoSch(new UbScheduler(1));
//                SPtr<WriteBoundaryConditionsCoProcessor> ppgeo(new WriteBoundaryConditionsCoProcessor(
//                    grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), comm));
//                ppgeo->process(0);
//                ppgeo.reset();
//            }
//
//            if (myid == 0)
//                UBLOG(logINFO, "Preprocess - end");
//        } else {
//            if (myid == 0) {
//                UBLOG(logINFO, "Parameters:");
//                UBLOG(logINFO, "uLb = " << uLB);
//                UBLOG(logINFO, "rho = " << rhoLB);
//                UBLOG(logINFO, "nuLb = " << nuLB);
//                UBLOG(logINFO, "Re = " << Re);
//                UBLOG(logINFO, "dx = " << dx);
//                UBLOG(logINFO, "number of levels = " << refineLevel + 1);
//                UBLOG(logINFO, "numOfThreads = " << numOfThreads);
//                UBLOG(logINFO, "path = " << pathname);
//            }
//
//            //rcp->restart((int)restartStep);
//            //grid->setTimeStep(restartStep);
//
//            if (myid == 0)
//                UBLOG(logINFO, "Restart - end");
//        }
//
//        grid->accept(bcVisitor);
//
//        //TwoDistributionsSetConnectorsBlockVisitor setConnsVisitor(comm);
//        //grid->accept(setConnsVisitor);
//
//        ThreeDistributionsDoubleGhostLayerSetConnectorsBlockVisitor setConnsVisitor(comm);
//        grid->accept(setConnsVisitor);
//
//        SPtr<UbScheduler> visSch(new UbScheduler(outTime));
//        SPtr<WriteMultiphaseQuantitiesCoProcessor> pp(new WriteMultiphaseQuantitiesCoProcessor(
//            grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
//        pp->process(0);
//        //SPtr<WriteMacroscopicQuantitiesCoProcessor> pp(new WriteMacroscopicQuantitiesCoProcessor(
//        //    grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
//
//        SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
//        SPtr<NUPSCounterCoProcessor> npr(new NUPSCounterCoProcessor(grid, nupsSch, numOfThreads, comm));
//
//        omp_set_num_threads(numOfThreads);
//
//        SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
//        SPtr<Calculator> calculator(new BasicCalculator(grid, stepGhostLayer, endTime));
//        calculator->addCoProcessor(npr);
//        calculator->addCoProcessor(pp);
//        //calculator->addCoProcessor(rcp);
//
//
//
//        if (myid == 0)
//            UBLOG(logINFO, "Simulation-start");
//        calculator->calculate();
//        if (myid == 0)
//            UBLOG(logINFO, "Simulation-end");
//    } catch (std::exception &e) {
//        cerr << e.what() << endl << flush;
//    } catch (std::string &s) {
//        cerr << s << endl;
//    } catch (...) {
//        cerr << "unknown exception" << endl;
//    }
//}
int main(int argc, char *argv[])
{
	//void* lmp;
 //   lammps_open_no_mpi(argc, argv, &lmp);
 //   // Sleep(30000);
 //   if (argv != NULL) {
 //       if (argv[1] != NULL) {
 //           run(string(argv[1]));
 //       } else {
 //           cout << "Configuration file is missing!" << endl;
 //       }
 //   }
	
	LAMMPS_NS::LAMMPS *lmp;
    // custom argument vector for LAMMPS library
    const char *lmpargv[] {"liblammps", "-log", "none"};
    int lmpargc = sizeof(lmpargv)/sizeof(const char *);

    // explicitly initialize MPI
    MPI_Init(&argc, &argv);

    // create LAMMPS instance
    lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv, MPI_COMM_WORLD);
    lmp->input->file("in.chute_wear");
    //lmp->input->one("run 1");
    
    //# Try extracting a global value
    //    print("")
    //    print("Attempting to get the number of atoms in simulation")
    //    numAtoms = lmp.extract_global("natoms", 0)
    //    print("natoms =", numAtoms)

    //    # Try extracting atom's positions
    //    print("")
    //    print("Attempting to get the atom's positions")
    //    pos = lmp.extract_atom("x",3)
    //    for k in range(0,numAtoms):
    //        print("Pos[%i] = [%f, %f, %f]" % (k, pos[k][0], pos[k][1], pos[k][2]))

    LAMMPS_NS::FixLbCouplingOnetoone 
        *couplingFix 
        = dynamic_cast<LAMMPS_NS::FixLbCouplingOnetoone*>
        (lmp->modify->find_fix_style("couple/lb/onetoone",0));

    cout << "test1\n";
    
    //double **t_liggghts = couplingFix->get_torque_ptr();
    cout << "test2\n";

    lmp->input->one("run 9 upto");

    for (int step = 0; step < 3; step++)
    {
        

        int numAtoms = lmp->atom->natoms;

        //double** pos = (double**)lmp->atom->extract("x");
        double** pos = lmp->atom->x;
        
        //double* forceX = lmp->atom->fx;

        for (int i = 0; i < numAtoms; i++)
        {
            double **f_liggghts = couplingFix->get_force_ptr();
            double** force = lmp->atom->f;
            cout << "Pos[" << i << "] = [" << pos[i][0] << ", " << pos[i][1] << ", " << pos[i][2] << "]\n";
            cout << "Force[" << i << "] = [" << f_liggghts[i][0] << ", " << f_liggghts[i][1] << ", " << f_liggghts[i][2] << "]\n";
            f_liggghts[i][2] += 500;
            cout << "Force[" << i << "] = [" << force[i][0] << ", " << force[i][1] << ", " << force[i][2] << "]\n";
        }

        couplingFix->comm_force_torque();

        lmp->input->one("run 10");
       
    }

    // delete LAMMPS instance
    delete lmp;

    // stop MPI environment
    MPI_Finalize();
    return 0;
}
