#include <iostream>
#include <string>
#include <memory>

#if defined(__unix__)
#include <stdio.h>
#include <stdlib.h>
#endif

#include "VirtualFluids.h"
#include "MultiphaseFlow/MultiphaseFlow.h"
#include "NonNewtonianFluids/NonNewtonianFluids.h"

using namespace std;

void run(string configname)
{
    using namespace vf::lbm::dir;

    try {
        vf::basics::ConfigurationFile config;
        config.load(configname);

        string pathname            = config.getValue<string>("pathname");
        int numOfThreads           = config.getValue<int>("numOfThreads");
        vector<int> blocknx        = config.getVector<int>("blocknx");
        vector<real> boundingBox = config.getVector<real>("boundingBox");
        real uLB             = config.getValue<real>("uLB");
        real nuL             = config.getValue<real>("nuL");
        real nuG             = config.getValue<real>("nuG");
        real densityRatio    = config.getValue<real>("densityRatio");
        real sigma           = config.getValue<real>("sigma");
        int interfaceThickness = config.getValue<int>("interfaceThickness");
        real radius          = config.getValue<real>("radius");
        real theta           = config.getValue<real>("contactAngle");
        //double gr              = config.getValue<double>("gravity");
        real phiL            = config.getValue<real>("phi_L");
        real phiH            = config.getValue<real>("phi_H");
        real tauH            = config.getValue<real>("Phase-field Relaxation");
        real mob             = config.getValue<real>("Mobility");

        real endTime     = config.getValue<real>("endTime");
        real outTime     = config.getValue<real>("outTime");
        real availMem    = config.getValue<real>("availMem");
        int refineLevel    = config.getValue<int>("refineLevel");
        real Re          = config.getValue<real>("Re");
        real dx          = config.getValue<real>("dx");
        bool logToFile     = config.getValue<bool>("logToFile");
        real restartStep = config.getValue<real>("restartStep");
        real cpStart     = config.getValue<real>("cpStart");
        real cpStep      = config.getValue<real>("cpStep");
        bool newStart      = config.getValue<bool>("newStart");
        //double rStep = config.getValue<double>("rStep");

        SPtr<vf::mpi::Communicator> comm = vf::mpi::MPICommunicator::getInstance();
        int myid                = comm->getProcessID();

        if (myid == 0)
            UBLOG(logINFO, "Droplet Test: Start!");

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
        
        std::string fileName = "./LastTimeStep" + std::to_string((int)boundingBox[1]) + ".txt";

//#if defined(__unix__)
//         double lastTimeStep = 0;
//         //if (!newStart) 
//         {
//             std::ifstream ifstr(fileName);
//             ifstr >> lastTimeStep;
//             restartStep = lastTimeStep;
//             if(endTime >= lastTimeStep)
//                endTime = lastTimeStep + rStep;
//             else
//                return;
//         }    
//#endif

        //Sleep(30000);

        // LBMReal dLB = 0; // = length[1] / dx;
        real rhoLB = 0.0;
        real nuLB  = nuL; //(uLB*dLB) / Re;

        //diameter of circular droplet
        real D  = 2.0*radius;

        //density retio
        real r_rho = densityRatio;

        //density of heavy fluid
        real rho_h = 1.0;
        //density of light fluid
        real rho_l = rho_h / r_rho;

        //kinimatic viscosity
        real nu_h = nuL;
        //LBMReal nu_l = nuG;
        //#dynamic viscosity
        real mu_h = rho_h * nu_h;
        
        //gravity
        real g_y = Re* Re* mu_h* mu_h / (rho_h * (rho_h - rho_l) * D * D * D);
        //Eotvos number
        real Eo = 100;
        //surface tension
        sigma = rho_h* g_y* D* D / Eo;

        //g_y = 0;

        real beta  = 12.0 * sigma / interfaceThickness;
        real kappa = 1.5 * interfaceThickness * sigma;

        if (myid == 0) {
                //UBLOG(logINFO, "uLb = " << uLB);
                //UBLOG(logINFO, "rho = " << rhoLB);
                UBLOG(logINFO, "D = " << D);
                UBLOG(logINFO, "nuL = " << nuL);
                UBLOG(logINFO, "nuG = " << nuG);
                UBLOG(logINFO, "Re = " << Re);
                UBLOG(logINFO, "Eo = " << Eo);
                UBLOG(logINFO, "g_y = " << g_y);
                UBLOG(logINFO, "sigma = " << sigma);
                UBLOG(logINFO, "dx = " << dx);
                UBLOG(logINFO, "Preprocess - start");
        }

        SPtr<LBMUnitConverter> conv(new LBMUnitConverter());

        //const int baseLevel = 0;

        
        SPtr<Rheology> thix = Rheology::getInstance();
        thix->setYieldStress(0);

        SPtr<LBMKernel> kernel;

        //kernel = SPtr<LBMKernel>(new MultiphaseScratchCumulantLBMKernel());
       // kernel = SPtr<LBMKernel>(new MultiphaseCumulantLBMKernel());
        //kernel = SPtr<LBMKernel>(new MultiphaseTwoPhaseFieldsPressureFilterLBMKernel());
        //kernel = SPtr<LBMKernel>(new MultiphasePressureFilterLBMKernel());
        //kernel = make_shared< MultiphaseScaleDistributionLBMKernel>();
        //kernel = SPtr<LBMKernel>(new MultiphaseSharpInterfaceLBMKernel());
        kernel = make_shared<MultiphaseSharpInterfaceLBMKernel>();

        mu::Parser fgr;
        fgr.SetExpr("-(rho-rho_l)*g_y");
        fgr.DefineConst("rho_l", rho_l);
        fgr.DefineConst("g_y", g_y);

        kernel->setWithForcing(true);
        kernel->setForcingX1(0.0);
        kernel->setForcingX2(fgr);
        kernel->setForcingX3(0.0);

        kernel->setPhiL(phiL);
        kernel->setPhiH(phiH);
        kernel->setPhaseFieldRelaxation(tauH);
        kernel->setMobility(mob);
        kernel->setInterfaceWidth(interfaceThickness);


        kernel->setCollisionFactorMultiphase(nuL, nuG);
        kernel->setDensityRatio(densityRatio);
        kernel->setMultiphaseModelParameters(beta, kappa);
        kernel->setContactAngle(theta);

        SPtr<BCSet> bcProc(new BCSet());
        // BCSetPtr bcProc(new ThinWallBCSet());

        kernel->setBCSet(bcProc);

        SPtr<BC> noSlipBC(new NoSlipBC());
        noSlipBC->setBCStrategy(SPtr<BCStrategy>(new MultiphaseNoSlipBCStrategy()));
        //////////////////////////////////////////////////////////////////////////////////
        // BC visitor
        MultiphaseBoundaryConditionsBlockVisitor bcVisitor;
        bcVisitor.addBC(noSlipBC);

        SPtr<Grid3D> grid(new Grid3D(comm));
        grid->setDeltaX(dx);
        grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);
        grid->setPeriodicX1(true);
        grid->setPeriodicX2(false);
        grid->setPeriodicX3(true);
        grid->setGhostLayerWidth(2);

        SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, DIR_MMM, MetisPartitioner::RECURSIVE));

        //////////////////////////////////////////////////////////////////////////
        // restart
        SPtr<UbScheduler> rSch(new UbScheduler(cpStep, cpStart));
        //SPtr<MPIIORestartSimulationObserver> rcp(new MPIIORestartSimulationObserver(grid, rSch, pathname, comm));
        SPtr<MPIIOMigrationSimulationObserver> rcp(new MPIIOMigrationSimulationObserver(grid, rSch, metisVisitor, pathname, comm));
        //SPtr<MPIIOMigrationBESimulationObserver> rcp(new MPIIOMigrationBESimulationObserver(grid, rSch, pathname, comm));
        // rcp->setNu(nuLB);
        // rcp->setNuLG(nuL, nuG);
        // rcp->setDensityRatio(densityRatio);

        rcp->setLBMKernel(kernel);
        rcp->setBCSet(bcProc);
        //////////////////////////////////////////////////////////////////////////

        if (newStart) {

            // bounding box
            real g_minX1 = boundingBox[0];
            real g_minX2 = boundingBox[2];
            real g_minX3 = boundingBox[4];

            real g_maxX1 = boundingBox[1];
            real g_maxX2 = boundingBox[3];
            real g_maxX3 = boundingBox[5];

            // geometry
            SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
            if (myid == 0)
                GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube",
                    WbWriterVtkXmlBinary::getInstance());



            GenBlocksGridVisitor genBlocks(gridCube);
            grid->accept(genBlocks);

            real dx2 = 2.0 * dx;
            GbCuboid3DPtr wallYmin(new GbCuboid3D(g_minX1 - dx2, g_minX2 - dx2, g_minX3 - dx2, g_maxX1 + dx2, g_minX2, g_maxX3 + dx2));
            GbSystem3D::writeGeoObject(wallYmin.get(), pathname + "/geo/wallYmin", WbWriterVtkXmlASCII::getInstance());
            GbCuboid3DPtr wallYmax(new GbCuboid3D(g_minX1 - dx2, g_maxX2, g_minX3 - dx2, g_maxX1 + dx2, g_maxX2 + dx2, g_maxX3 + dx2));
            GbSystem3D::writeGeoObject(wallYmax.get(), pathname + "/geo/wallYmax", WbWriterVtkXmlASCII::getInstance());

            SPtr<D3Q27Interactor> wallYminInt(new D3Q27Interactor(wallYmin, grid, noSlipBC, Interactor3D::SOLID));
            SPtr<D3Q27Interactor> wallYmaxInt(new D3Q27Interactor(wallYmax, grid, noSlipBC, Interactor3D::SOLID));
 
            SPtr<WriteBlocksSimulationObserver> ppblocks(new WriteBlocksSimulationObserver(
                grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

            InteractorsHelper intHelper(grid, metisVisitor, true);
            intHelper.addInteractor(wallYminInt);
            intHelper.addInteractor(wallYmaxInt);
            intHelper.selectBlocks();

            ppblocks->update(0);
            ppblocks.reset();

            unsigned long long numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
            int ghostLayer                    = 5;
            unsigned long long numberOfNodesPerBlock =
                (unsigned long long)(blocknx[0]) * (unsigned long long)(blocknx[1]) * (unsigned long long)(blocknx[2]);
            unsigned long long numberOfNodes = numberOfBlocks * numberOfNodesPerBlock;
            unsigned long long numberOfNodesPerBlockWithGhostLayer =
                numberOfBlocks * (blocknx[0] + ghostLayer) * (blocknx[1] + ghostLayer) * (blocknx[2] + ghostLayer);
            real needMemAll =
                real(numberOfNodesPerBlockWithGhostLayer * (27 * sizeof(real) + sizeof(int) + sizeof(float) * 4));
            real needMem = needMemAll / real(comm->getNumberOfProcesses());

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
            real x1c = 1.5;//
            //2.5 * D;             // (g_maxX1 - g_minX1-1)/2; //
            real x2c = 12.5 * D; //(g_maxX2 - g_minX2-1)/2;
            real x3c = 2.5 * D;
            //1.5; // 2.5 * D; //(g_maxX3 - g_minX3-1)/2;
            //LBMReal x3c = 2.5 * D;
            mu::Parser fct1;
            fct1.SetExpr("0.5+0.5*tanh(2*(sqrt(0*(x1-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)");
            fct1.DefineConst("x1c", x1c);
            fct1.DefineConst("x2c", x2c);
            fct1.DefineConst("x3c", x3c);
            fct1.DefineConst("radius", radius);
            fct1.DefineConst("interfaceThickness", interfaceThickness);

            mu::Parser fct2;
            fct2.SetExpr("0.5*uLB-uLB*0.5*tanh(2*(sqrt(0*(x1-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)");
            //fct2.SetExpr("uLB");
            fct2.DefineConst("uLB", uLB);
            fct2.DefineConst("x1c", x1c);
            fct2.DefineConst("x2c", x2c);
            fct2.DefineConst("x3c", x3c);
            fct2.DefineConst("radius", radius);
            fct2.DefineConst("interfaceThickness", interfaceThickness);

            //MultiphaseInitDistributionsBlockVisitor initVisitor(densityRatio);
            MultiphaseVelocityFormInitDistributionsBlockVisitor initVisitor;
            initVisitor.setPhi(fct1);
            initVisitor.setVx1(fct2);
            grid->accept(initVisitor);

            // boundary conditions grid
            {
                SPtr<UbScheduler> geoSch(new UbScheduler(1));
                SPtr<WriteBoundaryConditionsSimulationObserver> ppgeo(new WriteBoundaryConditionsSimulationObserver(
                    grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), comm));
                ppgeo->update(0);
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

        grid->accept(bcVisitor);

        //TwoDistributionsSetConnectorsBlockVisitor setConnsVisitor(comm);
        //grid->accept(setConnsVisitor);

        //ThreeDistributionsDoubleGhostLayerSetConnectorsBlockVisitor setConnsVisitor(comm);
        //grid->accept(setConnsVisitor);

        TwoDistributionsDoubleGhostLayerSetConnectorsBlockVisitor setConnsVisitor(comm);
        grid->accept(setConnsVisitor);

        SPtr<UbScheduler> visSch(new UbScheduler(outTime));
        real t_ast, t;
        t_ast = 2;
        t = (int)(t_ast/std::sqrt(g_y/D));
        visSch->addSchedule(t,t,t); //t=2
        t_ast = 3;
        t = (int)(t_ast/std::sqrt(g_y/D));        
        visSch->addSchedule(t,t,t); //t=3
        t_ast = 4;
        t = (int)(t_ast/std::sqrt(g_y/D));        
        visSch->addSchedule(t,t,t); //t=4
        t_ast = 5;
        t = (int)(t_ast/std::sqrt(g_y/D));        
        visSch->addSchedule(t,t,t); //t=5
        t_ast = 6;
        t = (int)(t_ast/std::sqrt(g_y/D)); 
        visSch->addSchedule(t,t,t); //t=6
        t_ast = 7;
        t = (int)(t_ast/std::sqrt(g_y/D));         
        visSch->addSchedule(t,t,t); //t=7
        t_ast = 9;
        t = (int)(t_ast/std::sqrt(g_y/D));         
        visSch->addSchedule(t,t,t); //t=9

        SPtr<WriteSharpInterfaceQuantitiesSimulationObserver> pp(new WriteSharpInterfaceQuantitiesSimulationObserver(
            grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
        if(grid->getTimeStep() == 0) 
            pp->update(0);

        SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
        SPtr<NUPSCounterSimulationObserver> npr(new NUPSCounterSimulationObserver(grid, nupsSch, numOfThreads, comm));

        omp_set_num_threads(numOfThreads);

        SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
        SPtr<Simulation> simulation(new Simulation(grid, stepGhostLayer, endTime));
        simulation->addSimulationObserver(npr);
        simulation->addSimulationObserver(pp);
        simulation->addSimulationObserver(rcp);


        if (myid == 0)
            UBLOG(logINFO, "Simulation-start");
        simulation->run();
        if (myid == 0)
            UBLOG(logINFO, "Simulation-end");
            
//#if defined(__unix__)
//         //if (!newStart) 
//         //{
//            if (myid == 0) 
//            {
//                std::ofstream ostr(fileName);
//                ostr << endTime;
//                cout << "start sbatch\n";
//                //system("./start.sh");
//                //system("echo test!");
//                std::string str = "sbatch startJob" + std::to_string((int)boundingBox[1]) + ".sh";
//                //system("sbatch startJob512.sh");
//                system(str.c_str());
//            }   
//            //MPI_Barrier((MPI_Comm)comm->getNativeCommunicator()); 
//         //}
//#endif

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
