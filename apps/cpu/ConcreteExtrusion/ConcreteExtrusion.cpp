#include <iostream>
#include <memory>
#include <string>

#include "MultiphaseFlow/MultiphaseFlow.h"
#include "NonNewtonianFluids/NonNewtonianFluids.h"
#include "VirtualFluids.h"

using namespace std;

void setInflowBC(real x1, real x2, real x3, real radius, int dir)
{
}

void run(string configname)
{
    using namespace vf::lbm::dir;

    try {

        // Sleep(30000);

        vf::basics::ConfigurationFile config;
        config.load(configname);

        string pathname = config.getValue<string>("pathname");
        string pathGeo = config.getValue<string>("pathGeo");
        int numOfThreads = config.getValue<int>("numOfThreads");
        vector<int> blocknx = config.getVector<int>("blocknx");
        real U_LB = config.getValue<real>("U_LB");
        int interfaceWidth = config.getValue<int>("interfaceWidth");
        real theta = config.getValue<real>("contactAngle");
        real D_LB = config.getValue<real>("D_LB");
        real phiL = config.getValue<real>("phi_L");
        real phiH = config.getValue<real>("phi_H");
        real tauH = config.getValue<real>("Phase-field Relaxation");
        real mob = config.getValue<real>("Mobility");

        real endTime = config.getValue<real>("endTime");
        real outTime = config.getValue<real>("outTime");
        real availMem = config.getValue<real>("availMem");

        bool logToFile = config.getValue<bool>("logToFile");
        real restartStep = config.getValue<real>("restartStep");
        real cpStart = config.getValue<real>("cpStart");
        real cpStep = config.getValue<real>("cpStep");
        bool newStart = config.getValue<bool>("newStart");

        int caseN = config.getValue<int>("case");

        real factorLx = config.getValue<real>("factorLx");
        real factorLy = config.getValue<real>("factorLy");
        real factorLz = config.getValue<real>("factorLz");


        SPtr<vf::mpi::Communicator> comm = vf::mpi::MPICommunicator::getInstance();
        int myid = comm->getProcessID();

        if (myid == 0) UBLOG(logINFO, "Jet Breakup: Start!");

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

        real rho_h = 0, rho_l = 0, r_rho = 0, mu_h = 0, /*mu_l,*/ Uo = 0, D = 0, sigma = 0;
        real Dg = 10;

        switch (caseN) {
            case 1:
                // density of heavy fluid (kg/m^3)
                rho_h = 848;
                // density of light fluid (kg/m^3)
                rho_l = 34.5;
                // density ratio
                r_rho = rho_h / rho_l;
                // dynamic viscosity of heavy fluid (Pa � s)
                mu_h = 2.87e-3;
                // dynamic viscosity of light fluid (Pa � s)
                // mu_l = 1.97e-5;
                // velocity (m/s)
                Uo = 100;
                // diameter of jet (m)
                D = 0.0001;
                // surface tension (N/m)
                sigma = 0.03;
                break;
            case 2:
                // density of heavy fluid (kg/m^3)
                rho_h = 848;
                // density of light fluid (kg/m^3)
                rho_l = 1.205;
                // density ratio
                r_rho = rho_h / rho_l;
                // dynamic viscosity of heavy fluid (Pa � s)
                mu_h = 2.87e-3;
                // dynamic viscosity of light fluid (Pa � s)
                // mu_l = 1.84e-5;
                // velocity (m/s)
                Uo = 200;
                // diameter of jet (m)
                D = 0.0001;
                // surface tension (N/m)
                sigma = 0.03;
                break;
            case 3:
                // density of heavy fluid (kg/m^3)
                rho_h = 1000;
                // density of light fluid (kg/m^3)
                rho_l = 1.0;
                // density ratio
                r_rho = rho_h / rho_l;
                // dynamic viscosity of heavy fluid (Pa � s)
                mu_h = 2.87e-3;
                // dynamic viscosity of light fluid (Pa � s)
                // mu_l = 1.84e-5;
                // velocity (m/s)
                Uo = 200;
                // diameter of jet (m)
                D = 0.0001;
                Dg = 10;
                // surface tension (N/m)
                sigma = 0.03;
                break;
        }

        real Re = rho_h * Uo * D / mu_h;
        real We = rho_h * Uo * Uo * D / sigma;

        real dx = Dg / D_LB;
        real nu_h = U_LB * D_LB / Re;
        real nu_l = nu_h;
        nu_h *= 0.1;

        real rho_h_LB = 1;
        // surface tension
        real sigma_LB = 0.0; //rho_h_LB *U_LB *U_LB *D_LB / We;

        // LBMReal dLB = 0; // = length[1] / dx;
        real rhoLB = 0.0;
        // LBMReal nuLB = nu_l; //(uLB*dLB) / Re;

        real beta = 12.0 * sigma_LB / interfaceWidth;
        real kappa = 1.5 * interfaceWidth * sigma_LB;

        double tau0 = 715.218181094648*1000.; // Pa
        double muConcrete = 2.1133054011798826; // [Pa s]
        real u = Uo; //[m/s]


        double Bm = (tau0 * D) / (muConcrete * u);
        double tau0_LB = 0.02;
        //Bm *nu_h *U_LB / (D / dx);

        SPtr<Rheology> rheo = Rheology::getInstance();
        rheo->setYieldStress(tau0_LB);


        if (myid == 0) {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "U_LB = " << U_LB);
            UBLOG(logINFO, "rho = " << rhoLB);
            UBLOG(logINFO, "nu_l = " << nu_l);
            UBLOG(logINFO, "nu_h = " << nu_h);
            UBLOG(logINFO, "Re = " << Re);
            UBLOG(logINFO, "We = " << We);
            UBLOG(logINFO, "dx = " << dx);
            UBLOG(logINFO, "sigma = " << sigma);
            UBLOG(logINFO, "sigma_LB = " << sigma_LB);
            UBLOG(logINFO, "tau0_LB = " << tau0_LB);
            UBLOG(logINFO, "density ratio = " << r_rho);
            // UBLOG(logINFO, "number of levels = " << refineLevel + 1);
            UBLOG(logINFO, "numOfThreads = " << numOfThreads);
            UBLOG(logINFO, "path = " << pathname);
        }

        // bounding box

        real Lx = factorLx * Dg;
        real Ly = factorLy * Dg;
        real Lz = factorLz * Dg;

        real g_minX1 = 0;
        real g_minX2 = -0.5 * Ly;
        real g_minX3 = -0.5 * Lz;

        // double g_maxX1 = 8.0*D;
        // double g_maxX2 = 2.5*D;
        // double g_maxX3 = 2.5*D;

        real g_maxX1 = Lx;
        real g_maxX2 = 0.5 * Ly;
        real g_maxX3 = 0.5 * Lz;

        SPtr<LBMUnitConverter> conv(new LBMUnitConverter());

        // const int baseLevel = 0;

        SPtr<LBMKernel> kernel;

        // kernel = SPtr<LBMKernel>(new MultiphaseScratchCumulantLBMKernel());
        // kernel = SPtr<LBMKernel>(new MultiphaseCumulantLBMKernel());
        // kernel = SPtr<LBMKernel>(new MultiphaseTwoPhaseFieldsCumulantLBMKernel());
        // kernel = SPtr<LBMKernel>(new MultiphaseTwoPhaseFieldsVelocityCumulantLBMKernel());
        // kernel = SPtr<LBMKernel>(new MultiphaseTwoPhaseFieldsPressureFilterLBMKernel());
        // kernel = SPtr<LBMKernel>(new MultiphasePressureFilterLBMKernel());
        // kernel = SPtr<LBMKernel>(new MultiphaseSimpleVelocityBaseExternalPressureLBMKernel());
        kernel = make_shared<MultiphaseScaleDistributionLBMKernel>();
        //kernel = make_shared<MultiphaseSharpInterfaceLBMKernel>();

        kernel->setWithForcing(false);
        kernel->setForcingX1(0.0);
        kernel->setForcingX2(0.0);
        kernel->setForcingX3(0.0);

        kernel->setPhiL(phiL);
        kernel->setPhiH(phiH);
        kernel->setPhaseFieldRelaxation(tauH);
        kernel->setMobility(mob);

        // nuL, nuG, densityRatio, beta, kappa, theta,

        kernel->setCollisionFactorMultiphase(nu_h, nu_l);
        kernel->setDensityRatio(r_rho);
        kernel->setMultiphaseModelParameters(beta, kappa);
        kernel->setContactAngle(theta);
        kernel->setInterfaceWidth(interfaceWidth);
        // dynamicPointerCast<MultiphasePressureFilterLBMKernel>(kernel)->setPhaseFieldBC(0.0);
        kernel->setSigma(sigma_LB);

        SPtr<BCSet> bcProc(new BCSet());
        // BCSetPtr bcProc(new ThinWallBCSet());

        kernel->setBCSet(bcProc);

        SPtr<Grid3D> grid(new Grid3D(comm));
         //grid->setPeriodicX1(true);
         grid->setPeriodicX2(true);
         grid->setPeriodicX3(true);
        grid->setGhostLayerWidth(2);

        SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, DIR_MMM, MetisPartitioner::RECURSIVE));

        //////////////////////////////////////////////////////////////////////////
        // restart
        SPtr<UbScheduler> rSch(new UbScheduler(cpStep, cpStart));
        // SPtr<MPIIORestartSimulationObserver> rcp(new MPIIORestartSimulationObserver(grid, rSch, pathname, comm));
        SPtr<MPIIOMigrationSimulationObserver> rcp(new MPIIOMigrationSimulationObserver(grid, rSch, metisVisitor, pathname, comm));
        // SPtr<MPIIOMigrationBESimulationObserver> rcp(new MPIIOMigrationBESimulationObserver(grid, rSch, pathname, comm));
        // rcp->setNu(nuLB);
        // rcp->setNuLG(nuL, nuG);
        // rcp->setDensityRatio(densityRatio);

        rcp->setLBMKernel(kernel);
        rcp->setBCSet(bcProc);
        //////////////////////////////////////////////////////////////////////////
        // BC Adapter
        //////////////////////////////////////////////////////////////////////////////
        mu::Parser fctF1;
        //u_actual(r) = U_max * (4 / (D^2) * (R^2 - r^2))
        //fctF1.SetExpr("vx1*(1-((x2-y0)^2+(x3-z0)^2)/(R^2))");
        //fctF1.SetExpr("vx1*(1-(sqrt((x2-y0)^2+(x3-z0)^2))/(R))");
        //fctF1.SetExpr("vy1*(1-(sqrt((x2-x0)^2+(x3-z0)^2)/R))^0.1");
        fctF1.SetExpr("vx1");
        //fctF1.DefineConst("vx1", 0);
        fctF1.DefineConst("vx1", U_LB);
        fctF1.DefineConst("R", 0.5*Dg);
        fctF1.DefineConst("y0", (g_minX2+g_maxX2)/2);
        fctF1.DefineConst("z0", (g_minX3+g_maxX3)/2);
        SPtr<BC> velBCF1(new MultiphaseVelocityBC(true, false, false, fctF1, phiH, 0.0, BCFunction::INFCONST));

        mu::Parser fctF2;
        fctF2.SetExpr("vx1");
        fctF2.DefineConst("vx1", U_LB);
        fctF2.SetExpr("vx1");
        fctF2.DefineConst("vx1", U_LB);

        // real startTime = 1;
        // SPtr<BC> velBCF1(new MultiphaseVelocityBC(true, false, false, fctF1, phiH, 0.0, startTime));
        //SPtr<BC> velBCF2(new MultiphaseVelocityBC(true, false, false, fctF2, phiH, startTime, endTime));
        SPtr<BC> velBCF2(new MultiphaseVelocityBC(false, true, false, fctF2, phiH, 0.0, BCFunction::INFCONST));



        SPtr<BC> noSlipBC(new NoSlipBC());
        noSlipBC->setBCStrategy(SPtr<BCStrategy>(new MultiphaseNoSlipBCStrategy()));

        SPtr<BC> denBC(new DensityBC(rhoLB));
        denBC->setBCStrategy(SPtr<BCStrategy>(new MultiphasePressureBCStrategy()));

        SPtr<BC> slipBC(new SlipBC());
        slipBC->setBCStrategy(SPtr<BCStrategy>(new MultiphaseSlipBCStrategy()));

        mu::Parser fctPhi_F1;
        fctPhi_F1.SetExpr("phiH");
        fctPhi_F1.DefineConst("phiH", phiH);

        mu::Parser fctPhi_F2;
        fctPhi_F2.SetExpr("phiL");
        fctPhi_F2.DefineConst("phiL", phiL);

        mu::Parser fctvel_F2_init;
        fctvel_F2_init.SetExpr("U");
        fctvel_F2_init.DefineConst("U", 0);

        velBCF1->setBCStrategy(SPtr<BCStrategy>(new MultiphaseVelocityBCStrategy()));
        velBCF2->setBCStrategy(SPtr<BCStrategy>(new MultiphaseVelocityBCStrategy()));


        //////////////////////////////////////////////////////////////////////////////////
        // BC visitor
        MultiphaseBoundaryConditionsBlockVisitor bcVisitor;
        bcVisitor.addBC(noSlipBC);
        //bcVisitor.addBC(slipBC);
        bcVisitor.addBC(denBC);
        bcVisitor.addBC(velBCF1);
        bcVisitor.addBC(velBCF2);

        // SPtr<D3Q27Interactor> inflowF1Int;
        // SPtr<D3Q27Interactor> cylInt;

        SPtr<D3Q27Interactor> inflowInt;

        //if (newStart) {

            //  if (newStart) {

            real inflowLength = 10;

            // geometry
            SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
            if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

            // if (myid == 0)
            //     UBLOG(logINFO, "Read geoFile:start");
            // SPtr<GbTriFaceMesh3D> cylinder = make_shared<GbTriFaceMesh3D>();
            // cylinder->readMeshFromSTLFileBinary(pathGeo + "/" + geoFile, false);
            // GbSystem3D::writeGeoObject(cylinder.get(), pathname + "/geo/Stlgeo", WbWriterVtkXmlBinary::getInstance());
            // if (myid == 0)
            //     UBLOG(logINFO, "Read geoFile:stop");
            //  inflow
            //  GbCuboid3DPtr geoInflowF1(new GbCuboid3D(g_minX1, g_minX2 - 0.5 * dx, g_minX3, g_maxX1, g_minX2 - 1.0 *
            //  dx, g_maxX3));
            // GbCuboid3DPtr geoInflowF1(new GbCuboid3D(g_minX1 * 0.5 - dx, g_minX2 - dx, g_minX3 * 0.5 - dx,
            //                                          g_maxX1 * 0.5 + dx, g_minX2, g_maxX3 * 0.5 + dx));
            // if (myid == 0)
            //     GbSystem3D::writeGeoObject(geoInflowF1.get(), pathname + "/geo/geoInflowF1",
            //                                WbWriterVtkXmlASCII::getInstance());

            //GbCylinder3DPtr geoInflow(new GbCylinder3D(g_minX1 - 2.0*dx, 0.0, 0.0, g_minX1, 0.0, 0.0, Dg / 2.0));

            GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1 - 2.0 * dx, g_minX2 + inflowLength, g_minX3 + (g_maxX3 - g_minX3) / 3.0, g_minX1 + 2.0 * dx, 3 * inflowLength + g_minX2, g_maxX3 - (g_maxX3 - g_minX3) / 3.0));

            if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), pathname + "/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

            GbCylinder3DPtr geoSolid(new GbCylinder3D(g_minX1 - 2.0 * dx, g_maxX2 / 2.0, g_maxX3 / 2.0, g_minX1+2.0*dx, g_maxX2 / 2.0, g_maxX3 / 2.0, 1.5*D / 2.0));
            if (myid == 0) GbSystem3D::writeGeoObject(geoSolid.get(), pathname + "/geo/geoSolid", WbWriterVtkXmlASCII::getInstance());

            SPtr<GbTriFaceMesh3D> meshInflowPipe = std::make_shared<GbTriFaceMesh3D>();
            if (myid == 0) UBLOG(logINFO, "Read meshInflowPipe:start");
            meshInflowPipe->readMeshFromSTLFileBinary(pathGeo + "/JetTube4.stl", false);
            //meshInflowPipe->readMeshFromSTLFileASCII(pathGeo + "/JetTubeScaled5.stl", false);
            if (myid == 0) UBLOG(logINFO, "Read meshInflowPipe:end");
            //meshInflowPipe->scale(1e-04, 1e-04, 1e-04);
            if (myid == 0) GbSystem3D::writeGeoObject(meshInflowPipe.get(), pathname + "/geo/meshInflowPipe", WbWriterVtkXmlBinary::getInstance());
            SPtr<Interactor3D> intrInflowPipe = std::make_shared<D3Q27TriFaceMeshInteractor>(meshInflowPipe, grid, noSlipBC, Interactor3D::SOLID, Interactor3D::POINTS);

            // GbCylinder3DPtr cylinder2(
            //    new GbCylinder3D(0.0, g_minX2 - 2.0 * dx / 2.0, 0.0, 0.0, g_minX2 + 4.0 * dx, 0.0, 8.0+2.0*dx));
            // if (myid == 0)
            //    GbSystem3D::writeGeoObject(cylinder2.get(), pathname + "/geo/cylinder2",
            //                               WbWriterVtkXmlASCII::getInstance());
            // outflow
            // GbCuboid3DPtr geoOutflow(new GbCuboid3D(-1.0, -1, -1.0, 121.0, 1.0, 121.0)); // For JetBreakup (Original)
            // GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_minX1, g_maxX2 - 40 * dx, g_minX3, g_maxX1, g_maxX2, g_maxX3));
            GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2 - 2.0 * dx, g_minX3 - 2.0 * dx, g_maxX1 + 2.0 * dx, g_maxX2 + 2.0 * dx, g_maxX3 + 2.0 * dx));
            if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), pathname + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

            // double blockLength = blocknx[0] * dx;

            if (myid == 0) {
                UBLOG(logINFO, "Preprocess - start");
            }

            grid->setDeltaX(dx);
            grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

            //grid->setPeriodicX1(false);
            //grid->setPeriodicX2(false);
            //grid->setPeriodicX3(false);

            GenBlocksGridVisitor genBlocks(gridCube);
            grid->accept(genBlocks);

            SPtr<WriteBlocksSimulationObserver> ppblocks(new WriteBlocksSimulationObserver(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

            // SPtr<Interactor3D> tubes(new D3Q27TriFaceMeshInteractor(cylinder, grid, noSlipBC,
            //                                                         Interactor3D::SOLID, Interactor3D::POINTS));

            // inflowF1Int =
            //    SPtr<D3Q27Interactor>(new D3Q27Interactor(cylinder1, grid, noSlipBC, Interactor3D::SOLID));
            // inflowF1Int->addBC(velBCF2);

            SPtr<D3Q27Interactor> outflowInt(new D3Q27Interactor(geoOutflow, grid, denBC, Interactor3D::SOLID));

            // Create boundary conditions geometry
            GbCuboid3DPtr wallXmin(new GbCuboid3D(g_minX1 - 2.0 * dx, g_minX2 - 2.0 * dx, g_minX3 - 2.0 * dx, g_minX1, g_maxX2 + 2.0 * dx, g_maxX3 + 2.0 * dx));
            GbSystem3D::writeGeoObject(wallXmin.get(), pathname + "/geo/wallXmin", WbWriterVtkXmlASCII::getInstance());
            GbCuboid3DPtr wallXmax(new GbCuboid3D(g_maxX1, g_minX2 - 2.0 * dx, g_minX3 - 2.0 * dx, g_maxX1 + 2.0 * dx, g_maxX2 + 2.0 * dx, g_maxX3 + 2.0 * dx));
            GbSystem3D::writeGeoObject(wallXmax.get(), pathname + "/geo/wallXmax", WbWriterVtkXmlASCII::getInstance());
            GbCuboid3DPtr wallZmin(new GbCuboid3D(g_minX1 - 2.0*dx, g_minX2 - 2.0*dx, g_minX3 - 2.0*dx, g_maxX1 + 2.0*dx, g_maxX2 + 2.0*dx, g_minX3));
            GbSystem3D::writeGeoObject(wallZmin.get(), pathname + "/geo/wallZmin", WbWriterVtkXmlASCII::getInstance());
            GbCuboid3DPtr wallZmax(new GbCuboid3D(g_minX1 - 2.0*dx, g_minX2 - 2.0*dx, g_maxX3, g_maxX1 + 2.0*dx, g_maxX2 + 2.0*dx, g_maxX3 + 2.0*dx));
            GbSystem3D::writeGeoObject(wallZmax.get(), pathname + "/geo/wallZmax", WbWriterVtkXmlASCII::getInstance());
            GbCuboid3DPtr wallYmin(new GbCuboid3D(g_minX1 - 2.0 * dx, g_minX2 - 2.0 * dx, g_minX3 - 2.0 * dx, g_maxX1 + 2.0 * dx, g_minX2, g_maxX3 + 2.0 * dx));
            GbSystem3D::writeGeoObject(wallYmin.get(), pathname + "/geo/wallYmin", WbWriterVtkXmlASCII::getInstance());
            GbCuboid3DPtr wallYmax(new GbCuboid3D(g_minX1 - 2.0 * dx, g_maxX2, g_minX3 - 2.0 * dx, g_maxX1 + 2.0 * dx, g_maxX2 + 2.0 * dx, g_maxX3 + 2.0 * dx));
            GbSystem3D::writeGeoObject(wallYmax.get(), pathname + "/geo/wallYmax", WbWriterVtkXmlASCII::getInstance());

            // Add boundary conditions to grid generator
            SPtr<D3Q27Interactor> wallXminInt(new D3Q27Interactor(wallXmin, grid, noSlipBC, Interactor3D::SOLID));
            SPtr<D3Q27Interactor> wallXmaxInt(new D3Q27Interactor(wallXmax, grid, velBCF2, Interactor3D::SOLID));
            SPtr<D3Q27Interactor> wallZminInt(new D3Q27Interactor(wallZmin, grid, noSlipBC, Interactor3D::SOLID));
            SPtr<D3Q27Interactor> wallZmaxInt(new D3Q27Interactor(wallZmax, grid, noSlipBC, Interactor3D::SOLID));
            SPtr<D3Q27Interactor> wallYminInt(new D3Q27Interactor(wallYmin, grid, noSlipBC, Interactor3D::SOLID));
            SPtr<D3Q27Interactor> wallYmaxInt(new D3Q27Interactor(wallYmax, grid, noSlipBC, Interactor3D::SOLID));


//////////////////////////
            
            GbCuboid3DPtr wallInfZmin(new GbCuboid3D(g_minX1 - 2.0 * dx, g_minX2 - 2.0 * dx, g_minX3 - 2.0 * dx, g_minX1 + inflowLength * dx, g_maxX2 + 2.0 * dx, (g_maxX3 + g_minX3) / 2.0 - inflowLength));
            GbSystem3D::writeGeoObject(wallInfZmin.get(), pathname + "/geo/wallInfZmin", WbWriterVtkXmlASCII::getInstance());
            GbCuboid3DPtr wallInfZmax(new GbCuboid3D(g_minX1 - 2.0 * dx, g_minX2 - 2.0 * dx, (g_maxX3 + g_minX3) / 2.0 + inflowLength, g_minX1 + inflowLength * dx, g_maxX2 + 2.0 * dx, g_maxX3 + 2.0 * dx));
            GbSystem3D::writeGeoObject(wallInfZmax.get(), pathname + "/geo/wallInfZmax", WbWriterVtkXmlASCII::getInstance());
            GbCuboid3DPtr wallInfYmin(new GbCuboid3D(g_minX1 - 2.0 * dx, g_minX2 - 2.0 * dx, g_minX3 - 2.0 * dx, g_minX1 + inflowLength * dx, g_minX2 + inflowLength, g_maxX3 + 2.0 * dx));
            GbSystem3D::writeGeoObject(wallInfYmin.get(), pathname + "/geo/wallInfYmin", WbWriterVtkXmlASCII::getInstance());
            GbCuboid3DPtr wallInfYmax(new GbCuboid3D(g_minX1 - 2.0 * dx, 3* inflowLength+g_minX2, g_minX3 - 2.0 * dx, g_minX1 + inflowLength * dx, g_maxX2 + 2.0 * dx, g_maxX3));
            GbSystem3D::writeGeoObject(wallInfYmax.get(), pathname + "/geo/wallInfYmax", WbWriterVtkXmlASCII::getInstance());

            // Add boundary conditions to grid generator
            SPtr<D3Q27Interactor> wallInfZminInt(new D3Q27Interactor(wallInfZmin, grid, noSlipBC, Interactor3D::SOLID));
            SPtr<D3Q27Interactor> wallInfZmaxInt(new D3Q27Interactor(wallInfZmax, grid, noSlipBC, Interactor3D::SOLID));
            SPtr<D3Q27Interactor> wallInfYminInt(new D3Q27Interactor(wallInfYmin, grid, noSlipBC, Interactor3D::SOLID));
            SPtr<D3Q27Interactor> wallInfYmaxInt(new D3Q27Interactor(wallInfYmax, grid, noSlipBC, Interactor3D::SOLID));
//////////////////////////////////////





            // cylInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(cylinder1, grid, velBCF1, Interactor3D::SOLID));
            // cylInt->addBC(velBCF2);
            //  SPtr<D3Q27Interactor> cyl2Int(new D3Q27Interactor(cylinder2, grid, noSlipBC,
            //  Interactor3D::SOLID));

            inflowInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, velBCF1, Interactor3D::SOLID));
            //inflowInt->addBC(velBCF2);

            GbCylinder3DPtr geoAirInflow(new GbCylinder3D(g_minX1 - 2.0*dx, g_maxX2 / 2.0, g_maxX3 / 2.0, g_minX1, g_maxX2 / 2.0, g_maxX3 / 2.0, D * 0.9 / 2.0));
            if (myid == 0) GbSystem3D::writeGeoObject(geoAirInflow.get(), pathname + "/geo/geoAirInflow", WbWriterVtkXmlASCII::getInstance());            
            SPtr<D3Q27Interactor> inflowAirInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoAirInflow, grid, velBCF2, Interactor3D::SOLID));

            SPtr<D3Q27Interactor> solidInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoSolid, grid, noSlipBC, Interactor3D::SOLID));


            InteractorsHelper intHelper(grid, metisVisitor, true);
            //intHelper.addInteractor(cylInt);
            //intHelper.addInteractor(tubes);
            //intHelper.addInteractor(outflowInt);
            // intHelper.addInteractor(cyl2Int);
            //intHelper.addInteractor(intrInflowPipe);
            intHelper.addInteractor(wallXminInt);
            intHelper.addInteractor(wallXmaxInt);
            //intHelper.addInteractor(wallZminInt);
            //intHelper.addInteractor(wallZmaxInt);
            //intHelper.addInteractor(wallYminInt);
            //intHelper.addInteractor(wallYmaxInt);

            intHelper.addInteractor(wallInfZminInt);
            intHelper.addInteractor(wallInfZmaxInt);
            intHelper.addInteractor(wallInfYminInt);
            intHelper.addInteractor(wallInfYmaxInt);

            intHelper.addInteractor(inflowInt);
            //intHelper.addInteractor(outflowInt);
            //intHelper.addInteractor(inflowAirInt);

            intHelper.selectBlocks();

            ppblocks->update(0);
            ppblocks.reset();

            unsigned long long numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
            int ghostLayer = 3;
            unsigned long long numberOfNodesPerBlock = (unsigned long long)(blocknx[0]) * (unsigned long long)(blocknx[1]) * (unsigned long long)(blocknx[2]);
            unsigned long long numberOfNodes = numberOfBlocks * numberOfNodesPerBlock;
            unsigned long long numberOfNodesPerBlockWithGhostLayer = numberOfBlocks * (blocknx[0] + ghostLayer) * (blocknx[1] + ghostLayer) * (blocknx[2] + ghostLayer);
            real needMemAll = real(numberOfNodesPerBlockWithGhostLayer * (27 * sizeof(real) + sizeof(int) + sizeof(float) * 4));
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

            MultiphaseSetKernelBlockVisitor kernelVisitor(kernel, nu_h, nu_l, availMem, needMem);

            grid->accept(kernelVisitor);

            //if (!newStart) {
            //    rcp->readBlocks((int)restartStep);
            //    grid->accept(metisVisitor);
            //    rcp->readDataSet((int)restartStep);
            //    grid->setTimeStep(restartStep);
            //}

            intHelper.setBC();

        if (newStart) {

                // initialization of distributions
                // mu::Parser fct1;
                // fct1.SetExpr("phiL");
                // fct1.DefineConst("phiL", phiL);
                real x1c = g_minX1 - Dg*5; // (g_maxX1 - g_minX1-1)/2; //
                real x2c = (g_minX2 + g_maxX2) / 2;
                real x3c = (g_minX3 + g_maxX3) / 2;

                mu::Parser fct1;
                fct1.SetExpr("0.5-0.5*tanh(2*(sqrt((x1-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)");
                //fct1.SetExpr("(0.5-0.5*tanh(2*(sqrt((x1/radius-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness))-(0.5-0.5*tanh(2*(sqrt((x1/radius-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius2)/interfaceThickness))");
                // fct1.SetExpr("x1 < 4*dx ? (0.5-0.5*tanh(2*(1-4*dx)*(sqrt((x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)) : 0");
                fct1.DefineConst("x1c", x1c);
                fct1.DefineConst("x2c", x2c);
                fct1.DefineConst("x3c", x3c);
                fct1.DefineConst("dx", dx);
                fct1.DefineConst("radius", 5.0*Dg + (inflowLength-1)*dx /* + 2. * dx*/);
                fct1.DefineConst("radius2", Dg / 4.);
                fct1.DefineConst("interfaceThickness", interfaceWidth * dx);

                mu::Parser fct2;
                // fct1.SetExpr("0.5-0.5*tanh(2*(sqrt((x1/radius-x1c)^2+(x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)");
                fct2.SetExpr("x1 < 4*dx ? (0.5-0.5*tanh(2*(1-4*dx)*(sqrt((x2-x2c)^2+(x3-x3c)^2)-radius)/interfaceThickness)) : 0");
                fct2.DefineConst("x1c", x1c);
                fct2.DefineConst("x2c", x2c);
                fct2.DefineConst("x3c", x3c);
                fct2.DefineConst("dx", dx);
                fct2.DefineConst("radius", 0.5 * D /* + 2. * dx*/);
                fct2.DefineConst("interfaceThickness", interfaceWidth * dx);

                MultiphaseVelocityFormInitDistributionsBlockVisitor initVisitor;
                //initVisitor.setPhi(fct1);
                grid->accept(initVisitor);
            }

            // boundary conditions grid
            {
                SPtr<UbScheduler> geoSch(new UbScheduler(1));
                SPtr<WriteBoundaryConditionsSimulationObserver> ppgeo(new WriteBoundaryConditionsSimulationObserver(grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), comm));
                ppgeo->update(0);
                ppgeo.reset();
            }

            if (myid == 0) UBLOG(logINFO, "Preprocess - end");
        //} else {
         if (!newStart) {
            rcp->restart((int)restartStep);
            grid->setTimeStep(restartStep);
            SPtr<WriteBlocksSimulationObserver> ppblocks(new WriteBlocksSimulationObserver(grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));
            ppblocks->update(10000);
            if (myid == 0) UBLOG(logINFO, "Restart - end");
        }

        //  TwoDistributionsSetConnectorsBlockVisitor setConnsVisitor(comm);
        //  grid->accept(setConnsVisitor);

        // ThreeDistributionsSetConnectorsBlockVisitor setConnsVisitor(comm);

        grid->accept(bcVisitor);

        //ThreeDistributionsDoubleGhostLayerSetConnectorsBlockVisitor setConnsVisitor(comm);
        TwoDistributionsDoubleGhostLayerSetConnectorsBlockVisitor setConnsVisitor(comm);
        grid->accept(setConnsVisitor);

        SPtr<UbScheduler> visSch(new UbScheduler(outTime));
        real t_ast, t;
        t_ast = 7.19;
        t = (int)(t_ast/(U_LB/(D_LB)));
        visSch->addSchedule(t,t,t); //t=7.19
        // SPtr<WriteMultiphaseQuantitiesSimulationObserver> pp(new WriteMultiphaseQuantitiesSimulationObserver(
        //     grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));
        SPtr<WriteSharpInterfaceQuantitiesSimulationObserver> pp(new WriteSharpInterfaceQuantitiesSimulationObserver(grid, visSch, pathname, WbWriterVtkXmlBinary::getInstance(), conv, comm));        
        pp->update(0);

        SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
        SPtr<NUPSCounterSimulationObserver> npr(new NUPSCounterSimulationObserver(grid, nupsSch, numOfThreads, comm));

        // SPtr<UbScheduler> timeBCSch(new UbScheduler(1, startTime, startTime));
        // auto timeDepBC = make_shared<TimeDependentBCSimulationObserver>(TimeDependentBCSimulationObserver(grid, timeBCSch));
        // timeDepBC->addInteractor(inflowInt);

#ifdef _OPENMP
        omp_set_num_threads(numOfThreads);
#endif

        SPtr<UbScheduler> stepGhostLayer(new UbScheduler(1));
        SPtr<Simulation> simulation(new Simulation(grid, stepGhostLayer, endTime));
        simulation->addSimulationObserver(npr);
        simulation->addSimulationObserver(pp);
        // simulation->addSimulationObserver(timeDepBC);
        // simulation->addSimulationObserver(timeDepBC);
        simulation->addSimulationObserver(rcp);

        if (myid == 0) UBLOG(logINFO, "Simulation-start");
        simulation->run();
        if (myid == 0) UBLOG(logINFO, "Simulation-end");
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
