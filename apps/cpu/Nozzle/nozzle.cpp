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

    //double g_minX1 = -1341.81e-3 + 10e-3;
    //double g_minX2 =  0.360872;
    //double g_minX3 = -210e-3;

    //double g_maxX1 = -1260.81e-3 - 10e-3;
    //double g_maxX2 =  0.416302;
    //double g_maxX3 = 210e-3;

    //int blockNX[3] = { 10, 10, 10 };

    double g_minX1 = -1.31431;
    double g_minX2 = 0.375582;
    double g_minX3 = -210e-3 - 1e-3;

    double g_maxX1 = -1.28831;
    double g_maxX2 = 0.401582;
    double g_maxX3 = 0.206;

    int blockNX[3] = { 26, 26, 26 };

    double dx = 1e-3;

    double uLB  = 0.00001;
    //double rhoLB = 0.0;

    // concrete 
    double d_part = 1e-3;
    double V = 0.4;     // flow rate [m^3/h]
    double D = 0.026;   // shotcrete inlet diameter [m]
    double R = D / 2.0; // radius [m]
    double A = UbMath::PI * R * R;
    double u = V / 3600 / A;
    double muConcrete = 2.1133054011798826; // [Pa s]
    double rhoAir = 1.2041;                // [kg/m^3]
    double tau0 = 715.218181094648; //
    double rhoConcrete = 2400; // [kg/m^3]
    double nu = muConcrete / rhoConcrete;

    //double Re_D = d_part * u / nu;
    //if (myid == 0) UBLOG(logINFO, "Re_D = " << Re_D);
    //
    SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(d_part, 1., 2400, d_part / dx, uLB);
    if (myid == 0) std::cout << units->toString() << std::endl;

    double interfaceThickness = 4.096;
    double sigma = 0.03;
    double Re = rhoConcrete * u * d_part / muConcrete;
    double We = rhoConcrete * u * u * d_part / sigma;

    

    double nu_h_LB = uLB * d_part * units->getFactorLentghWToLb() / Re;
    double nu_l_LB = 0;// = nu_h_LB;
    

    double rho_h_LB = 1;

    // surface tension
    double sigma_LB = rho_h_LB * uLB * uLB * d_part * units->getFactorLentghWToLb() / We;
    

    // LBMReal dLB = 0; // = length[1] / dx;
    LBMReal rhoLB = 0.0;
    //LBMReal nuLB = nu_l; //(uLB*dLB) / Re;

    double beta = 12.0 * sigma_LB / interfaceThickness;
    double kappa = 1.5 * interfaceThickness * sigma_LB;

    double phiL = 0.0;
    double phiH = 1.0;
    double tauH = 0.6; // Phase - field Relaxation
    double mob = 0.02; // Mobility
    //double nuL = 1e-2;
    //double nuG = 0.015811388300841892;
    double densityRatio =  rhoConcrete / rhoAir;
    //double sigma_old = 1.0850694444444444e-06;
    //
    //double beta_old = 12.0 * sigma / interfaceThickness;
    //double kappa_old = 1.5 * interfaceThickness * sigma;
    
    double theta = 110; //contact angle

    //https://civilsir.com/density-of-cement-sand-and-aggregate-in-kg-m3-list-of-material-density/

    // SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(r_p, 1.480, 2060, r_p/dx);
    // SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(r_p, LBMUnitConverter::AIR_20C, r_p / dx);
    //SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(d_part, 1., 1000, d_part / dx, std::abs(uLB));
    //SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(d_part, 1., 1000, d_part / dx, std::abs(uLB));
    //SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(d_part, 1., 2400, d_part / dx, uRef);
    

    //SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
    //noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));
    SPtr<BCAdapter> noSlipBCAdapter(new NoSlipBCAdapter());
    noSlipBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new MultiphaseNoSlipBCAlgorithm()));



  
    // concrete inflow boundary condition
    mu::Parser fct;
    fct.SetExpr("U");
    fct.DefineConst("U", -u*units->getFactorVelocityWToLb());
    if (myid == 0) VF_LOG_INFO("Concrete inflow velocity = {} m/s", u);
    if (myid == 0) VF_LOG_INFO("Concrete inflow velocity = {} dx/dt", u * units->getFactorVelocityWToLb());
    if (myid == 0) VF_LOG_INFO("Concrete Re = {}", Re);
        
    //    // Štigler, J. (2014). Analytical velocity profile in tube for laminar and turbulent flow. Engineering
    //    // Mechanics, 21(6), 371-379.
    //    double cx1 = -1.31431 + R;
    //    double cx2 = 0.375582 + R;
    //    //double cx3 = 0.20105 + R;
    //    double L = g_maxX1 - g_minX1;
    //    double p_concrete = 1e5; // Pa = 1 Bar
    //    double p1 = p_concrete * units->getFactorPressureWToLb();
    //    double p2 = 0.0;
    //    double drhoLB = 1.0 + rhoLB;
    //    double muLB = drhoLB * nuLB;
    //    double N = R * R / 2 * muLB * uLB * (p1 - p2) / L - 3;

    //    // mu::Parser fct;
    //    fct.SetExpr("U*(1-(((((x2-y0)^2+(x1-x0)^2)^0.5)/R)^NplusOne))");
    //    fct.DefineConst("x0", cx1);
    //    fct.DefineConst("y0", cx2);
    //    //fct.DefineConst("z0", cx3);
    //    fct.DefineConst("R", R);
    //    fct.DefineConst("U", uLB * ((N + 3) / (N + 1)));
    //    fct.DefineConst("NplusOne", N + 1.0);
    

    //SPtr<BCAdapter> inflowConcreteBCAdapter(new VelocityBCAdapter(false, false, true, fct, 0, BCFunction::INFCONST));
    //inflowConcreteBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityBCAlgorithm()));
    SPtr<BCAdapter> inflowConcreteBCAdapter(new MultiphaseVelocityBCAdapter(false, false, true, fct, phiH, 0, BCFunction::INFCONST));
    inflowConcreteBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new MultiphaseVelocityBCAlgorithm()));

    
        // air inflow boundary condition
        //  Štigler, J. (2014). Analytical velocity profile in tube for laminar and turbulent flow. Engineering
        //  Mechanics, 21(6), 371-379.
        // SPtr<LBMUnitConverter> unitsAir = std::make_shared<LBMUnitConverter>(d_part, LBMUnitConverter::AIR_20C, d_part / dx);
        //SPtr<LBMUnitConverter> unitsAir = std::make_shared<LBMUnitConverter>(d_part, 1., 1.2041, d_part / dx, uLB);
        //double V = 40;      // flow rate [m^3/h]
        //double D = 0.0166;  // air inlet diameter [m]
        //double R = D / 2.0; // radius [m]
        //double A = UbMath::PI * R * R;
        //double u = V / 3600 / A;
        //double uLB = u * unitsAir->getFactorVelocityWToLb();
        //// double cx1 = -1.2788 + R;
        //double cx2 = 0.3803 + R;
        //double cx3 = 0.1517 + R;
        //double L = g_maxX1 - g_minX1;
        //double p_air = 7e5; // Pa = 7 Bar
        //double p1 = p_air;
        //double p2 = 0.0;
        //double mu = 17.2e-6; // Pa s, air 20° C
        //double N = R * R / 2 * mu * u * (p1 - p2) / L - 3;
        //if (myid == 0) VF_LOG_INFO("Air inflow velocity = {} m/s", u);
        //if (myid == 0) VF_LOG_INFO("Air inflow velocity = {} dx/dt", uLB);
        //

        //double nu = mu / rhoConcrete;
        //double Re = d_part * u / nu;
        //if (myid == 0) VF_LOG_INFO("Re_air = {}", Re);

        //double nuLB = d_part * unitsAir->getFactorLentghWToLb() * uLB / Re;
        //if (myid == 0) VF_LOG_INFO("nuLB_air = {}", nuLB);
        //nu_l_LB = nuLB;
    

    SPtr<LBMUnitConverter> unitsAir = std::make_shared<LBMUnitConverter>(d_part, 1., 1.2041, d_part / dx, uLB);
    double V_air = 40;      // flow rate [m^3/h]
    double D_air = 0.00553; // air inlet diameter [m]
    double R_air = D_air / 2.0; // radius [m]
    double A_air = UbMath::PI * R_air * R_air;
    double u_air = V_air / 3600 / A_air;
    double uLB_air = u_air * unitsAir->getFactorVelocityWToLb();
    // double cx1 = -1.2788 + R;
    double cx2 = 0.385822 + R_air;
    double cx3 = 0.135562 + R_air;
    double L_air = 0.00747;
    double p_air = 7e5; // Pa = 7 Bar
    double p1 = p_air;
    double p2 = 1e5;
    double mu_air = 17.2e-6; // Pa s, air 20° C
    double rho_air = 1.2041;  // [kg/m^3]
    double N = R_air * R_air / 2 * mu_air * u_air * (p1 - p2) / L_air - 3;
    if (myid == 0) VF_LOG_INFO("Air inflow velocity = {} m/s", u_air);
    if (myid == 0) VF_LOG_INFO("Air inflow velocity = {} dx/dt", uLB_air);

    double nu_air = mu_air / rho_air;
    double Re_air = d_part * u_air / nu_air;
    if (myid == 0) VF_LOG_INFO("Air Re = {}", Re_air);

    double nuLB_air = d_part * unitsAir->getFactorLentghWToLb() * uLB_air / Re_air;
    if (myid == 0) VF_LOG_INFO("nuLB_air = {}", nuLB_air);
    nu_l_LB = nuLB_air;

    if (myid == 0) VF_LOG_INFO("nu_h = {}", nu_h_LB);
    if (myid == 0) VF_LOG_INFO("nu_l = {}", nu_l_LB);
    if (myid == 0) VF_LOG_INFO("sigma_LB = {}", sigma_LB);

    

    mu::Parser fctVx1;
    //fctVx1.SetExpr("U");
    //fctVx1.DefineConst("U", uLB_air);
    mu::Parser fctVx2;
    fctVx2.SetExpr("U");
    fctVx2.DefineConst("U", 0);
    mu::Parser fctVx3;
    //fctVx3.SetExpr("U");
    //fctVx3.DefineConst("U", -uLB_air);
    
    fctVx1.SetExpr("U*(1-(((((x2-y0)^2+(x3-z0)^2)^0.5)/R)^NplusOne))");
    //fct.DefineConst("x0", cx1);
    fctVx1.DefineConst("y0", cx2);
    fctVx1.DefineConst("z0", cx3);
    fctVx1.DefineConst("R", R);
    fctVx1.DefineConst("U", uLB_air * ((N + 3) / (N + 1)));
    fctVx1.DefineConst("NplusOne", N + 1.0);

    fctVx3.SetExpr("U*(1-(((((x2-y0)^2+(x3-z0)^2)^0.5)/R)^NplusOne))");
    // fc3.DefineConst("x0", cx1);
    fctVx3.DefineConst("y0", cx2);
    fctVx3.DefineConst("z0", cx3);
    fctVx3.DefineConst("R", R);
    fctVx3.DefineConst("U", -uLB_air * ((N + 3) / (N + 1)));
    fctVx3.DefineConst("NplusOne", N + 1.0);
    

    //SPtr<BCAdapter> inflowAirBCAdapter(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
    //inflowAirBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityBCAlgorithm()));
    SPtr<BCAdapter> inflowAirBCAdapter(new MultiphaseVelocityBCAdapter(true, false, true, fctVx1, fctVx3, fctVx3, phiL, 0, BCFunction::INFCONST));
    inflowAirBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new MultiphaseVelocityBCAlgorithm()));

    SPtr<BCAdapter> outflowBCAdapter(new DensityBCAdapter(rhoLB));
    //outflowBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NonEqDensityBCAlgorithm()));
    //SPtr<BCAdapter> outflowBCAdapter(new DensityBCAdapter(rhoLB));
    outflowBCAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new MultiphaseNonReflectingOutflowBCAlgorithm()));
    //////////////////////////////////////////////////////////////////////////////////
    // BC visitor
    //BoundaryConditionsBlockVisitor bcVisitor;♣
    MultiphaseBoundaryConditionsBlockVisitor bcVisitor;
    bcVisitor.addBC(noSlipBCAdapter);
    bcVisitor.addBC(inflowConcreteBCAdapter);
    bcVisitor.addBC(inflowAirBCAdapter);
    bcVisitor.addBC(outflowBCAdapter);

    // SPtr<LBMKernel> kernel   = make_shared<IBcumulantK17LBMKernel>();
    // SPtr<LBMKernel> kernel   = make_shared<CumulantK17LBMKernel>();
    // SPtr<LBMKernel> kernel = make_shared<MultiphaseTwoPhaseFieldsPressureFilterLBMKernel>();
    SPtr<LBMKernel> kernel = make_shared<MultiphaseSimpleVelocityBaseExternalPressureLBMKernel>();

    kernel->setWithForcing(true);
    kernel->setForcingX1(0.0);
    kernel->setForcingX2(0.0);
    kernel->setForcingX3(0.0);

    kernel->setPhiL(phiL);
    kernel->setPhiH(phiH);
    kernel->setPhaseFieldRelaxation(tauH);
    kernel->setMobility(mob);
    kernel->setInterfaceWidth(interfaceThickness);

    kernel->setCollisionFactorMultiphase(nu_h_LB, nu_l_LB);
    kernel->setDensityRatio(densityRatio);
    kernel->setMultiphaseModelParameters(beta, kappa);
    kernel->setContactAngle(theta);

    SPtr<BCProcessor> bcProc = make_shared<BCProcessor>();
    kernel->setBCProcessor(bcProc);

    SPtr<Grid3D> grid = make_shared<Grid3D>(comm);
    grid->setPeriodicX1(false);
    grid->setPeriodicX2(false);
    grid->setPeriodicX3(false);
    grid->setDeltaX(dx);
    grid->setBlockNX(blockNX[0], blockNX[1], blockNX[2]);
    grid->setGhostLayerWidth(2);

    string geoPath = "d:/Projects/TRR277/Project/WP4/NozzleGeo";

    string outputPath = "d:/temp/NozzleFlowTest_Multiphase2";
    UbSystem::makeDirectory(outputPath);
    UbSystem::makeDirectory(outputPath + "/liggghts");

    //if (myid == 0) {
    //    stringstream logFilename;
    //    logFilename << outputPath + "/logfile" + UbSystem::toString(UbSystem::getTimeStamp()) + ".txt";
    //    UbLog::output_policy::setStream(logFilename.str());
    //}

    SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, vf::lbm::dir::DIR_MMM, MetisPartitioner::RECURSIVE));
    
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
    SPtr<Interactor3D> intrNozzleVolcanNozzle2 = std::make_shared<D3Q27TriFaceMeshInteractor>(meshNozzleVolcanNozzle2, grid, noSlipBCAdapter, Interactor3D::SOLID, Interactor3D::POINTS);
    ///////////////////////////////////////////////////////////
    //box
    SPtr<D3Q27Interactor> intrBox = SPtr<D3Q27Interactor>(new D3Q27Interactor(gridCube, grid, noSlipBCAdapter, Interactor3D::INVERSESOLID));
    ///////////////////////////////////////////////////////////
    //inflow
    GbCylinder3DPtr geoInflow(new GbCylinder3D(-1.30181+0.0005, 0.390872-0.00229, 0.20105, -1.30181+0.0005, 0.390872-0.00229, 0.23, 0.013));
    if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), outputPath + "/geo/geoInflow", WbWriterVtkXmlBinary::getInstance());
    SPtr<D3Q27Interactor> intrInflow = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, inflowConcreteBCAdapter, Interactor3D::SOLID));
    ///////////////////////////////////////////////////////////
    //outflow
    GbCylinder3DPtr geoOutflow(new GbCylinder3D(-1.30181+0.0005, 0.390872-0.00229, -0.22, -1.30181+0.0005, 0.390872-0.00229, -0.21, 0.013));
    //GbCylinder3DPtr geoOutflow(new GbCylinder3D(-1.30181+0.0005, 0.390872-0.00229, g_minX3, -1.30181+0.0005, 0.390872-0.00229, -0.21, 0.013));
    if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), outputPath + "/geo/geoOutflow", WbWriterVtkXmlBinary::getInstance());
    SPtr<D3Q27Interactor> intrOutflow = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, outflowBCAdapter, Interactor3D::SOLID));
    ///////////////////////////////////////////////////////////
    //SPtr<GbTriFaceMesh3D> geoAirInlet = std::make_shared<GbTriFaceMesh3D>();
    //if (myid == 0) UBLOG(logINFO, "Read Air_Inlet:start");
    //geoAirInlet->readMeshFromSTLFileASCII(geoPath + "/Air_Inlet.stl", true);
    //if (myid == 0) UBLOG(logINFO, "Read Air_Inlet:end");
    //if (myid == 0) GbSystem3D::writeGeoObject(geoAirInlet.get(), outputPath + "/geo/geoAirInlet", WbWriterVtkXmlBinary::getInstance());
    //SPtr<Interactor3D> intrAirInlet = std::make_shared<D3Q27TriFaceMeshInteractor>(geoAirInlet, grid, inflowAirBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES);
    /////////////////////////////////////////////////////////////
    //Fluid area
    GbCylinder3DPtr geoFluidArea(new GbCylinder3D(-1.30181+0.0005, 0.390872-0.00229, g_minX3, -1.30181+0.0005, 0.390872-0.00229, g_maxX3, 0.013));
    if (myid == 0) GbSystem3D::writeGeoObject(geoFluidArea.get(), outputPath + "/geo/geoFluidArea", WbWriterVtkXmlBinary::getInstance());
    SPtr<D3Q27Interactor> intrFluidArea = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoFluidArea, grid, noSlipBCAdapter, Interactor3D::INVERSESOLID));
    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    GbCylinder3DPtr geoAirInflow(new GbCylinder3D(-1.31431 - 0.0005, 0.388587, 0.1383275, -1.31431, 0.388587, 0.1383275, 0.002765));
    if (myid == 0) GbSystem3D::writeGeoObject(geoAirInflow.get(), outputPath + "/geo/geoAirInlet", WbWriterVtkXmlBinary::getInstance());
    SPtr<Interactor3D> intrAirInflow = std::make_shared<D3Q27Interactor>(geoAirInflow, grid, inflowAirBCAdapter, Interactor3D::SOLID, Interactor3D::EDGES);
    ///////////////////////////////////////////////////////////

    InteractorsHelper intHelper(grid, metisVisitor, true);
    
    intHelper.addInteractor(intrFluidArea);
    intHelper.addInteractor(intrNozzleVolcanNozzle2);
    //intHelper.addInteractor(intrBox);
    intHelper.addInteractor(intrInflow);
    intHelper.addInteractor(intrAirInflow);
    intHelper.addInteractor(intrOutflow);
    

    //intHelper.addInteractor(intrNozzleAirDistributor);
    //intHelper.addInteractor(intrNozzleAirInlet);
    //intHelper.addInteractor(intrNozzleSpacer);
    //intHelper.addInteractor(intrNozzleAccDistributor);
    //intHelper.addInteractor(intrNozzleAccInlet);
    //intHelper.addInteractor(intrNozzleVolcanNozzle1);
    


    intHelper.selectBlocks();

    SPtr<CoProcessor> ppblocks = make_shared<WriteBlocksCoProcessor>(
         grid, SPtr<UbScheduler>(new UbScheduler(1)), outputPath, WbWriterVtkXmlBinary::getInstance(), comm);
     ppblocks->process(0);
     ppblocks.reset();

     if (myid == 0) UBLOG(logINFO, Utilities::toString(grid, comm->getNumberOfProcesses()));


    //SetKernelBlockVisitor kernelVisitor(kernel, nuLB, comm->getNumberOfProcesses());
     MultiphaseSetKernelBlockVisitor kernelVisitor(kernel, nu_h_LB, nu_l_LB, 1e9, 1);
    grid->accept(kernelVisitor);

    intHelper.setBC();

    //InitDistributionsBlockVisitor initVisitor;
    //grid->accept(initVisitor);

    double x1c = -1.31431 + R;
    double x2c = 0.375582 + R;
    double x3c = 0.20105;

    mu::Parser fct1;
    //fct1.SetExpr(" 0.5 - 0.5 * tanh(2 * (sqrt((x1 - x1c) ^ 2 + (x2 - x2c) ^ 2 + (x3 - x3c) ^ 2) - radius) / interfaceThickness)");
    fct1.SetExpr(" 0.5 - 0.5 * tanh(2 * (sqrt((x1 - x1c) ^ 2 + (x2 - x2c) ^ 2 + (x3 - x3c) ^ 2) - radius) / interfaceThickness)");
    fct1.DefineConst("x1c", x1c);
    fct1.DefineConst("x2c", x2c);
    fct1.DefineConst("x3c", x3c);
    fct1.DefineConst("radius", R);
    fct1.DefineConst("interfaceThickness", interfaceThickness * dx);

    MultiphaseVelocityFormInitDistributionsBlockVisitor initVisitor;
    initVisitor.setPhi(fct1);
    grid->accept(initVisitor);

  
    string inFile1 = "d:/Projects/VirtualFluids_Develop/apps/cpu/Nozzle/in.nozzle";
    //string inFile2 = "d:/Projects/VirtualFluids_LIGGGHTS_coupling/apps/cpu/LiggghtsApp/in2.lbdem";
    MPI_Comm mpi_comm = *(MPI_Comm*)(comm->getNativeCommunicator());
    LiggghtsCouplingWrapper wrapper(argv, mpi_comm);

    double v_frac = 0.1;
    double dt_phys   = units->getFactorTimeLbToW();
    int demSubsteps = 10;
    double dt_dem   = dt_phys / (double)demSubsteps;
    int vtkSteps    = 1000;
    string demOutDir = outputPath + "/liggghts";

    //wrapper.execCommand("echo none");

    //wrapper.execFile((char*)inFile1.c_str());

    //// set timestep and output directory
    wrapper.setVariable("t_step", dt_dem);
    wrapper.setVariable("dmp_stp", vtkSteps * demSubsteps);
    wrapper.setVariable("dmp_dir", demOutDir);

    //wrapper.execFile((char *)inFile1.c_str());
    //wrapper.runUpto(demSubsteps - 1);
    //wrapper.runUpto(1000);

    SPtr<UbScheduler> lScheduler = make_shared<UbScheduler>(1); 
    SPtr<LiggghtsCouplingCoProcessor> lcCoProcessor =
        make_shared<LiggghtsCouplingCoProcessor>(grid, lScheduler, comm, wrapper, demSubsteps, units);

    // boundary conditions grid
    {
        SPtr<UbScheduler> geoSch(new UbScheduler(1));
        SPtr<WriteBoundaryConditionsCoProcessor> ppgeo(new WriteBoundaryConditionsCoProcessor(grid, geoSch, outputPath, WbWriterVtkXmlBinary::getInstance(), comm));
        ppgeo->process(0);
        ppgeo.reset();
    }

    grid->accept(bcVisitor);

    //OneDistributionSetConnectorsBlockVisitor setConnsVisitor(comm);
    //TwoDistributionsDoubleGhostLayerSetConnectorsBlockVisitor setConnsVisitor(comm);
    ThreeDistributionsDoubleGhostLayerSetConnectorsBlockVisitor setConnsVisitor(comm);
    grid->accept(setConnsVisitor);

    int numOfThreads          = 18;
    omp_set_num_threads(numOfThreads);

    SPtr<UbScheduler> nupsSch = std::make_shared<UbScheduler>(10, 10, 100);
    SPtr<NUPSCounterCoProcessor> nupsCoProcessor = make_shared<NUPSCounterCoProcessor>(grid, nupsSch, numOfThreads, comm);

    //// write data for visualization of macroscopic quantities
    SPtr < UbScheduler> visSch(new UbScheduler(vtkSteps));
    //SPtr<UbScheduler> visSch(new UbScheduler(1, 8700, 8800));
   // visSch->addSchedule(1, 8700, 8800);
    SPtr<WriteMultiphaseQuantitiesCoProcessor> writeMQCoProcessor(
        new WriteMultiphaseQuantitiesCoProcessor(grid, visSch, outputPath, WbWriterVtkXmlASCII::getInstance(),
                                                  SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));
    writeMQCoProcessor->process(0);

    int endTime = 1000000;
    SPtr<Calculator> calculator(new BasicCalculator(grid, lScheduler, endTime));
    calculator->addCoProcessor(nupsCoProcessor);
   // calculator->addCoProcessor(lcCoProcessor);
    calculator->addCoProcessor(writeMQCoProcessor);

    if (myid == 0) UBLOG(logINFO, "Simulation-start");
    calculator->calculate();
    if (myid == 0) UBLOG(logINFO, "Simulation-end");


    return 0;
}
