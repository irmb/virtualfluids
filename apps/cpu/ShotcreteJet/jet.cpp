#include <iostream>
#include <memory>
#include <string>

#include "VirtualFluids.h"

//#include "LiggghtsCoupling/LiggghtsCoupling.h"

#include "MultiphaseFlow/MultiphaseFlow.h"

#include "NonNewtonianFluids/NonNewtonianFluids.h"

using namespace std;

int main(int argc, char *argv[])
{
    //Sleep(30000);
    string configname;
    if (argv != NULL) {
        if (argv[1] != NULL) {
            configname = string(argv[1]);
        } else {
            cout << "Configuration file is missing!" << endl;
            return 0;
        }
    }


    try {

        vf::basics::ConfigurationFile config;
        config.load(configname);

        string outputPath = config.getValue<string>("outputPath");
        string geoPath = config.getValue<string>("geoPath");
        bool logToFile = config.getValue<bool>("logToFile");
        int vtkSteps = config.getValue<int>("vtkSteps");
        bool newStart = config.getValue<bool>("newStart");
        double cpStep  = config.getValue<double>("cpStep"); 
        double cpStart = config.getValue<double>("cpStart");
        double restartStep = config.getValue<int>("restartStep");
        int endTime = config.getValue<int>("endTime");

        std::shared_ptr<vf::mpi::Communicator> comm = vf::mpi::MPICommunicator::getInstance();
        int myid = comm->getProcessID();

        if (myid == 0) UBLOG(logINFO, "Jet Breakup: Start!");

        if (logToFile) {
#if defined(__unix__)
            if (myid == 0) {
                const char *str = outputPath.c_str();
                mkdir(str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            }
#endif

            if (myid == 0) {
                stringstream logFilename;
                logFilename << outputPath + "/logfile" + UbSystem::toString(UbSystem::getTimeStamp()) + ".txt";
                UbLog::output_policy::setStream(logFilename.str());
            }
        }        

        //bool newStart = false;
        // bounding box
        // double g_minX1 = -1341.81e-3;
        // double g_minX2 =  348.087e-3;
        // double g_minX3 = -210e-3;

        // double g_maxX1 = -1260.81e-3;
        // double g_maxX2 =  429.087e-3;
        // double g_maxX3 =  214.5e-3;

        // double g_minX1 = -1341.81e-3 + 10e-3;
        // double g_minX2 =  0.360872;
        // double g_minX3 = -210e-3;

        // double g_maxX1 = -1260.81e-3 - 10e-3;
        // double g_maxX2 =  0.416302;
        // double g_maxX3 = 210e-3;

        // int blockNX[3] = { 10, 10, 10 };

        //int gridNZ = 3;

        //double g_minX1 = -1.31431;
        //double g_minX2 = 0.375582;
        //double g_minX3 = -0.21 + 0.035 * 8.0; //-0.21; //-210e-3 - 0.2 - 6e-3; //- 1e-3;

        //double g_maxX1 = -1.28831;
        //double g_maxX2 = 0.401582;
        //double g_maxX3 = 0.175;//0.21;

        double dx = 1e-3;

        double g_maxX3_box = -0.065;        

        double g_minX1 = -1.49631;
        double g_minX2 = 0.193582;
        double g_minX3 = g_maxX3_box - 0.03;//-0.095; //-0.215; 

        double g_maxX1 = -1.10631;
        double g_maxX2 = 0.583582;
        double g_maxX3 = 0.175;

      

        //int blockNX[3] = { 26, 26, 35 };
        int blockNX[3] = { 15, 15, 15 };


        double uLB_ref = 0.0001;
        // double rhoLB = 0.0;

        // concrete
        double d_part = 1e-3;
        double V = 0.4; // flow rate [m^3/h]
        double D = 0.026;     // shotcrete inlet diameter [m]
        double R = D / 2.0;   // radius [m]
        double A = UbMath::PI * R * R;
        double u = V / 3600 / A;
        double muConcrete = 2.1133054011798826; // [Pa s]
        double rhoAir = 1.2041;                 // [kg/m^3]
        double tau0 = 715.218181094648;         // Pa
        double rhoConcrete = 2400;              // [kg/m^3]
        double nu = muConcrete / rhoConcrete;

        // double Re_D = d_part * u / nu;
        // if (myid == 0) UBLOG(logINFO, "Re_D = " << Re_D);
        //
        SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(d_part, 1., 2400, d_part / dx, uLB_ref);
        if (myid == 0) std::cout << units->toString() << std::endl;

        double interfaceThickness = 3; // 4.096;
        double sigma = 0.3; //0.03;
        double Re = rhoConcrete * u * d_part / muConcrete;
        double We = rhoConcrete * u * u * d_part / sigma;

        double u_LB_con = u * units->getFactorVelocityWToLb();
        double nu_h_LB = nu * units->getFactorViscosityWToLb(); // uLB_ref * d_part * units->getFactorLentghWToLb() / Re;
        double nu_l_LB = 0;                                     // = nu_h_LB;

        double rho_h_LB = 1;

        // surface tension
        double sigma_LB = rho_h_LB *u_LB_con *u_LB_con *d_part * units->getFactorLentghWToLb() / We;

        // LBMReal dLB = 0; // = length[1] / dx;
        LBMReal rhoLB = 0.0;
        // LBMReal nuLB = nu_l; //(uLB_ref*dLB) / Re;

        double beta = 12.0 * sigma_LB / interfaceThickness;
        double kappa = 1.5 * interfaceThickness * sigma_LB;

        double phiL = 0.0;
        double phiH = 1.0;
        double tauH = 0.6; // Phase - field Relaxation
        double mob = 0.02; // Mobility
        // double nuL = 1e-2;
        // double nuG = 0.015811388300841892;
        double densityRatio = rhoConcrete / rhoAir;
        // double sigma_old = 1.0850694444444444e-06;
        //
        // double beta_old = 12.0 * sigma / interfaceThickness;
        // double kappa_old = 1.5 * interfaceThickness * sigma;

        double theta = 110; // contact angle

        // https://civilsir.com/density-of-cement-sand-and-aggregate-in-kg-m3-list-of-material-density/

        // SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(r_p, 1.480, 2060, r_p/dx);
        // SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(r_p, LBMUnitConverter::AIR_20C, r_p / dx);
        // SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(d_part, 1., 1000, d_part / dx, std::abs(uLB_ref));
        // SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(d_part, 1., 1000, d_part / dx, std::abs(uLB_ref));
        // SPtr<LBMUnitConverter> units = std::make_shared<LBMUnitConverter>(d_part, 1., 2400, d_part / dx, uRef);

        double Bm = (tau0 * d_part) / (muConcrete * u);
        double tau0_LB = Bm * nu_h_LB * u_LB_con / (d_part * units->getFactorLentghWToLb());

        SPtr<Rheology> thix = Rheology::getInstance();
        thix->setYieldStress(tau0_LB);

        if (myid == 0) VF_LOG_INFO("Yield stress = {} Pa", tau0);
        if (myid == 0) VF_LOG_INFO("Yield stress LB = {} ", tau0_LB);

        //SPtr<BC> noSlipBC(new NoSlipBC());
        //noSlipBC->setBCStrategy(SPtr<BCStrategy>(new NoSlipBCStrategy()));
        SPtr<BC> noSlipBC(new NoSlipBC());
        noSlipBC->setBCStrategy(SPtr<BCStrategy>(new MultiphaseNoSlipBCStrategy()));

        // concrete inflow boundary condition
        mu::Parser fct;
        fct.SetExpr("U");
        fct.DefineConst("U", -u_LB_con);
        if (myid == 0) VF_LOG_INFO("Concrete inflow velocity = {} m/s", u);
        if (myid == 0) VF_LOG_INFO("Concrete inflow velocity = {} dx/dt", u_LB_con);
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
        //    double N = R * R / 2 * muLB * uLB_ref * (p1 - p2) / L - 3;

        //    // mu::Parser fct;
        //    fct.SetExpr("U*(1-(((((x2-y0)^2+(x1-x0)^2)^0.5)/R)^NplusOne))");
        //    fct.DefineConst("x0", cx1);
        //    fct.DefineConst("y0", cx2);
        //    //fct.DefineConst("z0", cx3);
        //    fct.DefineConst("R", R);
        //    fct.DefineConst("U", uLB_ref * ((N + 3) / (N + 1)));
        //    fct.DefineConst("NplusOne", N + 1.0);

        //SPtr<BC> inflowConcreteBC(new VelocityBC(false, false, true, fct, 0, BCFunction::INFCONST));
        //inflowConcreteBC->setBCStrategy(SPtr<BCStrategy>(new VelocityBCStrategy()));
        SPtr<BC> inflowConcreteBC(new MultiphaseVelocityBC(false, false, true, fct, phiH, 0, BCFunction::INFCONST));
        inflowConcreteBC->setBCStrategy(SPtr<BCStrategy>(new MultiphaseVelocityBCStrategy()));

        // air inflow boundary condition
        //  Štigler, J. (2014). Analytical velocity profile in tube for laminar and turbulent flow. Engineering
        //  Mechanics, 21(6), 371-379.
        // SPtr<LBMUnitConverter> unitsAir = std::make_shared<LBMUnitConverter>(d_part, LBMUnitConverter::AIR_20C, d_part / dx);
        // SPtr<LBMUnitConverter> unitsAir = std::make_shared<LBMUnitConverter>(d_part, 1., 1.2041, d_part / dx, uLB_ref);
        // double V = 40;      // flow rate [m^3/h]
        // double D = 0.0166;  // air inlet diameter [m]
        // double R = D / 2.0; // radius [m]
        // double A = UbMath::PI * R * R;
        // double u = V / 3600 / A;
        // double uLB_ref = u * unitsAir->getFactorVelocityWToLb();
        //// double cx1 = -1.2788 + R;
        // double cx2 = 0.3803 + R;
        // double cx3 = 0.1517 + R;
        // double L = g_maxX1 - g_minX1;
        // double p_air = 7e5; // Pa = 7 Bar
        // double p1 = p_air;
        // double p2 = 0.0;
        // double mu = 17.2e-6; // Pa s, air 20° C
        // double N = R * R / 2 * mu * u * (p1 - p2) / L - 3;
        // if (myid == 0) VF_LOG_INFO("Air inflow velocity = {} m/s", u);
        // if (myid == 0) VF_LOG_INFO("Air inflow velocity = {} dx/dt", uLB_ref);
        //

        // double nu = mu / rhoConcrete;
        // double Re = d_part * u / nu;
        // if (myid == 0) VF_LOG_INFO("Re_air = {}", Re);

        // double nuLB = d_part * unitsAir->getFactorLentghWToLb() * uLB_ref / Re;
        // if (myid == 0) VF_LOG_INFO("nuLB_air = {}", nuLB);
        // nu_l_LB = nuLB;
  
        SPtr<LBMUnitConverter> unitsAir = std::make_shared<LBMUnitConverter>(d_part, 1., 1.2041, d_part / dx, uLB_ref);
        //SPtr<LBMUnitConverter> unitsAir = std::make_shared<LBMUnitConverter>(d_part, LBMUnitConverter::AIR_20C, d_part / dx);
        double V_air = 40./6.;  // flow rate [m^3/h] //10.
        double D_air = 0.00553;     // air inlet diameter [m]
        double R_air = D_air / 2.0; // radius [m]
        double A_air = UbMath::PI * (R_air * R_air);
        double u_air = V_air / 3600 / A_air;
        double uLB_air = u_air * unitsAir->getFactorVelocityWToLb();
        // double cx1 = -1.2788 + R;
        double cx2 = 0.385822 + R_air;
        double cx3 = 0.135562 + R_air;
        //double L_air = 0.00747;
        double p_air = 7e5; // Pa = 7 Bar
        //double p1 = p_air;
        //double p2 = 1e5;
        double mu_air = 17.2e-6; // Pa s, air 20° C
        double rho_air = 1.2041; // [kg/m^3]
        double Re_inlet = D_air * u_air * rho_air / mu_air;
        double lambda = 0.3164 / pow(Re_inlet, 0.25);
        double deltaP = (lambda / (2. * R_air)) * (rho_air * pow(u_air, 2) / 2.); // Darcy friction factor (Rohrreibungszahl)
        double N = pow(R_air, 2) / (2. * mu_air * u_air) * deltaP - 3.;
        // double N = R_air * R_air / 2 * mu_air * u_air * (p1 - p2) / L_air - 3;
        if (myid == 0) VF_LOG_INFO("Air inflow velocity = {} m/s", u_air);
        if (myid == 0) VF_LOG_INFO("Air inflow velocity = {} dx/dt", uLB_air);

        double nu_air = mu_air / rho_air;
        double Re_air = d_part * u_air / nu_air;
        if (myid == 0) VF_LOG_INFO("Air Re = {}", Re_air);

        double nuLB_air = nu_air * unitsAir->getFactorViscosityWToLb(); // d_part * unitsAir->getFactorLentghWToLb() * uLB_air / Re_air;
        if (myid == 0) VF_LOG_INFO("nuLB_air = {}", nuLB_air);
        nu_l_LB = nuLB_air;

        if (myid == 0) VF_LOG_INFO("nu_h = {}", nu_h_LB);
        if (myid == 0) VF_LOG_INFO("nu_l = {}", nu_l_LB);
        if (myid == 0) VF_LOG_INFO("sigma_LB = {}", sigma_LB);

        double p_air_LB = p_air * unitsAir->getFactorPressureWToLb();
        if (myid == 0) VF_LOG_INFO("p_air_LB = {}", p_air_LB);

        // mu::Parser fctVx1;
        ////fctVx1.SetExpr("U");
        ////fctVx1.DefineConst("U", uLB_air);
        // mu::Parser fctVx2;
        // fctVx2.SetExpr("U");
        // fctVx2.DefineConst("U", 0);
        // mu::Parser fctVx3;
        ////fctVx3.SetExpr("U");
        ////fctVx3.DefineConst("U", -uLB_air);

        double cx1 = 0;
        double alpha = 0;
        double gamma = 0;
        double U = uLB_air;// * ((N + 3.) / (N + 1.));

        mu::Parser fctVx1;
        //fctVx1.SetExpr("U*cos(alpha*_pi/180)*(1-(((((x1-x0)^2+(x2-y0)^2+(x3-z0)^2)^0.5)/R)^NplusOne))");
        //fctVx1.SetExpr("U*cos(alpha*_pi/180)*(1-(((((x1-x0)^2+(x2-y0)^2+(x3-z0)^2)^0.5)/R)))");
        // fctVx1.SetExpr("(((x1-x0)^2+(x2-y0)^2+(x3-z0)^2)^0.5/R)^NplusOne");
        fctVx1.SetExpr("U*cos(alpha*_pi/180)");
        fctVx1.DefineConst("x0", cx1);
        fctVx1.DefineConst("y0", cx2);
        fctVx1.DefineConst("z0", cx3);
        fctVx1.DefineConst("R", R_air);
        fctVx1.DefineConst("U", U); //* ((N + 3.) / (N + 1.)));
        fctVx1.DefineConst("NplusOne", N + 1.0);
        fctVx1.DefineConst("alpha", alpha);
        fctVx1.DefineConst("gamma", gamma);

        mu::Parser fctVx2;
        //fctVx2.SetExpr("U*sin(alpha*_pi/180)*(1-(((((x1-x0)^2+(x2-y0)^2+(x3-z0)^2)^0.5)/R)^NplusOne))");
        //fctVx2.SetExpr("U*sin(alpha*_pi/180)*(1-(((((x1-x0)^2+(x2-y0)^2+(x3-z0)^2)^0.5)/R)))");
        fctVx2.SetExpr("U*sin(alpha*_pi/180)");
        fctVx2.DefineConst("x0", cx1);
        fctVx2.DefineConst("x0", cx1);
        fctVx2.DefineConst("y0", cx2);
        fctVx2.DefineConst("z0", cx3);
        fctVx2.DefineConst("R", R_air);
        fctVx2.DefineConst("U", U); //* ((N + 3.) / (N + 1.)));
        fctVx2.DefineConst("NplusOne", N + 1.0);
        fctVx2.DefineConst("alpha", alpha);

        mu::Parser fctVx3;
        //fctVx3.SetExpr("U*(1-(((((x1-x0)^2+(x2-y0)^2+(x3-z0)^2)^0.5)/R)^NplusOne))");
        //fctVx3.SetExpr("U*cos(alpha*_pi/180)*(1-(((((x1-x0)^2+(x2-y0)^2+(x3-z0)^2)^0.5)/R)))");
        fctVx3.SetExpr("U");
        fctVx3.DefineConst("x0", cx1);
        fctVx3.DefineConst("y0", cx2);
        fctVx3.DefineConst("z0", cx3);
        fctVx3.DefineConst("R", R_air);
        fctVx3.DefineConst("U", -U); //* ((N + 3.) / (N + 1.)));
        fctVx3.DefineConst("NplusOne", N + 1.0);
        fctVx3.DefineConst("alpha", alpha);
        fctVx3.DefineConst("gamma", gamma);

        // SPtr<BC> inflowAirBC1(new VelocityBC(true, false, false, fct, 0, BCFunction::INFCONST));
        // inflowAirBC1->setBCStrategy(SPtr<BCStrategy>(new VelocityBCStrategy()));
        // t = U * sin(alpha * _pi / 180) * (1 - (((((x1 - x0) ^ 2 + (x2 - y0) ^ 2 + (x3 - z0) ^ 2) ^ 0.5) / R) ^ NplusOne));
        cx1 = -1.31416;
        cx2 = 0.388684;
        cx3 = 0.138177;
        alpha = 0;
        gamma = 225;
        fctVx1.DefineConst("x0", cx1);
        fctVx1.DefineConst("y0", cx2);
        fctVx1.DefineConst("z0", cx3);
        fctVx1.DefineConst("alpha", alpha);
        fctVx1.DefineConst("gamma", gamma);
        fctVx2.DefineConst("x0", cx1);
        fctVx2.DefineConst("y0", cx2);
        fctVx2.DefineConst("z0", cx3);
        fctVx2.DefineConst("alpha", alpha);
        fctVx2.DefineConst("gamma", gamma);
        fctVx3.DefineConst("x0", cx1);
        fctVx3.DefineConst("y0", cx2);
        fctVx3.DefineConst("z0", cx3);
        fctVx3.DefineConst("alpha", alpha);
        fctVx3.DefineConst("gamma", gamma);

        //SPtr<BC> inflowAirBC1(new VelocityBC(true, true, true, fctVx1, fctVx2, fctVx3, 0, BCFunction::INFCONST));
        //inflowAirBC1->setBCStrategy(SPtr<BCStrategy>(new VelocityBCStrategy()));        
        SPtr<BC> inflowAirBC1(new MultiphaseVelocityBC(true, true, true, fctVx1, fctVx2, fctVx3, phiL, 0, BCFunction::INFCONST));
        inflowAirBC1->setBCStrategy(SPtr<BCStrategy>(new MultiphaseVelocityBCStrategy()));

        fctVx1.DefineVar("x1", &cx1);
        fctVx1.DefineVar("x2", &cx2);
        fctVx1.DefineVar("x3", &cx3);
        fctVx2.DefineVar("x1", &cx1);
        fctVx2.DefineVar("x2", &cx2);
        fctVx2.DefineVar("x3", &cx3);
        fctVx3.DefineVar("x1", &cx1);
        fctVx3.DefineVar("x2", &cx2);
        fctVx3.DefineVar("x3", &cx3);

        if (myid == 0) {

            VF_LOG_INFO("fctVx1 = {}", fctVx1.Eval());
            VF_LOG_INFO("fctVx2 = {}", fctVx2.Eval());
            VF_LOG_INFO("fctVx3 = {}", fctVx3.Eval());
            VF_LOG_INFO("N = {}", N);
            VF_LOG_INFO("NplusOne = {}", N + 1.0);
            // return 0;
        }
        cx1 = -1.31303;
        cx2 = 0.377234;
        cx3 = 0.138174;
        alpha = 60;
        fctVx1.DefineConst("x0", cx1);
        fctVx1.DefineConst("y0", cx2);
        fctVx1.DefineConst("z0", cx3);
        fctVx1.DefineConst("alpha", alpha);
        fctVx2.DefineConst("x0", cx1);
        fctVx2.DefineConst("y0", cx2);
        fctVx2.DefineConst("z0", cx3);
        fctVx2.DefineConst("alpha", alpha);
        fctVx3.DefineConst("x0", cx1);
        fctVx3.DefineConst("y0", cx2);
        fctVx3.DefineConst("z0", cx3);
        //SPtr<BC> inflowAirBC2(new VelocityBC(true, true, true, fctVx1, fctVx2, fctVx3, 0, BCFunction::INFCONST));
        //inflowAirBC2->setBCStrategy(SPtr<BCStrategy>(new VelocityBCStrategy()));
        SPtr<BC> inflowAirBC2(new MultiphaseVelocityBC(true, true, true, fctVx1, fctVx2, fctVx3, phiL, 0, BCFunction::INFCONST));
        inflowAirBC2->setBCStrategy(SPtr<BCStrategy>(new MultiphaseVelocityBCStrategy()));

        cx1 = -1.2948374155694822;
        cx2 = 0.37733728717266285;
        cx3 = 0.13840460401111598;
        alpha = 120;
        fctVx1.DefineConst("x0", cx1);
        fctVx1.DefineConst("y0", cx2);
        fctVx1.DefineConst("z0", cx3);
        fctVx1.DefineConst("alpha", alpha);
        fctVx2.DefineConst("x0", cx1);
        fctVx2.DefineConst("y0", cx2);
        fctVx2.DefineConst("z0", cx3);
        fctVx2.DefineConst("alpha", alpha);
        fctVx3.DefineConst("x0", cx1);
        fctVx3.DefineConst("y0", cx2);
        fctVx3.DefineConst("z0", cx3);
        //SPtr<BC> inflowAirBC3(new VelocityBC(true, true, true, fctVx1, fctVx2, fctVx3, 0, BCFunction::INFCONST));
        //inflowAirBC3->setBCStrategy(SPtr<BCStrategy>(new VelocityBCStrategy()));
        SPtr<BC> inflowAirBC3(new MultiphaseVelocityBC(true, true, true, fctVx1, fctVx2, fctVx3, phiL, 0, BCFunction::INFCONST));
        inflowAirBC3->setBCStrategy(SPtr<BCStrategy>(new MultiphaseVelocityBCStrategy()));

        cx1 = -1.28847;
        cx2 = 0.3885;
        cx3 = 0.1385;
        alpha = 180;
        fctVx1.DefineConst("x0", cx1);
        fctVx1.DefineConst("y0", cx2);
        fctVx1.DefineConst("z0", cx3);
        fctVx1.DefineConst("alpha", alpha);
        fctVx2.DefineConst("x0", cx1);
        fctVx2.DefineConst("y0", cx2);
        fctVx2.DefineConst("z0", cx3);
        fctVx2.DefineConst("alpha", alpha);
        fctVx3.DefineConst("x0", cx1);
        fctVx3.DefineConst("y0", cx2);
        fctVx3.DefineConst("z0", cx3);
        //SPtr<BC> inflowAirBC4(new VelocityBC(true, true, true, fctVx1, fctVx2, fctVx3, 0, BCFunction::INFCONST));
        //inflowAirBC4->setBCStrategy(SPtr<BCStrategy>(new VelocityBCStrategy()));
        SPtr<BC> inflowAirBC4(new MultiphaseVelocityBC(true, true, true, fctVx1, fctVx2, fctVx3, phiL, 0, BCFunction::INFCONST));
        inflowAirBC4->setBCStrategy(SPtr<BCStrategy>(new MultiphaseVelocityBCStrategy()));

        cx1 = -1.294771417778694;
        cx2 = 0.399787947463142;
        cx3 = 0.1383429692754194;
        alpha = 240;
        fctVx1.DefineConst("x0", cx1);
        fctVx1.DefineConst("y0", cx2);
        fctVx1.DefineConst("z0", cx3);
        fctVx1.DefineConst("alpha", alpha);
        fctVx2.DefineConst("x0", cx1);
        fctVx2.DefineConst("y0", cx2);
        fctVx2.DefineConst("z0", cx3);
        fctVx2.DefineConst("alpha", alpha);
        fctVx3.DefineConst("x0", cx1);
        fctVx3.DefineConst("y0", cx2);
        fctVx3.DefineConst("z0", cx3);
        //SPtr<BC> inflowAirBC5(new VelocityBC(true, true, true, fctVx1, fctVx2, fctVx3, 0, BCFunction::INFCONST));
        //inflowAirBC5->setBCStrategy(SPtr<BCStrategy>(new VelocityBCStrategy()));
        SPtr<BC> inflowAirBC5(new MultiphaseVelocityBC(true, true, true, fctVx1, fctVx2, fctVx3, phiL, 0, BCFunction::INFCONST));
        inflowAirBC5->setBCStrategy(SPtr<BCStrategy>(new MultiphaseVelocityBCStrategy()));

        cx1 = -1.3077338898450492;
        cx2 = 0.3998516560596088;
        cx3 = 0.13843501416896437;
        alpha = 300;
        fctVx1.DefineConst("x0", cx1);
        fctVx1.DefineConst("y0", cx2);
        fctVx1.DefineConst("z0", cx3);
        fctVx1.DefineConst("alpha", alpha);
        fctVx2.DefineConst("x0", cx1);
        fctVx2.DefineConst("y0", cx2);
        fctVx2.DefineConst("z0", cx3);
        fctVx2.DefineConst("alpha", alpha);
        fctVx3.DefineConst("x0", cx1);
        fctVx3.DefineConst("y0", cx2);
        fctVx3.DefineConst("z0", cx3);
        //SPtr<BC> inflowAirBC6(new VelocityBC(true, true, true, fctVx1, fctVx2, fctVx3, 0, BCFunction::INFCONST));
        //inflowAirBC6->setBCStrategy(SPtr<BCStrategy>(new VelocityBCStrategy()));
        SPtr<BC> inflowAirBC6(new MultiphaseVelocityBC(true, true, true, fctVx1, fctVx2, fctVx3, phiL, 0, BCFunction::INFCONST));
        inflowAirBC6->setBCStrategy(SPtr<BCStrategy>(new MultiphaseVelocityBCStrategy()));

        // Pressure BC for air inlet
        // SPtr<BC> inflowAirBC1(new DensityBC(p_air_LB));
        // inflowAirBC1->setBCStrategy(SPtr<BCStrategy>(new MultiphasePressureBCStrategy()));

        SPtr<BC> outflowBC(new DensityBC(rhoLB));
        //outflowBC->setBCStrategy(SPtr<BCStrategy>(new NonEqDensityBCStrategy()));
        outflowBC->setBCStrategy(SPtr<BCStrategy>(new MultiphasePressureBCStrategy()));
        
        // SPtr<BC> outflowBC(new DensityBC(rhoLB));
        //outflowBC->setBCStrategy(SPtr<BCStrategy>(new MultiphaseNonReflectingOutflowBCStrategy()));
        //////////////////////////////////////////////////////////////////////////////////
        // BC visitor
        //BoundaryConditionsBlockVisitor bcVisitor;
        MultiphaseBoundaryConditionsBlockVisitor bcVisitor;
        bcVisitor.addBC(noSlipBC);
        bcVisitor.addBC(inflowConcreteBC);
        bcVisitor.addBC(inflowAirBC1);
        bcVisitor.addBC(outflowBC);

        // SPtr<LBMKernel> kernel   = make_shared<IBcumulantK17LBMKernel>();
         //SPtr<LBMKernel> kernel   = make_shared<CumulantK17LBMKernel>();
        // SPtr<LBMKernel> kernel = make_shared<MultiphaseTwoPhaseFieldsPressureFilterLBMKernel>();
        // SPtr<LBMKernel> kernel = make_shared<MultiphaseSimpleVelocityBaseExternalPressureLBMKernel>();
        SPtr<LBMKernel> kernel = make_shared<MultiphaseSharpInterfaceLBMKernel>();
        //SPtr<LBMKernel> kernel = make_shared<MultiphaseScaleDistributionLBMKernel>();
        //SPtr<LBMKernel> kernel = make_shared<IBcumulantK17LBMKernel>();
        //SPtr<LBMKernel> kernel = make_shared<IBsharpInterfaceLBMKernel>();

        kernel->setWithForcing(false);
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
        kernel->setSigma(sigma_LB);

        SPtr<BCSet> bcProc = make_shared<BCSet>();
        kernel->setBCSet(bcProc);

        SPtr<Grid3D> grid = make_shared<Grid3D>(comm);
        grid->setPeriodicX1(false);
        grid->setPeriodicX2(false);
        grid->setPeriodicX3(false);
        grid->setDeltaX(dx);
        grid->setBlockNX(blockNX[0], blockNX[1], blockNX[2]);
        grid->setGhostLayerWidth(2);

        //string geoPath = "/home/niikonst/NozzleGeo";

        //string outputPath = "/scratch/projects/nii00154/ShotcreteJet2";
        UbSystem::makeDirectory(outputPath);
        UbSystem::makeDirectory(outputPath + "/liggghts");

        // if (myid == 0) {
        //     stringstream logFilename;
        //     logFilename << outputPath + "/logfile" + UbSystem::toString(UbSystem::getTimeStamp()) + ".txt";
        //     UbLog::output_policy::setStream(logFilename.str());
        // }

        /////////////////////////////////////////////////////////////////////
        //LIGGGHTS things
        /////////////////////////////////////////////////////////////////////
        // string inFile1 = "d:/Projects/TRR277/Project/WP4/Config/in.nozzle";
        // // string inFile2 = "d:/Projects/VirtualFluids_LIGGGHTS_coupling/apps/cpu/LiggghtsApp/in2.lbdem";
        // MPI_Comm mpi_comm = *(MPI_Comm *)(comm->getNativeCommunicator());
        // LiggghtsCouplingWrapper wrapper(argv, mpi_comm);

        // double v_frac = 0.1;
        // double dt_phys = units->getFactorTimeLbToW();
        // int demSubsteps = 10;
        // double dt_dem = dt_phys / (double)demSubsteps;
         //int vtkSteps = 1000;
        // string demOutDir = outputPath + "/liggghts";

        // // wrapper.execCommand("echo none");

        // // wrapper.execFile((char*)inFile1.c_str());

        // //// set timestep and output directory
        // //////wrapper.setVariable("t_step", dt_dem);
        // //////wrapper.setVariable("dmp_stp", vtkSteps * demSubsteps);
        // //////wrapper.setVariable("dmp_dir", demOutDir);

        // //!!!!//wrapper.execFile((char *)inFile1.c_str());
        // //wrapper.runUpto(demSubsteps - 1);
        // // wrapper.runUpto(1000);

        // //LatticeDecomposition lDec((g_maxX1 - g_minX1) / dx, (g_maxX2 - g_minX2) / dx, (g_maxX3 - g_minX3) / dx, wrapper.lmp, grid);

        SPtr<UbScheduler> lScheduler = make_shared<UbScheduler>(1);
        // SPtr<LiggghtsCouplingSimulationObserver> lcSimulationObserver = make_shared<LiggghtsCouplingSimulationObserver>(grid, lScheduler, comm, wrapper, demSubsteps, unitsAir);
        // //SPtr<Grid3DVisitor> partVisitor = make_shared<LiggghtsPartitioningGridVisitor>(std::ceil((g_maxX1 - g_minX1) / dx), std::ceil((g_maxX2 - g_minX2) / dx), std::ceil((g_maxX3 - g_minX3) / dx), wrapper.lmp);
        // //SPtr<Grid3DVisitor> partVisitor = make_shared<LiggghtsPartitioningGridVisitor>(blockNX[0], blockNX[1], blockNX[2] * gridNZ, wrapper.lmp);
        
        // /////////////////////////////////////////////////////////////////////
        // /////////////////////////////////////////////////////////////////////
        
        SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, vf::lbm::dir::DIR_MMM, MetisPartitioner::KWAY));

        //////////////////////////////////////////////////////////////////////////
        // restart
        //double cpStep  = 15000; 
        //double cpStart = 10;
        SPtr<UbScheduler> rSch(new UbScheduler(cpStep, cpStart));
        SPtr<MPIIOMigrationSimulationObserver> rcp(new MPIIOMigrationSimulationObserver(grid, rSch, metisVisitor, outputPath, comm));
        rcp->setLBMKernel(kernel);
        rcp->setBCSet(bcProc);
        //////////////////////////////////////////////////////////////////////////

//if (newStart) {
        SPtr<GbObject3D> gridCube = make_shared<GbCuboid3D>(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3);
        if (myid == 0) GbSystem3D::writeGeoObject(gridCube.get(), outputPath + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

        GenBlocksGridVisitor genBlocks(gridCube);
        grid->accept(genBlocks);

        // geo
        //////////////////////////////////////////////////////////
        //int accuracy = Interactor3D::EDGES;
        ///////////////////////////////////
        SPtr<GbTriFaceMesh3D> meshNozzleAirDistributor = std::make_shared<GbTriFaceMesh3D>();
        if (myid == 0) UBLOG(logINFO, "Read meshNozzleAirDistributor:start");
        meshNozzleAirDistributor->readMeshFromSTLFileASCII(geoPath + "/01_Nozzle_Air_Distributor.stl", false);
        if (myid == 0) UBLOG(logINFO, "Read meshNozzleAirDistributor:end");
        if (myid == 0) GbSystem3D::writeGeoObject(meshNozzleAirDistributor.get(), outputPath + "/geo/meshNozzleAirDistributor", WbWriterVtkXmlBinary::getInstance());
        SPtr<Interactor3D> intrNozzleAirDistributor = std::make_shared<D3Q27TriFaceMeshInteractor>(meshNozzleAirDistributor, grid, noSlipBC, Interactor3D::SOLID, Interactor3D::POINTS);
        /////////////////////////////////////////////////////////////
        //SPtr<GbTriFaceMesh3D> meshNozzleAirInlet = std::make_shared<GbTriFaceMesh3D>();
        //if (myid == 0) UBLOG(logINFO, "Read meshNozzleAirInlet:start");
        //meshNozzleAirInlet->readMeshFromSTLFileASCII(geoPath + "/02_Nozzle_Air_Inlet.stl", false);
        //if (myid == 0) UBLOG(logINFO, "Read meshNozzleAirInlet:end");
        //if (myid == 0) GbSystem3D::writeGeoObject(meshNozzleAirInlet.get(), outputPath + "/geo/meshNozzleAirInlet", WbWriterVtkXmlBinary::getInstance());
        //SPtr<Interactor3D> intrNozzleAirInlet = std::make_shared<D3Q27TriFaceMeshInteractor>(meshNozzleAirInlet, grid, noSlipBC, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy);
        /////////////////////////////////////////////////////////////
        //SPtr<GbTriFaceMesh3D> meshNozzleSpacer = std::make_shared<GbTriFaceMesh3D>();
        //if (myid == 0) UBLOG(logINFO, "Read meshNozzleSpacer:start");
        //meshNozzleSpacer->readMeshFromSTLFileASCII(geoPath + "/03_Nozzle_Spacer.stl", true);
        //if (myid == 0) UBLOG(logINFO, "Read meshNozzleSpacer:end");
        //if (myid == 0) GbSystem3D::writeGeoObject(meshNozzleSpacer.get(), outputPath + "/geo/meshNozzleSpacer", WbWriterVtkXmlBinary::getInstance());
        //SPtr<Interactor3D> intrNozzleSpacer = std::make_shared<D3Q27TriFaceMeshInteractor>(meshNozzleSpacer, grid, noSlipBC, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy);
        /////////////////////////////////////////////////////////////
        //SPtr<GbTriFaceMesh3D> meshNozzleAccDistributor = std::make_shared<GbTriFaceMesh3D>();
        //if (myid == 0) UBLOG(logINFO, "Read meshNozzleAccDistributor:start");
        //meshNozzleAccDistributor->readMeshFromSTLFileASCII(geoPath + "/04_Nozzle_Acc_Distributor.stl", false);
        //if (myid == 0) UBLOG(logINFO, "Read meshNozzleAccDistributor:end");
        //if (myid == 0) GbSystem3D::writeGeoObject(meshNozzleAccDistributor.get(), outputPath + "/geo/meshNozzleAccDistributor", WbWriterVtkXmlBinary::getInstance());
        //SPtr<Interactor3D> intrNozzleAccDistributor = std::make_shared<D3Q27TriFaceMeshInteractor>(meshNozzleAccDistributor, grid, noSlipBC, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy);
        /////////////////////////////////////////////////////////////
        //SPtr<GbTriFaceMesh3D> meshNozzleAccInlet = std::make_shared<GbTriFaceMesh3D>();
        //if (myid == 0) UBLOG(logINFO, "Read meshNozzleAccInlet:start");
        //meshNozzleAccInlet->readMeshFromSTLFileASCII(geoPath + "/05_Nozzle_Acc_Inlet.stl", false);
        //if (myid == 0) UBLOG(logINFO, "Read meshNozzleAccInlet:end");
        //if (myid == 0) GbSystem3D::writeGeoObject(meshNozzleAccInlet.get(), outputPath + "/geo/meshNozzleAccInlet", WbWriterVtkXmlBinary::getInstance());
        //SPtr<Interactor3D> intrNozzleAccInlet = std::make_shared<D3Q27TriFaceMeshInteractor>(meshNozzleAccInlet, grid, noSlipBC, Interactor3D::SOLID, (Interactor3D::Accuracy)accuracy);
        /////////////////////////////////////////////////////////////
        //SPtr<GbTriFaceMesh3D> meshNozzleVolcanNozzle1 = std::make_shared<GbTriFaceMesh3D>();
        //if (myid == 0) UBLOG(logINFO, "Read meshNozzleVolcanNozzle1:start");
        //meshNozzleVolcanNozzle1->readMeshFromSTLFileBinary(geoPath + "/06_1_Nozzle_Volcan_Nozzle.stl", true);
        //if (myid == 0) UBLOG(logINFO, "Read meshNozzleVolcanNozzle1:end");
        //if (myid == 0) GbSystem3D::writeGeoObject(meshNozzleVolcanNozzle1.get(), outputPath + "/geo/meshNozzleVolcanNozzle1", WbWriterVtkXmlBinary::getInstance());
        //SPtr<Interactor3D> intrNozzleVolcanNozzle1 = std::make_shared<D3Q27TriFaceMeshInteractor>(meshNozzleVolcanNozzle1, grid, noSlipBC, Interactor3D::SOLID, Interactor3D::EDGES);
        ///////////////////////////////////////////////////////////
        SPtr<GbTriFaceMesh3D> meshNozzleVolcanNozzle2 = std::make_shared<GbTriFaceMesh3D>();
        if (myid == 0) UBLOG(logINFO, "Read meshNozzleVolcanNozzle2:start");
        //meshNozzleVolcanNozzle2->readMeshFromSTLFileBinary(geoPath + "/06_2_Nozzle_Volcan_Nozzle.stl", true);
        meshNozzleVolcanNozzle2->readMeshFromSTLFileBinary(geoPath + "/Nozzle_Volcan_Nozzle_Shift.stl", true);
        if (myid == 0) UBLOG(logINFO, "Read meshNozzleVolcanNozzle2:end");
        if (myid == 0) GbSystem3D::writeGeoObject(meshNozzleVolcanNozzle2.get(), outputPath + "/geo/meshNozzleVolcanNozzle2", WbWriterVtkXmlBinary::getInstance());
        SPtr<Interactor3D> intrNozzleVolcanNozzle2 = std::make_shared<D3Q27TriFaceMeshInteractor>(meshNozzleVolcanNozzle2, grid, noSlipBC, Interactor3D::SOLID, Interactor3D::POINTS);
        ///////////////////////////////////////////////////////////
        // box
        SPtr<D3Q27Interactor> intrBox = SPtr<D3Q27Interactor>(new D3Q27Interactor(gridCube, grid, noSlipBC, Interactor3D::INVERSESOLID));
        ///////////////////////////////////////////////////////////
        // inflow
        //GbCylinder3DPtr geoInflow(new GbCylinder3D(-1.30181 + 0.0005, 0.390872 - 0.00229, 0.20105, -1.30181 + 0.0005, 0.390872 - 0.00229, 0.23, 0.013));
        GbCylinder3DPtr geoInflow(new GbCylinder3D(-1.30181 + 0.0005, 0.390872 - 0.00229, g_maxX3 - 2.0 * dx, -1.30181 + 0.0005, 0.390872 - 0.00229, g_maxX3+2.0*dx, 0.013));
        if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), outputPath + "/geo/geoInflow", WbWriterVtkXmlBinary::getInstance());
        SPtr<D3Q27Interactor> intrInflow = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, inflowConcreteBC, Interactor3D::SOLID));
        ///////////////////////////////////////////////////////////
        // outflow
        //GbCylinder3DPtr geoOutflow(new GbCylinder3D(-1.30181 + 0.0005, 0.390872 - 0.00229, -0.22, -1.30181 + 0.0005, 0.390872 - 0.00229, -0.21, 0.013));
        //GbCylinder3DPtr geoOutflow(new GbCylinder3D(-1.30181 + 0.0005, 0.390872 - 0.00229, -0.426, -1.30181 + 0.0005, 0.390872 - 0.00229, -0.415, 0.013));
        GbCylinder3DPtr geoOutflow(new GbCylinder3D(-1.30181 + 0.0005, 0.390872 - 0.00229, g_minX3, -1.30181 + 0.0005, 0.390872 - 0.00229, g_minX3+2.*dx, 0.013));
        if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow.get(), outputPath + "/geo/geoOutflow", WbWriterVtkXmlBinary::getInstance());
        SPtr<D3Q27Interactor> intrOutflow = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, outflowBC, Interactor3D::SOLID));
        ///////////////////////////////////////////////////////////
        // SPtr<GbTriFaceMesh3D> geoAirInlet = std::make_shared<GbTriFaceMesh3D>();
        // if (myid == 0) UBLOG(logINFO, "Read Air_Inlet:start");
        // geoAirInlet->readMeshFromSTLFileASCII(geoPath + "/Air_Inlet.stl", true);
        // if (myid == 0) UBLOG(logINFO, "Read Air_Inlet:end");
        // if (myid == 0) GbSystem3D::writeGeoObject(geoAirInlet.get(), outputPath + "/geo/geoAirInlet", WbWriterVtkXmlBinary::getInstance());
        // SPtr<Interactor3D> intrAirInlet = std::make_shared<D3Q27TriFaceMeshInteractor>(geoAirInlet, grid, inflowAirBC1, Interactor3D::SOLID, Interactor3D::EDGES);
        /////////////////////////////////////////////////////////////
        // Fluid area
        //GbCylinder3DPtr geoFluidArea(new GbCylinder3D(-1.30181 + 0.0005, 0.390872 - 0.00229, g_minX3, -1.30181 + 0.0005, 0.390872 - 0.00229, g_maxX3, 0.013));
        GbCylinder3DPtr geoFluidArea(new GbCylinder3D(-1.30181 + 0.0005, 0.390872 - 0.00229, g_maxX3_box, -1.30181 + 0.0005, 0.390872 - 0.00229, g_maxX3, 0.013));
        if (myid == 0) GbSystem3D::writeGeoObject(geoFluidArea.get(), outputPath + "/geo/geoFluidArea", WbWriterVtkXmlBinary::getInstance());
        SPtr<D3Q27Interactor> intrFluidArea = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoFluidArea, grid, noSlipBC, Interactor3D::INVERSESOLID));

        //SPtr<GbTriFaceMesh3D> meshFluidArea = std::make_shared<GbTriFaceMesh3D>();
        //if (myid == 0) UBLOG(logINFO, "Read geoFluidArea:start");
        //meshFluidArea->readMeshFromSTLFileBinary(geoPath + "/FluidArea.stl", true);
        //if (myid == 0) UBLOG(logINFO, "Read geoFluidArea:end");
        //if (myid == 0) GbSystem3D::writeGeoObject(meshFluidArea.get(), outputPath + "/geo/meshFluidArea", WbWriterVtkXmlBinary::getInstance());
        //SPtr<Interactor3D> intrFluidArea = std::make_shared<D3Q27TriFaceMeshInteractor>(meshFluidArea, grid, noSlipBC, Interactor3D::SOLID, Interactor3D::EDGES);

        ///////////////////////////////////////////////////////////
        SPtr<GbTriFaceMesh3D> meshInflowPipe = std::make_shared<GbTriFaceMesh3D>();
        if (myid == 0) UBLOG(logINFO, "Read geoFluidArea:start");
        meshInflowPipe->readMeshFromSTLFileBinary(geoPath + "/InflowPipe.stl", true);
        if (myid == 0) UBLOG(logINFO, "Read geoFluidArea:end");
        if (myid == 0) GbSystem3D::writeGeoObject(meshInflowPipe.get(), outputPath + "/geo/meshInflowPipe", WbWriterVtkXmlBinary::getInstance());
        SPtr<Interactor3D> intrInflowPipe = std::make_shared<D3Q27TriFaceMeshInteractor>(meshInflowPipe, grid, noSlipBC, Interactor3D::SOLID, Interactor3D::EDGES);
        ///////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////
        SPtr<GbTriFaceMesh3D> meshInflowPipe2 = std::make_shared<GbTriFaceMesh3D>();
        if (myid == 0) UBLOG(logINFO, "Read geoFluidArea:start");
        meshInflowPipe2->readMeshFromSTLFileBinary(geoPath + "/LongTube2.stl", true);
        if (myid == 0) UBLOG(logINFO, "Read geoFluidArea:end");
        if (myid == 0) GbSystem3D::writeGeoObject(meshInflowPipe2.get(), outputPath + "/geo/LongTube", WbWriterVtkXmlBinary::getInstance());
        SPtr<Interactor3D> intrInflowPipe2 = std::make_shared<D3Q27TriFaceMeshInteractor>(meshInflowPipe2, grid, noSlipBC, Interactor3D::SOLID, Interactor3D::EDGES);
        ///////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////
        //outflows 
        //////////////////////////////////////////////////////////
        SPtr<GbObject3D> geoOutflow1 = make_shared<GbCuboid3D>(g_minX1 - 2.0 * dx, g_minX2 - 2.0 * dx, g_minX3 - 2.0 * dx, g_maxX1 + 2.0 * dx, g_maxX2 + 2.0 * dx, g_minX3);
        if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow1.get(), outputPath + "/geo/geoOutflow1", WbWriterVtkXmlBinary::getInstance());
        SPtr<D3Q27Interactor> intrOutflow1 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow1, grid, outflowBC, Interactor3D::SOLID));

        SPtr<GbObject3D> geoOutflow2 = make_shared<GbCuboid3D>(g_minX1 - 2.0 * dx, g_minX2 - 2.0 * dx, g_minX3 - 2.0 * dx, g_minX1, g_maxX2 + 2.0 * dx, g_maxX3 + 2.0 * dx);
        if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow2.get(), outputPath + "/geo/geoOutflow2", WbWriterVtkXmlBinary::getInstance());
        SPtr<D3Q27Interactor> intrOutflow2 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow2, grid, outflowBC, Interactor3D::SOLID));

        SPtr<GbObject3D> geoOutflow3 = make_shared<GbCuboid3D>(g_maxX1, g_minX2 - 2.0 * dx, g_minX3 - 2.0 * dx, g_maxX1 + 2.0 * dx, g_maxX2 + 2.0 * dx, g_maxX3 + 2.0 * dx);
        if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow3.get(), outputPath + "/geo/geoOutflow3", WbWriterVtkXmlBinary::getInstance());
        SPtr<D3Q27Interactor> intrOutflow3 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow3, grid, outflowBC, Interactor3D::SOLID));

        SPtr<GbObject3D> geoOutflow4 = make_shared<GbCuboid3D>(g_minX1 - 2.0 * dx, g_minX2 - 2.0 * dx, g_minX3 - 2.0 * dx, g_maxX1 + 2.0 * dx, g_minX2, g_maxX3 + 2.0 * dx);
        if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow4.get(), outputPath + "/geo/geoOutflow4", WbWriterVtkXmlBinary::getInstance());
        SPtr<D3Q27Interactor> intrOutflow4 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow4, grid, outflowBC, Interactor3D::SOLID));

        SPtr<GbObject3D> geoOutflow5 = make_shared<GbCuboid3D>(g_minX1 - 2.0 * dx, g_maxX2, g_minX3 - 2.0 * dx, g_maxX1 + 2.0 * dx, g_maxX2 + 2.0 * dx, g_maxX3 + 2.0 * dx);
        if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow5.get(), outputPath + "/geo/geoOutflow5", WbWriterVtkXmlBinary::getInstance());
        SPtr<D3Q27Interactor> intrOutflow5 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow5, grid, outflowBC, Interactor3D::SOLID));

        //SPtr<GbObject3D> geoOutflow6 = make_shared<GbCuboid3D>(g_minX1, g_minX2, g_maxX3_box, g_maxX1, g_maxX2, g_maxX3_box + 2.0 * dx);
        //if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow6.get(), outputPath + "/geo/geoOutflow6", WbWriterVtkXmlBinary::getInstance());
        //SPtr<D3Q27Interactor> intrOutflow6 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow6, grid, outflowBC, Interactor3D::SOLID));

        SPtr<GbTriFaceMesh3D> geoOutflow6 = std::make_shared<GbTriFaceMesh3D>();
        if (myid == 0) UBLOG(logINFO, "Read geoOutflow6:start");
        geoOutflow6->readMeshFromSTLFileBinary(geoPath + "/OutflowTop.stl", true);
        if (myid == 0) UBLOG(logINFO, "Read geoOutflow6:end");
        if (myid == 0) GbSystem3D::writeGeoObject(geoOutflow6.get(), outputPath + "/geo/geoOutflow6", WbWriterVtkXmlBinary::getInstance());
        SPtr<Interactor3D> intrOutflow6 = std::make_shared<D3Q27TriFaceMeshInteractor>(geoOutflow6, grid, outflowBC, Interactor3D::SOLID, Interactor3D::POINTS);

        ///////////////////////////////////////////////////////////
        GbCylinder3DPtr geoAirInflow(new GbCylinder3D(-1.31431 - 0.0005, 0.388587, 0.1383275, -1.31431, 0.388587, 0.1383275, 0.002765));
        if (myid == 0) GbSystem3D::writeGeoObject(geoAirInflow.get(), outputPath + "/geo/geoAirInlet", WbWriterVtkXmlBinary::getInstance());
        SPtr<Interactor3D> intrAirInflow = std::make_shared<D3Q27Interactor>(geoAirInflow, grid, inflowAirBC1, Interactor3D::SOLID, Interactor3D::EDGES);
        ///////////////////////////////////////////////////////////
        SPtr<GbTriFaceMesh3D> meshAirInlet1 = std::make_shared<GbTriFaceMesh3D>();
        if (myid == 0) UBLOG(logINFO, "Read meshAirInlet1:start");
        meshAirInlet1->readMeshFromSTLFileASCII(geoPath + "/Air_Inlet_1.stl", true);
        if (myid == 0) UBLOG(logINFO, "Read meshAirInlet1:end");
        if (myid == 0) GbSystem3D::writeGeoObject(meshAirInlet1.get(), outputPath + "/geo/meshAirInlet1", WbWriterVtkXmlBinary::getInstance());
        SPtr<Interactor3D> intAirInlet1 = std::make_shared<D3Q27TriFaceMeshInteractor>(meshAirInlet1, grid, inflowAirBC1, Interactor3D::SOLID, Interactor3D::POINTS);
        ///////////////////////////////////////////////////////////
        SPtr<GbTriFaceMesh3D> meshAirInlet2 = std::make_shared<GbTriFaceMesh3D>();
        if (myid == 0) UBLOG(logINFO, "Read meshAirInlet2:start");
        meshAirInlet2->readMeshFromSTLFileASCII(geoPath + "/Air_Inlet_2.stl", true);
        if (myid == 0) UBLOG(logINFO, "Read meshAirInlet2:end");
        if (myid == 0) GbSystem3D::writeGeoObject(meshAirInlet2.get(), outputPath + "/geo/meshAirInlet2", WbWriterVtkXmlBinary::getInstance());
        SPtr<Interactor3D> intAirInlet2 = std::make_shared<D3Q27TriFaceMeshInteractor>(meshAirInlet2, grid, inflowAirBC2, Interactor3D::SOLID, Interactor3D::POINTS);
        ///////////////////////////////////////////////////////////
        SPtr<GbTriFaceMesh3D> meshAirInlet3 = std::make_shared<GbTriFaceMesh3D>();
        if (myid == 0) UBLOG(logINFO, "Read meshAirInlet3:start");
        meshAirInlet3->readMeshFromSTLFileASCII(geoPath + "/Air_Inlet_3.stl", true);
        if (myid == 0) UBLOG(logINFO, "Read meshAirInlet3:end");
        if (myid == 0) GbSystem3D::writeGeoObject(meshAirInlet3.get(), outputPath + "/geo/meshAirInlet3", WbWriterVtkXmlBinary::getInstance());
        SPtr<Interactor3D> intAirInlet3 = std::make_shared<D3Q27TriFaceMeshInteractor>(meshAirInlet3, grid, inflowAirBC3, Interactor3D::SOLID, Interactor3D::POINTS);
        ///////////////////////////////////////////////////////////
        SPtr<GbTriFaceMesh3D> meshAirInlet4 = std::make_shared<GbTriFaceMesh3D>();
        if (myid == 0) UBLOG(logINFO, "Read meshAirInlet4:start");
        meshAirInlet4->readMeshFromSTLFileASCII(geoPath + "/Air_Inlet_4.stl", true);
        if (myid == 0) UBLOG(logINFO, "Read meshAirInlet4:end");
        if (myid == 0) GbSystem3D::writeGeoObject(meshAirInlet4.get(), outputPath + "/geo/meshAirInlet4", WbWriterVtkXmlBinary::getInstance());
        SPtr<Interactor3D> intAirInlet4 = std::make_shared<D3Q27TriFaceMeshInteractor>(meshAirInlet4, grid, inflowAirBC4, Interactor3D::SOLID, Interactor3D::POINTS);
        ///////////////////////////////////////////////////////////
        SPtr<GbTriFaceMesh3D> meshAirInlet5 = std::make_shared<GbTriFaceMesh3D>();
        if (myid == 0) UBLOG(logINFO, "Read meshAirInlet5:start");
        meshAirInlet5->readMeshFromSTLFileASCII(geoPath + "/Air_Inlet_5.stl", true);
        if (myid == 0) UBLOG(logINFO, "Read meshAirInlet5:end");
        if (myid == 0) GbSystem3D::writeGeoObject(meshAirInlet5.get(), outputPath + "/geo/meshAirInlet5", WbWriterVtkXmlBinary::getInstance());
        SPtr<Interactor3D> intAirInlet5 = std::make_shared<D3Q27TriFaceMeshInteractor>(meshAirInlet5, grid, inflowAirBC5, Interactor3D::SOLID, Interactor3D::POINTS);
        ///////////////////////////////////////////////////////////
        SPtr<GbTriFaceMesh3D> meshAirInlet6 = std::make_shared<GbTriFaceMesh3D>();
        if (myid == 0) UBLOG(logINFO, "Read meshAirInlet6:start");
        meshAirInlet6->readMeshFromSTLFileASCII(geoPath + "/Air_Inlet_6.stl", true);
        if (myid == 0) UBLOG(logINFO, "Read meshAirInlet6:end");
        if (myid == 0) GbSystem3D::writeGeoObject(meshAirInlet6.get(), outputPath + "/geo/meshAirInlet6", WbWriterVtkXmlBinary::getInstance());
        SPtr<Interactor3D> intAirInlet6 = std::make_shared<D3Q27TriFaceMeshInteractor>(meshAirInlet6, grid, inflowAirBC6, Interactor3D::SOLID, Interactor3D::POINTS);
        ///////////////////////////////////////////////////////////
        
        SPtr<GbObject3D> geoBox1 = make_shared<GbCuboid3D>(g_minX1, g_minX2, g_maxX3_box, g_minX1+12.0*blockNX[0]*dx, g_maxX2, g_maxX3);
        if (myid == 0) GbSystem3D::writeGeoObject(geoBox1.get(), outputPath + "/geo/geoBox1", WbWriterVtkXmlBinary::getInstance());
        SPtr<D3Q27Interactor> intrGeoBox1 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoBox1, grid, outflowBC, Interactor3D::SOLID));    

        SPtr<GbObject3D> geoBox2 = make_shared<GbCuboid3D>(g_minX1 + 14.0 * blockNX[0] * dx, g_minX2, g_maxX3_box, g_maxX1, g_maxX2, g_maxX3);
        if (myid == 0) GbSystem3D::writeGeoObject(geoBox2.get(), outputPath + "/geo/geoBox2", WbWriterVtkXmlBinary::getInstance());
        SPtr<D3Q27Interactor> intrGeoBox2 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoBox2, grid, outflowBC, Interactor3D::SOLID));

        SPtr<GbObject3D> geoBox3 = make_shared<GbCuboid3D>(g_minX1, g_minX2, g_maxX3_box, g_maxX1, g_minX2 + 12.0 * blockNX[0] * dx, g_maxX3);
        if (myid == 0) GbSystem3D::writeGeoObject(geoBox3.get(), outputPath + "/geo/geoBox3", WbWriterVtkXmlBinary::getInstance());
        SPtr<D3Q27Interactor> intrGeoBox3 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoBox3, grid, outflowBC, Interactor3D::SOLID));

        SPtr<GbObject3D> geoBox4 = make_shared<GbCuboid3D>(g_minX1, g_minX2 + 14.0 * blockNX[0] * dx, g_maxX3_box, g_maxX1, g_maxX2, g_maxX3);
        if (myid == 0) GbSystem3D::writeGeoObject(geoBox4.get(), outputPath + "/geo/geoBox4", WbWriterVtkXmlBinary::getInstance());
        SPtr<D3Q27Interactor> intrGeoBox4 = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoBox4, grid, outflowBC, Interactor3D::SOLID));
        
        InteractorsHelper intHelper1(grid, metisVisitor, true);
        //intHelper1.addInteractor(intrFluidArea);
        intHelper1.addInteractor(intrGeoBox1);
        intHelper1.addInteractor(intrGeoBox2);
        intHelper1.addInteractor(intrGeoBox3);
        intHelper1.addInteractor(intrGeoBox4);
        intHelper1.selectBlocks();

        //MultiphaseSetKernelBlockVisitor kernelVisitor(kernel, nu_h_LB, nu_l_LB, 1e9, 1);
        //grid->accept(kernelVisitor);

        //intHelper1.setBC();

        //SPtr<GbObject3D> gridCube2 = make_shared<GbCuboid3D>(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3_box);
        //if (myid == 0) GbSystem3D::writeGeoObject(gridCube2.get(), outputPath + "/geo/gridCube2", WbWriterVtkXmlBinary::getInstance());
        //GenBlocksGridVisitor genBlocks2(gridCube2);
        //grid->accept(genBlocks2);

        MultiphaseSetKernelBlockVisitor kernelVisitor2(kernel, nu_h_LB, nu_l_LB, 1e9, 1 ); // ,MultiphaseSetKernelBlockVisitor::AddKernel);
         grid->accept(kernelVisitor2);

        vector<SPtr<Block3D>> blocks;
         grid->getBlocksByCuboid(0, g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3, blocks);

        for (auto block : blocks) {
            if (block) {
                block->setActive(true);
                SPtr<BCArray3D> bcArray = block->getKernel()->getBCSet()->getBCArray();

                int minX1 = 0;
                int minX2 = 0;
                int minX3 = 0;

                int maxX1 = (int)(bcArray->getNX1()) - 1;
                int maxX2 = (int)(bcArray->getNX2()) - 1;
                int maxX3 = (int)(bcArray->getNX3()) - 1;

                for (int ix3 = minX3; ix3 <= maxX3; ix3++) {
                    for (int ix2 = minX2; ix2 <= maxX2; ix2++) {
                        for (int ix1 = minX1; ix1 <= maxX1; ix1++) {
                            bcArray->setFluid(ix1, ix2, ix3);
                        }
                    }
                }
            }
            }

        if (!newStart)
        {
            rcp->readBlocks((int)restartStep);
            grid->accept(metisVisitor);
            rcp->readDataSet((int)restartStep);
            grid->setTimeStep(restartStep);
        }

        InteractorsHelper intHelper2(grid, metisVisitor, false);
        intHelper2.addInteractor(intrInflowPipe);
        intHelper2.addInteractor(intrInflowPipe2);
        intHelper2.addInteractor(intrNozzleAirDistributor);
        //intHelper2.addInteractor(intrFluidArea);
        intHelper2.addInteractor(intrNozzleVolcanNozzle2);
        // intHelper.addInteractor(intrBox);
        intHelper2.addInteractor(intrInflow);
        //// intHelper.addInteractor(intrAirInflow);
        intHelper2.addInteractor(intAirInlet1);
        intHelper2.addInteractor(intAirInlet2);
        intHelper2.addInteractor(intAirInlet3);
        intHelper2.addInteractor(intAirInlet4);
        intHelper2.addInteractor(intAirInlet5);
        intHelper2.addInteractor(intAirInlet6);
        intHelper2.addInteractor(intrOutflow1);
        intHelper2.addInteractor(intrOutflow2);
        intHelper2.addInteractor(intrOutflow3);
        intHelper2.addInteractor(intrOutflow4);
        intHelper2.addInteractor(intrOutflow5);
        intHelper2.addInteractor(intrOutflow6);

        // intHelper.addInteractor(intrNozzleAirDistributor);
        // intHelper.addInteractor(intrNozzleAirInlet);
        // intHelper.addInteractor(intrNozzleSpacer);
        // intHelper.addInteractor(intrNozzleAccDistributor);
        // intHelper.addInteractor(intrNozzleAccInlet);
        // intHelper.addInteractor(intrNozzleVolcanNozzle1);

        intHelper2.selectBlocks();

         // if (myid == 0) UBLOG(logINFO, Utilities::toString(grid, comm->getNumberOfProcesses()));

        //SetKernelBlockVisitor kernelVisitor(kernel, nu_l_LB, comm->getNumberOfProcesses());


        intHelper2.setBC();



        //GenBlocksGridVisitor genBlocks2(gridCube2);
        //grid->accept(genBlocks2);

        //grid->accept(metisVisitor);

        //MultiphaseSetKernelBlockVisitor kernelVisitor2(kernel, nu_h_LB, nu_l_LB, 1e9, 1, MultiphaseSetKernelBlockVisitor::AddKernel);
        //grid->accept(kernelVisitor2);

        //SetBcBlocksBlockVisitor v1(intrOutflow1);
        //grid->accept(v1);
        //intrOutflow1->initInteractor();


        //SetBcBlocksBlockVisitor v2(intrOutflow2);
        //grid->accept(v2);
        //intrOutflow2->initInteractor();

        //SetBcBlocksBlockVisitor v3(intrOutflow3);
        //grid->accept(v3);
        //intrOutflow3->initInteractor();

        //SetBcBlocksBlockVisitor v4(intrOutflow4);
        //grid->accept(v4);
        //intrOutflow4->initInteractor();

        //SetBcBlocksBlockVisitor v5(intrOutflow5);
        //grid->accept(v5);
        //intrOutflow5->initInteractor();

        //SetBcBlocksBlockVisitor v6(intrOutflow6);
        //grid->accept(v6);
        //intrOutflow6->initInteractor();

        //SetBcBlocksBlockVisitor v7(intrNozzleVolcanNozzle2);
        //grid->accept(v7);
        //intrNozzleVolcanNozzle2->initInteractor();


       SPtr<SimulationObserver> ppblocks = make_shared<WriteBlocksSimulationObserver>(grid, SPtr<UbScheduler>(new UbScheduler(1)), outputPath, WbWriterVtkXmlBinary::getInstance(), comm);
        ppblocks->update(0);
        ppblocks.reset();

        // boundary conditions grid
        {
            SPtr<UbScheduler> geoSch(new UbScheduler(1));
            SPtr<WriteBoundaryConditionsSimulationObserver> ppgeo(new WriteBoundaryConditionsSimulationObserver(grid, geoSch, outputPath, WbWriterVtkXmlBinary::getInstance(), comm));
            ppgeo->update(0);
            ppgeo.reset();
        }



        if (newStart)
        {
            double x1c = -1.31431 + R;
            double x2c = 0.375582 + R;
            double Ri = 5;
            double x3c = 0.136 + Ri;

            mu::Parser fct1;
            // fct1.SetExpr(" 0.5 - 0.5 * tanh(2 * (sqrt((x1 - x1c) ^ 2 + (x2 - x2c) ^ 2 + (x3 - x3c) ^ 2) - radius) / interfaceThickness)");
            fct1.SetExpr(" 0.5 - 0.5 * tanh(2 * (sqrt((x1 - x1c) ^ 2 + (x2 - x2c) ^ 2 + (x3 - x3c) ^ 2) - radius) / interfaceThickness)");
            fct1.DefineConst("x1c", x1c);
            fct1.DefineConst("x2c", x2c);
            fct1.DefineConst("x3c", x3c);
            fct1.DefineConst("radius", Ri);
            fct1.DefineConst("interfaceThickness", interfaceThickness * dx);

            MultiphaseVelocityFormInitDistributionsBlockVisitor initVisitor;
            initVisitor.setPhi(fct1);
            grid->accept(initVisitor);
        }
        //else
        //{
        //    //rcp->restart((int)restartStep);
        //    rcp->readBlocks((int)restartStep);
        //    grid->accept(metisVisitor);
        //    rcp->readDataSet((int)restartStep);
        //    grid->setTimeStep(restartStep);
        //}


       

//}
//else
//{
//            //double restartStep = 10;
//            rcp->restart((int)restartStep);
//            grid->setTimeStep(restartStep);
//
//            //GbCylinder3DPtr geoInflow(new GbCylinder3D(-1.30181 + 0.0005, 0.390872 - 0.00229, g_maxX3 - 2.0 * dx, -1.30181 + 0.0005, 0.390872 - 0.00229, g_maxX3 + 2.0 * dx, 0.013));
//            //if (myid == 0) GbSystem3D::writeGeoObject(geoInflow.get(), outputPath + "/geo/geoInflow", WbWriterVtkXmlBinary::getInstance());
//            //SPtr<D3Q27Interactor> intrInflow = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, inflowConcreteBC, Interactor3D::SOLID));
//            //SetBcBlocksBlockVisitor v1(intrInflow);
//            //grid->accept(v1);
//            //intrInflow->initInteractor();
//
//            if (myid == 0)  UBLOG(logINFO, "Restart - end");
//}
        

        grid->accept(bcVisitor);
        //OneDistributionSetConnectorsBlockVisitor setConnsVisitor(comm);
        TwoDistributionsDoubleGhostLayerSetConnectorsBlockVisitor setConnsVisitor(comm);
        // ThreeDistributionsDoubleGhostLayerSetConnectorsBlockVisitor setConnsVisitor(comm);
        grid->accept(setConnsVisitor);

        int numOfThreads = 18;
        omp_set_num_threads(numOfThreads);

        SPtr<UbScheduler> nupsSch = std::make_shared<UbScheduler>(10, 10, 100);
        SPtr<NUPSCounterSimulationObserver> nupsSimulationObserver = make_shared<NUPSCounterSimulationObserver>(grid, nupsSch, numOfThreads, comm);

        //// write data for visualization of macroscopic quantities
        SPtr<UbScheduler> visSch(new UbScheduler(vtkSteps));
        // SPtr<UbScheduler> visSch(new UbScheduler(1, 8700, 8800));
        // visSch->addSchedule(1, 8700, 8800);
        SPtr<WriteSharpInterfaceQuantitiesSimulationObserver> writeMQSimulationObserver(new WriteSharpInterfaceQuantitiesSimulationObserver(grid, visSch, outputPath, WbWriterVtkXmlBinary::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));
        //writeMQSimulationObserver->update(10);

        //SPtr<WriteMacroscopicQuantitiesSimulationObserver> writeMQSimulationObserver(new WriteMacroscopicQuantitiesSimulationObserver(grid, visSch, outputPath, WbWriterVtkXmlBinary::getInstance(), SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));
        //writeMQSimulationObserver->update(0);

        //int endTime = 20;
        SPtr<Simulation> simulation(new Simulation(grid, lScheduler, endTime));
        simulation->addSimulationObserver(nupsSimulationObserver);
        //!!!//simulation->addSimulationObserver(lcSimulationObserver);
        simulation->addSimulationObserver(writeMQSimulationObserver);
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
    return 0;
}
