#include "CompressibleOffsetMomentsInterpolationProcessor.h"
#include "D3Q27System.h"

#include <lbm/Scaling.h>
#include <lbm/Interpolation_FC.h>
#include <lbm/Interpolation_CF.h>

//using namespace UbMath;
using namespace vf::basics::constant;

//////////////////////////////////////////////////////////////////////////
CompressibleOffsetMomentsInterpolationProcessor::CompressibleOffsetMomentsInterpolationProcessor(real omegaC, real omegaF)
   : omegaC(omegaC), omegaF(omegaF)
{
}
//////////////////////////////////////////////////////////////////////////
InterpolationProcessorPtr CompressibleOffsetMomentsInterpolationProcessor::clone()
{
   InterpolationProcessorPtr iproc = InterpolationProcessorPtr (new CompressibleOffsetMomentsInterpolationProcessor(this->omegaC, this->omegaF));

   dynamicPointerCast<CompressibleOffsetMomentsInterpolationProcessor>(iproc)->OxxPyyPzzC = this->OxxPyyPzzC;
   dynamicPointerCast<CompressibleOffsetMomentsInterpolationProcessor>(iproc)->OxxPyyPzzF = this->OxxPyyPzzF;

   return iproc;
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetMomentsInterpolationProcessor::setOmegas( real omegaC, real omegaF )
{
   this->omegaC = omegaC;
   this->omegaF = omegaF;

   real dtC = (c3o1 *shearViscosity)/((c1o1/omegaC)-c1o2);
   real dtF = (c3o1 *shearViscosity)/((c1o1/omegaF)-c1o2);

   if (bulkViscosity != 0)
   {
      this->OxxPyyPzzC = LBMSystem::calcOmega2(bulkViscosity, dtC);
      this->OxxPyyPzzF = LBMSystem::calcOmega2(bulkViscosity, dtF);
   }
   else
   {
      this->OxxPyyPzzC = c1o1;
      this->OxxPyyPzzF = c1o1;
   }
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetMomentsInterpolationProcessor::setOffsets(real xoff, real yoff, real zoff)
{
   this->xoff = xoff;
   this->yoff = yoff;
   this->zoff = zoff;     
   this->xoff_sq = xoff * xoff;
   this->yoff_sq = yoff * yoff;
   this->zoff_sq = zoff * zoff;
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetMomentsInterpolationProcessor::interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF, real xoff, real yoff, real zoff)
{
    setOffsets(xoff, yoff, zoff);

    vf::lbm::MomentsOnSourceNode moments_PPP;
    vf::lbm::MomentsOnSourceNode moments_MPP;
    vf::lbm::MomentsOnSourceNode moments_PMP;
    vf::lbm::MomentsOnSourceNode moments_MMP;
    vf::lbm::MomentsOnSourceNode moments_PPM;
    vf::lbm::MomentsOnSourceNode moments_MPM;
    vf::lbm::MomentsOnSourceNode moments_PMM;
    vf::lbm::MomentsOnSourceNode moments_MMM;

    vf::lbm::calculateMomentsOnSourceNodes(icellC.BSW, omegaC, moments_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(icellC.TSW, omegaC, moments_MMP);
    vf::lbm::calculateMomentsOnSourceNodes(icellC.TNW, omegaC, moments_MPP);
    vf::lbm::calculateMomentsOnSourceNodes(icellC.BNW, omegaC, moments_MPM);
    vf::lbm::calculateMomentsOnSourceNodes(icellC.BSE, omegaC, moments_PMM);
    vf::lbm::calculateMomentsOnSourceNodes(icellC.TNE, omegaC, moments_PPP);
    vf::lbm::calculateMomentsOnSourceNodes(icellC.TSE, omegaC, moments_PMP);
    vf::lbm::calculateMomentsOnSourceNodes(icellC.BNE, omegaC, moments_PPM);

    vf::lbm::Coefficients coefficients;
    calculateCoefficients(xoff, yoff, zoff, coefficients, moments_PPP, moments_MPP, moments_PMP, moments_MMP, moments_PPM, moments_MPM, moments_PMM, moments_MMM);

    interpolate_CF(icellF.BSW, omegaF, coefficients, -0.25, -0.25, -0.25);
    interpolate_CF(icellF.BNE, omegaF, coefficients,  0.25,  0.25, -0.25);
    interpolate_CF(icellF.TNW, omegaF, coefficients, -0.25,  0.25,  0.25);
    interpolate_CF(icellF.TSE, omegaF, coefficients,  0.25, -0.25,  0.25);
    interpolate_CF(icellF.BNW, omegaF, coefficients, -0.25,  0.25, -0.25);
    interpolate_CF(icellF.BSE, omegaF, coefficients,  0.25, -0.25, -0.25);
    interpolate_CF(icellF.TSW, omegaF, coefficients, -0.25, -0.25,  0.25);
    interpolate_CF(icellF.TNE, omegaF, coefficients,  0.25,  0.25,  0.25);
}

void CompressibleOffsetMomentsInterpolationProcessor::interpolate_CF(real* fine, real omegaF, vf::lbm::Coefficients &coefficients, real x, real y, real z) 
{
    real eps_new = c1o2; // ratio of grid resolutions
    real f[27];

    vf::lbm::interpolate_cf(x, y, z, f, coefficients, eps_new, omegaF);

    fine[DIR_P00] = f[DIR_P00];
    fine[DIR_M00] = f[DIR_M00];
    fine[DIR_0P0] = f[DIR_0P0];
    fine[DIR_0M0] = f[DIR_0M0];
    fine[DIR_00P] = f[DIR_00P];
    fine[DIR_00M] = f[DIR_00M];
    fine[DIR_PP0] = f[DIR_PP0];
    fine[DIR_MM0] = f[DIR_MM0];
    fine[DIR_PM0] = f[DIR_PM0];
    fine[DIR_MP0] = f[DIR_MP0];
    fine[DIR_P0P] = f[DIR_P0P];
    fine[DIR_M0M] = f[DIR_M0M];
    fine[DIR_P0M] = f[DIR_P0M];
    fine[DIR_M0P] = f[DIR_M0P];
    fine[DIR_0PP] = f[DIR_0PP];
    fine[DIR_0MM] = f[DIR_0MM];
    fine[DIR_0PM] = f[DIR_0PM];
    fine[DIR_0MP] = f[DIR_0MP];
    fine[DIR_000] = f[DIR_000];
    fine[DIR_PPP] = f[DIR_PPP];
    fine[DIR_PMP] = f[DIR_PMP];
    fine[DIR_PPM] = f[DIR_PPM];
    fine[DIR_PMM] = f[DIR_PMM];
    fine[DIR_MPP] = f[DIR_MPP];
    fine[DIR_MMP] = f[DIR_MMP];
    fine[DIR_MPM] = f[DIR_MPM];
    fine[DIR_MMM] = f[DIR_MMM];
}


//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetMomentsInterpolationProcessor::interpolateFineToCoarse(D3Q27ICell& icellF, real* icellC, real xoff, real yoff, real zoff)
{
   setOffsets(xoff, yoff, zoff);
   calcInterpolation(icellF, icellC);
}

void CompressibleOffsetMomentsInterpolationProcessor::calcInterpolation(const D3Q27ICell& icell, real* icellC)
{
    const real eps_new = c2o1; // ratio of grid resolutions

    vf::lbm::MomentsOnSourceNode moments_PPP;
    vf::lbm::MomentsOnSourceNode moments_MPP;
    vf::lbm::MomentsOnSourceNode moments_PMP;
    vf::lbm::MomentsOnSourceNode moments_MMP;
    vf::lbm::MomentsOnSourceNode moments_PPM;
    vf::lbm::MomentsOnSourceNode moments_MPM;
    vf::lbm::MomentsOnSourceNode moments_PMM;
    vf::lbm::MomentsOnSourceNode moments_MMM;

    vf::lbm::calculateMomentsOnSourceNodes(icell.BSW, omegaF, moments_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(icell.TSW, omegaF, moments_MMP);
    vf::lbm::calculateMomentsOnSourceNodes(icell.TNW, omegaF, moments_MPP);
    vf::lbm::calculateMomentsOnSourceNodes(icell.BNW, omegaF, moments_MPM);
    vf::lbm::calculateMomentsOnSourceNodes(icell.BSE, omegaF, moments_PMM);
    vf::lbm::calculateMomentsOnSourceNodes(icell.TNE, omegaF, moments_PPP);
    vf::lbm::calculateMomentsOnSourceNodes(icell.TSE, omegaF, moments_PMP);
    vf::lbm::calculateMomentsOnSourceNodes(icell.BNE, omegaF, moments_PPM);

    vf::lbm::Coefficients coefficients;
    calculateCoefficients(xoff, yoff, zoff, coefficients, moments_PPP, moments_MPP, moments_PMP, moments_MMP, moments_PPM, moments_MPM, moments_PMM, moments_MMM);

    vf::lbm::interpolate_fc(icellC, eps_new, omegaC, coefficients);
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetMomentsInterpolationProcessor::calcMoments(const real* const f, real omega, real& press, real& vx1, real& vx2, real& vx3, 
                                                    real& kxy, real& kyz, real& kxz, real& kxxMyy, real& kxxMzz)
{
   using namespace D3Q27System;
   using namespace vf::lbm::dir;

   real drho = 0.0;
   D3Q27System::calcCompMacroscopicValues(f,drho,vx1,vx2,vx3);
   
   press = drho; //interpolate rho!

   kxy   = -3.*omega*((((f[DIR_MMP]+f[DIR_PPM])-(f[DIR_MPP]+f[DIR_PMM]))+((f[DIR_MMM]+f[DIR_PPP])-(f[DIR_MPM]+f[DIR_PMP])))+((f[DIR_MM0]+f[DIR_PP0])-(f[DIR_MP0]+f[DIR_PM0]))/(c1o1 + drho)-(vx1*vx2));// might not be optimal MG 25.2.13
   kyz   = -3.*omega*((((f[DIR_MMM]+f[DIR_PPP])-(f[DIR_PMP]+f[DIR_MPM]))+((f[DIR_PMM]+f[DIR_MPP])-(f[DIR_MMP]+f[DIR_PPM])))+((f[DIR_0MM]+f[DIR_0PP])-(f[DIR_0MP]+f[DIR_0PM]))/(c1o1 + drho)-(vx2*vx3));
   kxz   = -3.*omega*((((f[DIR_MPM]+f[DIR_PMP])-(f[DIR_MMP]+f[DIR_PPM]))+((f[DIR_MMM]+f[DIR_PPP])-(f[DIR_PMM]+f[DIR_MPP])))+((f[DIR_M0M]+f[DIR_P0P])-(f[DIR_M0P]+f[DIR_P0M]))/(c1o1 + drho)-(vx1*vx3));
   kxxMyy = -3./2.*omega*((((f[DIR_M0M]+f[DIR_P0P])-(f[DIR_0MM]+f[DIR_0PP]))+((f[DIR_M0P]+f[DIR_P0M])-(f[DIR_0MP]+f[DIR_0PM])))+((f[DIR_M00]+f[DIR_P00])-(f[DIR_0M0]+f[DIR_0P0]))/(c1o1 + drho)-(vx1*vx1-vx2*vx2));
   kxxMzz = -3./2.*omega*((((f[DIR_MP0]+f[DIR_PM0])-(f[DIR_0MM]+f[DIR_0PP]))+((f[DIR_MM0]+f[DIR_PP0])-(f[DIR_0MP]+f[DIR_0PM])))+((f[DIR_M00]+f[DIR_P00])-(f[DIR_00M]+f[DIR_00P]))/(c1o1 + drho)-(vx1*vx1-vx3*vx3));
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetMomentsInterpolationProcessor::calcInterpolatedCoefficiets(const D3Q27ICell& icell, real omega, real eps_new)
{
   real        vx1_SWT,vx2_SWT,vx3_SWT;
   real        vx1_NWT,vx2_NWT,vx3_NWT;
   real        vx1_NET,vx2_NET,vx3_NET;
   real        vx1_SET,vx2_SET,vx3_SET;
   real        vx1_SWB,vx2_SWB,vx3_SWB;
   real        vx1_NWB,vx2_NWB,vx3_NWB;
   real        vx1_NEB,vx2_NEB,vx3_NEB;
   real        vx1_SEB,vx2_SEB,vx3_SEB;

   real        kxyFromfcNEQ_SWT, kyzFromfcNEQ_SWT, kxzFromfcNEQ_SWT, kxxMyyFromfcNEQ_SWT, kxxMzzFromfcNEQ_SWT;
   real        kxyFromfcNEQ_NWT, kyzFromfcNEQ_NWT, kxzFromfcNEQ_NWT, kxxMyyFromfcNEQ_NWT, kxxMzzFromfcNEQ_NWT;
   real        kxyFromfcNEQ_NET, kyzFromfcNEQ_NET, kxzFromfcNEQ_NET, kxxMyyFromfcNEQ_NET, kxxMzzFromfcNEQ_NET;
   real        kxyFromfcNEQ_SET, kyzFromfcNEQ_SET, kxzFromfcNEQ_SET, kxxMyyFromfcNEQ_SET, kxxMzzFromfcNEQ_SET;
   real        kxyFromfcNEQ_SWB, kyzFromfcNEQ_SWB, kxzFromfcNEQ_SWB, kxxMyyFromfcNEQ_SWB, kxxMzzFromfcNEQ_SWB;
   real        kxyFromfcNEQ_NWB, kyzFromfcNEQ_NWB, kxzFromfcNEQ_NWB, kxxMyyFromfcNEQ_NWB, kxxMzzFromfcNEQ_NWB;
   real        kxyFromfcNEQ_NEB, kyzFromfcNEQ_NEB, kxzFromfcNEQ_NEB, kxxMyyFromfcNEQ_NEB, kxxMzzFromfcNEQ_NEB;
   real        kxyFromfcNEQ_SEB, kyzFromfcNEQ_SEB, kxzFromfcNEQ_SEB, kxxMyyFromfcNEQ_SEB, kxxMzzFromfcNEQ_SEB;

//     real drho_PPP, vx1_PPP, vx2_PPP, vx3_PPP;
//     real drho_MPP, vx1_MPP, vx2_MPP, vx3_MPP;
//     real drho_PMP, vx1_PMP, vx2_PMP, vx3_PMP;
//     real drho_MMP, vx1_MMP, vx2_MMP, vx3_MMP;
//     real drho_PPM, vx1_PPM, vx2_PPM, vx3_PPM;
//     real drho_MPM, vx1_MPM, vx2_MPM, vx3_MPM;
//     real drho_PMM, vx1_PMM, vx2_PMM, vx3_PMM;
//     real drho_MMM, vx1_MMM, vx2_MMM, vx3_MMM;

//     real kxyFromfcNEQ_PPP, kyzFromfcNEQ_PPP, kxzFromfcNEQ_PPP, kxxMyyFromfcNEQ_PPP, kxxMzzFromfcNEQ_PPP;
//     real kxyFromfcNEQ_MPP, kyzFromfcNEQ_MPP, kxzFromfcNEQ_MPP, kxxMyyFromfcNEQ_MPP, kxxMzzFromfcNEQ_MPP;
//     real kxyFromfcNEQ_PMP, kyzFromfcNEQ_PMP, kxzFromfcNEQ_PMP, kxxMyyFromfcNEQ_PMP, kxxMzzFromfcNEQ_PMP;
//     real kxyFromfcNEQ_MMP, kyzFromfcNEQ_MMP, kxzFromfcNEQ_MMP, kxxMyyFromfcNEQ_MMP, kxxMzzFromfcNEQ_MMP;
//     real kxyFromfcNEQ_PPM, kyzFromfcNEQ_PPM, kxzFromfcNEQ_PPM, kxxMyyFromfcNEQ_PPM, kxxMzzFromfcNEQ_PPM;
//     real kxyFromfcNEQ_MPM, kyzFromfcNEQ_MPM, kxzFromfcNEQ_MPM, kxxMyyFromfcNEQ_MPM, kxxMzzFromfcNEQ_MPM;
//     real kxyFromfcNEQ_PMM, kyzFromfcNEQ_PMM, kxzFromfcNEQ_PMM, kxxMyyFromfcNEQ_PMM, kxxMzzFromfcNEQ_PMM;
//     real kxyFromfcNEQ_MMM, kyzFromfcNEQ_MMM, kxzFromfcNEQ_MMM, kxxMyyFromfcNEQ_MMM, kxxMzzFromfcNEQ_MMM;


//    vf::lbm::calculateMomentsOnSourceNodes(icell.TSW,omega,press_SWT,vx1_SWT,vx2_SWT,vx3_SWT, kxyFromfcNEQ_SWT, kyzFromfcNEQ_SWT, kxzFromfcNEQ_SWT, kxxMyyFromfcNEQ_SWT, kxxMzzFromfcNEQ_SWT);
//    vf::lbm::calculateMomentsOnSourceNodes(icell.TNW,omega,press_NWT,vx1_NWT,vx2_NWT,vx3_NWT, kxyFromfcNEQ_NWT, kyzFromfcNEQ_NWT, kxzFromfcNEQ_NWT, kxxMyyFromfcNEQ_NWT, kxxMzzFromfcNEQ_NWT);
//    vf::lbm::calculateMomentsOnSourceNodes(icell.TNE,omega,press_NET,vx1_NET,vx2_NET,vx3_NET, kxyFromfcNEQ_NET, kyzFromfcNEQ_NET, kxzFromfcNEQ_NET, kxxMyyFromfcNEQ_NET, kxxMzzFromfcNEQ_NET);
//    vf::lbm::calculateMomentsOnSourceNodes(icell.TSE,omega,press_SET,vx1_SET,vx2_SET,vx3_SET, kxyFromfcNEQ_SET, kyzFromfcNEQ_SET, kxzFromfcNEQ_SET, kxxMyyFromfcNEQ_SET, kxxMzzFromfcNEQ_SET);
//    vf::lbm::calculateMomentsOnSourceNodes(icell.BSW,omega,press_SWB,vx1_SWB,vx2_SWB,vx3_SWB, kxyFromfcNEQ_SWB, kyzFromfcNEQ_SWB, kxzFromfcNEQ_SWB, kxxMyyFromfcNEQ_SWB, kxxMzzFromfcNEQ_SWB);
//    vf::lbm::calculateMomentsOnSourceNodes(icell.BNW,omega,press_NWB,vx1_NWB,vx2_NWB,vx3_NWB, kxyFromfcNEQ_NWB, kyzFromfcNEQ_NWB, kxzFromfcNEQ_NWB, kxxMyyFromfcNEQ_NWB, kxxMzzFromfcNEQ_NWB);
//    vf::lbm::calculateMomentsOnSourceNodes(icell.BNE,omega,press_NEB,vx1_NEB,vx2_NEB,vx3_NEB, kxyFromfcNEQ_NEB, kyzFromfcNEQ_NEB, kxzFromfcNEQ_NEB, kxxMyyFromfcNEQ_NEB, kxxMzzFromfcNEQ_NEB);
//    vf::lbm::calculateMomentsOnSourceNodes(icell.BSE,omega,press_SEB,vx1_SEB,vx2_SEB,vx3_SEB, kxyFromfcNEQ_SEB, kyzFromfcNEQ_SEB, kxzFromfcNEQ_SEB, kxxMyyFromfcNEQ_SEB, kxxMzzFromfcNEQ_SEB);



//     vf::lbm::calculateMomentsOnSourceNodes(icell.BSW, omegaF, drho_MMM, vx1_MMM, vx2_MMM, vx3_MMM,kxyFromfcNEQ_MMM, kyzFromfcNEQ_MMM, kxzFromfcNEQ_MMM, kxxMyyFromfcNEQ_MMM, kxxMzzFromfcNEQ_MMM);
//     vf::lbm::calculateMomentsOnSourceNodes(icell.TSW, omegaF, drho_MMP, vx1_MMP, vx2_MMP, vx3_MMP,kxyFromfcNEQ_MMP, kyzFromfcNEQ_MMP, kxzFromfcNEQ_MMP, kxxMyyFromfcNEQ_MMP, kxxMzzFromfcNEQ_MMP);
//     vf::lbm::calculateMomentsOnSourceNodes(icell.TNW, omegaF, drho_MPP, vx1_MPP, vx2_MPP, vx3_MPP,kxyFromfcNEQ_MPP, kyzFromfcNEQ_MPP, kxzFromfcNEQ_MPP, kxxMyyFromfcNEQ_MPP, kxxMzzFromfcNEQ_MPP);
//     vf::lbm::calculateMomentsOnSourceNodes(icell.BNW, omegaF, drho_MPM, vx1_MPM, vx2_MPM, vx3_MPM,kxyFromfcNEQ_MPM, kyzFromfcNEQ_MPM, kxzFromfcNEQ_MPM, kxxMyyFromfcNEQ_MPM, kxxMzzFromfcNEQ_MPM);
//     vf::lbm::calculateMomentsOnSourceNodes(icell.BSE, omegaF, drho_PMM, vx1_PMM, vx2_PMM, vx3_PMM,kxyFromfcNEQ_PMM, kyzFromfcNEQ_PMM, kxzFromfcNEQ_PMM, kxxMyyFromfcNEQ_PMM, kxxMzzFromfcNEQ_PMM);
//     vf::lbm::calculateMomentsOnSourceNodes(icell.TNE, omegaF, drho_PPP, vx1_PPP, vx2_PPP, vx3_PPP,kxyFromfcNEQ_PPP, kyzFromfcNEQ_PPP, kxzFromfcNEQ_PPP, kxxMyyFromfcNEQ_PPP, kxxMzzFromfcNEQ_PPP);
//     vf::lbm::calculateMomentsOnSourceNodes(icell.TSE, omegaF, drho_PMP, vx1_PMP, vx2_PMP, vx3_PMP,kxyFromfcNEQ_PMP, kyzFromfcNEQ_PMP, kxzFromfcNEQ_PMP, kxxMyyFromfcNEQ_PMP, kxxMzzFromfcNEQ_PMP);
//     vf::lbm::calculateMomentsOnSourceNodes(icell.BNE, omegaF, drho_PPM, vx1_PPM, vx2_PPM, vx3_PPM,kxyFromfcNEQ_PPM, kyzFromfcNEQ_PPM, kxzFromfcNEQ_PPM, kxxMyyFromfcNEQ_PPM, kxxMzzFromfcNEQ_PPM);



//     a_000 = c1o64 * (
//             c2o1 * (
//             ((kxyFromfcNEQ_MMM - kxyFromfcNEQ_PPP) + (kxyFromfcNEQ_MMP - kxyFromfcNEQ_PPM)) + ((kxyFromfcNEQ_PMM - kxyFromfcNEQ_MPP) + (kxyFromfcNEQ_PMP - kxyFromfcNEQ_MPM)) + 
//             ((kxzFromfcNEQ_MMM - kxzFromfcNEQ_PPP) + (kxzFromfcNEQ_PPM - kxzFromfcNEQ_MMP)) + ((kxzFromfcNEQ_PMM - kxzFromfcNEQ_MPP) + (kxzFromfcNEQ_MPM - kxzFromfcNEQ_PMP)) + 
//             ((vx2_PPP + vx2_MMM) + (vx2_PPM + vx2_MMP)) - ((vx2_MPP + vx2_PMM) + (vx2_MPM + vx2_PMP)) + 
//             ((vx3_PPP + vx3_MMM) - (vx3_PPM + vx3_MMP)) + ((vx3_PMP + vx3_MPM) - (vx3_MPP + vx3_PMM))) + 
//             c8o1 * (((vx1_PPP + vx1_MMM) + (vx1_PPM + vx1_MMP)) + ((vx1_MPP + vx1_PMM) + (vx1_PMP + vx1_MPM))) +
//             ((kxxMyyFromfcNEQ_MMM - kxxMyyFromfcNEQ_PPP) + (kxxMyyFromfcNEQ_MMP - kxxMyyFromfcNEQ_PPM)) + 
//             ((kxxMyyFromfcNEQ_MPP - kxxMyyFromfcNEQ_PMM) + (kxxMyyFromfcNEQ_MPM - kxxMyyFromfcNEQ_PMP)) +
//             ((kxxMzzFromfcNEQ_MMM - kxxMzzFromfcNEQ_PPP) + (kxxMzzFromfcNEQ_MMP - kxxMzzFromfcNEQ_PPM)) + 
//             ((kxxMzzFromfcNEQ_MPP - kxxMzzFromfcNEQ_PMM) + (kxxMzzFromfcNEQ_MPM - kxxMzzFromfcNEQ_PMP)));
//     b_000 = c1o64 * (
//             c2o1 * (
//             ((kxxMyyFromfcNEQ_PPP - kxxMyyFromfcNEQ_MMM) + (kxxMyyFromfcNEQ_PPM - kxxMyyFromfcNEQ_MMP)) + 
//             ((kxxMyyFromfcNEQ_MPP - kxxMyyFromfcNEQ_PMM) + (kxxMyyFromfcNEQ_MPM - kxxMyyFromfcNEQ_PMP)) + 
//             ((kxyFromfcNEQ_MMM - kxyFromfcNEQ_PPP) + (kxyFromfcNEQ_MMP - kxyFromfcNEQ_PPM)) + 
//             ((kxyFromfcNEQ_MPP - kxyFromfcNEQ_PMM) + (kxyFromfcNEQ_MPM - kxyFromfcNEQ_PMP)) + 
//             ((kyzFromfcNEQ_MMM - kyzFromfcNEQ_PPP) + (kyzFromfcNEQ_PPM - kyzFromfcNEQ_MMP)) + 
//             ((kyzFromfcNEQ_PMM - kyzFromfcNEQ_MPP) + (kyzFromfcNEQ_MPM - kyzFromfcNEQ_PMP)) + 
//             ((vx1_PPP + vx1_MMM) + (vx1_PPM + vx1_MMP)) - ((vx1_MPM + vx1_MPP) + (vx1_PMM + vx1_PMP)) + 
//             ((vx3_PPP + vx3_MMM) - (vx3_PPM + vx3_MMP)) + ((vx3_MPP + vx3_PMM) - (vx3_MPM + vx3_PMP))) + 
//             c8o1 * (((vx2_PPP + vx2_MMM) + (vx2_PPM + vx2_MMP)) + ((vx2_MPP + vx2_PMM) + (vx2_MPM + vx2_PMP))) + 
//             ((kxxMzzFromfcNEQ_MMM - kxxMzzFromfcNEQ_PPP) + (kxxMzzFromfcNEQ_MMP - kxxMzzFromfcNEQ_PPM)) +
//             ((kxxMzzFromfcNEQ_PMM - kxxMzzFromfcNEQ_MPP) + (kxxMzzFromfcNEQ_PMP - kxxMzzFromfcNEQ_MPM)));
//     c_000 = c1o64 * ( 
//             c2o1 * (
//             ((kxxMzzFromfcNEQ_PPP - kxxMzzFromfcNEQ_MMM) + (kxxMzzFromfcNEQ_MMP - kxxMzzFromfcNEQ_PPM)) + 
//             ((kxxMzzFromfcNEQ_MPP - kxxMzzFromfcNEQ_PMM) + (kxxMzzFromfcNEQ_PMP - kxxMzzFromfcNEQ_MPM)) + 
//             ((kxzFromfcNEQ_MMM - kxzFromfcNEQ_PPP) + (kxzFromfcNEQ_MMP - kxzFromfcNEQ_PPM)) + 
//             ((kxzFromfcNEQ_MPP - kxzFromfcNEQ_PMM) + (kxzFromfcNEQ_MPM - kxzFromfcNEQ_PMP)) + 
//             ((kyzFromfcNEQ_MMM - kyzFromfcNEQ_PPP) + (kyzFromfcNEQ_MMP - kyzFromfcNEQ_PPM)) + 
//             ((kyzFromfcNEQ_PMM - kyzFromfcNEQ_MPP) + (kyzFromfcNEQ_PMP - kyzFromfcNEQ_MPM)) + 
//             ((vx1_PPP + vx1_MMM) - (vx1_MMP + vx1_PPM)) + ((vx1_MPM + vx1_PMP) - (vx1_MPP + vx1_PMM)) + 
//             ((vx2_PPP + vx2_MMM) - (vx2_MMP + vx2_PPM)) + ((vx2_MPP + vx2_PMM) - (vx2_MPM + vx2_PMP))) + 
//             c8o1 * (((vx3_PPP + vx3_MMM) + (vx3_PPM + vx3_MMP)) + ((vx3_PMM + vx3_MPP) + (vx3_PMP + vx3_MPM))) +
//             ((kxxMyyFromfcNEQ_MMM - kxxMyyFromfcNEQ_PPP) + (kxxMyyFromfcNEQ_PPM - kxxMyyFromfcNEQ_MMP)) + 
//             ((kxxMyyFromfcNEQ_PMM - kxxMyyFromfcNEQ_MPP) + (kxxMyyFromfcNEQ_MPM - kxxMyyFromfcNEQ_PMP)));

//     a_100 = c1o4 * (((vx1_PPP - vx1_MMM) + (vx1_PPM - vx1_MMP)) + ((vx1_PMM - vx1_MPP) + (vx1_PMP - vx1_MPM)));
//     b_100 = c1o4 * (((vx2_PPP - vx2_MMM) + (vx2_PPM - vx2_MMP)) + ((vx2_PMM - vx2_MPP) + (vx2_PMP - vx2_MPM)));
//     c_100 = c1o4 * (((vx3_PPP - vx3_MMM) + (vx3_PPM - vx3_MMP)) + ((vx3_PMM - vx3_MPP) + (vx3_PMP - vx3_MPM)));

//     a_200 = c1o16 * ( 
//             c2o1 * (
//             ((vx2_PPP + vx2_MMM) + (vx2_PPM - vx2_MPP)) + ((vx2_MMP - vx2_PMM) - (vx2_MPM + vx2_PMP)) + 
//             ((vx3_PPP + vx3_MMM) - (vx3_PPM + vx3_MPP)) + ((vx3_MPM + vx3_PMP) - (vx3_MMP + vx3_PMM))) + 
//             ((kxxMyyFromfcNEQ_PPP - kxxMyyFromfcNEQ_MMM) + (kxxMyyFromfcNEQ_PPM - kxxMyyFromfcNEQ_MMP)) + 
//             ((kxxMyyFromfcNEQ_PMM - kxxMyyFromfcNEQ_MPP) + (kxxMyyFromfcNEQ_PMP - kxxMyyFromfcNEQ_MPM)) + 
//             ((kxxMzzFromfcNEQ_PPP - kxxMzzFromfcNEQ_MMM) + (kxxMzzFromfcNEQ_PPM - kxxMzzFromfcNEQ_MMP)) + 
//             ((kxxMzzFromfcNEQ_PMM - kxxMzzFromfcNEQ_MPP) + (kxxMzzFromfcNEQ_PMP - kxxMzzFromfcNEQ_MPM)));
//     b_200 = c1o8 * (
//             c2o1 * (
//             -((vx1_PPP + vx1_MMM) + (vx1_PPM + vx1_MMP)) + ((vx1_MPP + vx1_PMM) + (vx1_MPM + vx1_PMP))) +
//             ((kxyFromfcNEQ_PPP - kxyFromfcNEQ_MMM) + (kxyFromfcNEQ_PPM - kxyFromfcNEQ_MMP)) + 
//             ((kxyFromfcNEQ_PMM - kxyFromfcNEQ_MPP) + (kxyFromfcNEQ_PMP - kxyFromfcNEQ_MPM)));
//     c_200 = c1o8 * (
//             c2o1 * (
//             ((vx1_PPM + vx1_MMP) - (vx1_PPP + vx1_MMM)) + ((vx1_MPP + vx1_PMM) - (vx1_MPM + vx1_PMP))) +
//             ((kxzFromfcNEQ_PPP - kxzFromfcNEQ_MMM) + (kxzFromfcNEQ_PPM - kxzFromfcNEQ_MMP)) + 
//             ((kxzFromfcNEQ_PMM - kxzFromfcNEQ_MPP) + (kxzFromfcNEQ_PMP - kxzFromfcNEQ_MPM)));

//     a_010 = c1o4 * (((vx1_PPP - vx1_MMM) + (vx1_PPM - vx1_MMP)) + ((vx1_MPP - vx1_PMM) + (vx1_MPM - vx1_PMP)));
//     b_010 = c1o4 * (((vx2_PPP - vx2_MMM) + (vx2_PPM - vx2_MMP)) + ((vx2_MPP - vx2_PMM) + (vx2_MPM - vx2_PMP)));
//     c_010 = c1o4 * (((vx3_PPP - vx3_MMM) + (vx3_PPM - vx3_MMP)) + ((vx3_MPP - vx3_PMM) + (vx3_MPM - vx3_PMP)));

//     a_020 = c1o8 * (
//             c2o1 * (-((vx2_PPP + vx2_MMM) + (vx2_MMP + vx2_PPM)) + ((vx2_MPP + vx2_PMM) + (vx2_MPM + vx2_PMP))) +
//             ((kxyFromfcNEQ_PPP - kxyFromfcNEQ_MMM) + (kxyFromfcNEQ_PPM - kxyFromfcNEQ_MMP)) + 
//             ((kxyFromfcNEQ_MPP - kxyFromfcNEQ_PMM) + (kxyFromfcNEQ_MPM - kxyFromfcNEQ_PMP)));
//     b_020 = c1o16 * (
//             c2o1 * (
//             ((kxxMyyFromfcNEQ_MMM - kxxMyyFromfcNEQ_PPP) + (kxxMyyFromfcNEQ_MMP - kxxMyyFromfcNEQ_PPM)) +
//             ((kxxMyyFromfcNEQ_PMM - kxxMyyFromfcNEQ_MPP) + (kxxMyyFromfcNEQ_PMP - kxxMyyFromfcNEQ_MPM)) +
//             ((vx1_PPP + vx1_MMM) + (vx1_PPM + vx1_MMP)) - ((vx1_MPP + vx1_PMM) + (vx1_PMP + vx1_MPM)) + 
//             ((vx3_PPP + vx3_MMM) - (vx3_PPM + vx3_MMP)) + ((vx3_MPP + vx3_PMM) - (vx3_MPM + vx3_PMP))) +
//             ((kxxMzzFromfcNEQ_PPP - kxxMzzFromfcNEQ_MMM) + (kxxMzzFromfcNEQ_PPM - kxxMzzFromfcNEQ_MMP)) + 
//             ((kxxMzzFromfcNEQ_MPP - kxxMzzFromfcNEQ_PMM) + (kxxMzzFromfcNEQ_MPM - kxxMzzFromfcNEQ_PMP)));
//     c_020 = c1o8 * (
//             c2o1 * (((vx2_MMP + vx2_PPM) - (vx2_PPP + vx2_MMM)) + ((vx2_PMP + vx2_MPM) - (vx2_MPP + vx2_PMM))) +
//             ((kyzFromfcNEQ_PPP - kyzFromfcNEQ_MMM) + (kyzFromfcNEQ_PPM - kyzFromfcNEQ_MMP)) +
//             ((kyzFromfcNEQ_MPP - kyzFromfcNEQ_PMM) + (kyzFromfcNEQ_MPM - kyzFromfcNEQ_PMP)));

//     a_001 = c1o4 * (((vx1_PPP - vx1_MMM) + (vx1_MMP - vx1_PPM)) + ((vx1_MPP - vx1_PMM) + (vx1_PMP - vx1_MPM)));
//     b_001 = c1o4 * (((vx2_PPP - vx2_MMM) + (vx2_MMP - vx2_PPM)) + ((vx2_MPP - vx2_PMM) + (vx2_PMP - vx2_MPM)));
//     c_001 = c1o4 * (((vx3_PPP - vx3_MMM) + (vx3_MMP - vx3_PPM)) + ((vx3_MPP - vx3_PMM) + (vx3_PMP - vx3_MPM)));

//     a_002 = c1o8 * (
//             c2o1 * (((vx3_PPM + vx3_MMP) - (vx3_PPP + vx3_MMM)) + ((vx3_MPP + vx3_PMM) - (vx3_PMP + vx3_MPM))) +
//                     ((kxzFromfcNEQ_PPP - kxzFromfcNEQ_MMM) + (kxzFromfcNEQ_MMP - kxzFromfcNEQ_PPM)) +
//                     ((kxzFromfcNEQ_PMP - kxzFromfcNEQ_MPM) + (kxzFromfcNEQ_MPP - kxzFromfcNEQ_PMM)));
//     b_002 = c1o8 * (
//             c2o1 * (((vx3_PPM + vx3_MMP) - (vx3_PPP + vx3_MMM)) + ((vx3_MPM + vx3_PMP) - (vx3_PMM + vx3_MPP))) + 
//                     ((kyzFromfcNEQ_PPP - kyzFromfcNEQ_MMM) + (kyzFromfcNEQ_MMP - kyzFromfcNEQ_PPM)) + 
//                     ((kyzFromfcNEQ_PMP - kyzFromfcNEQ_MPM) + (kyzFromfcNEQ_MPP - kyzFromfcNEQ_PMM)));
//     c_002 = c1o16 * (
//             c2o1 * (
//             ((kxxMzzFromfcNEQ_MMM - kxxMzzFromfcNEQ_PPP) + (kxxMzzFromfcNEQ_PPM - kxxMzzFromfcNEQ_MMP)) + 
//             ((kxxMzzFromfcNEQ_MPM - kxxMzzFromfcNEQ_PMP) + (kxxMzzFromfcNEQ_PMM - kxxMzzFromfcNEQ_MPP)) + 
//             ((vx1_PPP + vx1_MMM) - (vx1_MMP + vx1_PPM)) + ((vx1_MPM + vx1_PMP) - (vx1_PMM + vx1_MPP)) + 
//             ((vx2_PPP + vx2_MMM) - (vx2_MMP + vx2_PPM)) + ((vx2_PMM + vx2_MPP) - (vx2_MPM + vx2_PMP))) + 
//             ((kxxMyyFromfcNEQ_PPP - kxxMyyFromfcNEQ_MMM) + (kxxMyyFromfcNEQ_MMP - kxxMyyFromfcNEQ_PPM)) +
//             ((kxxMyyFromfcNEQ_PMP - kxxMyyFromfcNEQ_MPM) + (kxxMyyFromfcNEQ_MPP - kxxMyyFromfcNEQ_PMM)));

//     a_110 = c1o2 * (((vx1_PPP + vx1_MMM) + (vx1_MMP + vx1_PPM)) - ((vx1_MPM + vx1_PMP) + (vx1_PMM + vx1_MPP)));
//     b_110 = c1o2 * (((vx2_PPP + vx2_MMM) + (vx2_MMP + vx2_PPM)) - ((vx2_MPM + vx2_PMP) + (vx2_PMM + vx2_MPP)));
//     c_110 = c1o2 * (((vx3_PPP + vx3_MMM) + (vx3_MMP + vx3_PPM)) - ((vx3_MPM + vx3_PMP) + (vx3_PMM + vx3_MPP)));

//     a_101 = c1o2 * (((vx1_PPP + vx1_MMM) - (vx1_MMP + vx1_PPM)) + ((vx1_MPM + vx1_PMP) - (vx1_PMM + vx1_MPP)));
//     b_101 = c1o2 * (((vx2_PPP + vx2_MMM) - (vx2_MMP + vx2_PPM)) + ((vx2_MPM + vx2_PMP) - (vx2_PMM + vx2_MPP)));
//     c_101 = c1o2 * (((vx3_PPP + vx3_MMM) - (vx3_MMP + vx3_PPM)) + ((vx3_MPM + vx3_PMP) - (vx3_PMM + vx3_MPP)));
    
//     a_011 = c1o2 * (((vx1_PPP + vx1_MMM) - (vx1_MMP + vx1_PPM)) + ((vx1_PMM + vx1_MPP) - (vx1_MPM + vx1_PMP)));
//     b_011 = c1o2 * (((vx2_PPP + vx2_MMM) - (vx2_MMP + vx2_PPM)) + ((vx2_PMM + vx2_MPP) - (vx2_MPM + vx2_PMP)));
//     c_011 = c1o2 * (((vx3_PPP + vx3_MMM) - (vx3_MMP + vx3_PPM)) + ((vx3_PMM + vx3_MPP) - (vx3_MPM + vx3_PMP)));

//     //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
//     kxyAverage    = c0o1;
//     kyzAverage    = c0o1;
//     kxzAverage    = c0o1;
//     kxxMyyAverage = c0o1;
//     kxxMzzAverage = c0o1;

//     ////////////////////////////////////////////////////////////////////////////////
//     //! - Set the relative position of the offset cell {-1, 0, 1}
//     //!
//     // real xoff    = offsetFC.xOffFC[nodeIndex];
//     // real yoff    = offsetFC.yOffFC[nodeIndex];
//     // real zoff    = offsetFC.zOffFC[nodeIndex];
     
//     real xoff_sq = xoff * xoff;
//     real yoff_sq = yoff * yoff;
//     real zoff_sq = zoff * zoff;

//     //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//     //!- Calculate coefficients for the polynomial interpolation of the pressure
//     //! 
//     real LaplaceRho = 
//         ((xoff != c0o1) || (yoff != c0o1) || (zoff != c0o1))
//         ? c0o1 : -c3o1 * (a_100 * a_100 + b_010 * b_010 + c_001 * c_001) - c6o1 * (b_100 * a_010 + c_100 * a_001 + c_010 * b_001);
//     d_000 =  c1o8 * ((((drho_PPP + drho_MMM) + (drho_PPM + drho_MMP)) + ((drho_PMM + drho_MPP) + (drho_PMP + drho_MPM))) - c2o1 * LaplaceRho);
//     d_100 = c1o4 * (((drho_PPP - drho_MMM) + (drho_PPM - drho_MMP)) + ((drho_PMM - drho_MPP) + (drho_PMP - drho_MPM)));
//     d_010 = c1o4 * (((drho_PPP - drho_MMM) + (drho_PPM - drho_MMP)) + ((drho_MPP - drho_PMM) + (drho_MPM - drho_PMP)));
//     d_001 = c1o4 * (((drho_PPP - drho_MMM) + (drho_MMP - drho_PPM)) + ((drho_MPP - drho_PMM) + (drho_PMP - drho_MPM)));
//     d_110 = c1o2 * (((drho_PPP + drho_MMM) + (drho_PPM + drho_MMP)) - ((drho_PMM + drho_MPP) + (drho_PMP + drho_MPM)));
//     d_101 = c1o2 * (((drho_PPP + drho_MMM) - (drho_PPM + drho_MMP)) + ((drho_PMP + drho_MPM) - (drho_PMM + drho_MPP)));
//     d_011 = c1o2 * (((drho_PPP + drho_MMM) - (drho_PPM + drho_MMP)) + ((drho_PMM + drho_MPP) - (drho_PMP + drho_MPM)));


//     //////////////////////////////////////////////////////////////////////////
//     //! - Extrapolation for refinement in to the wall (polynomial coefficients)
//     //!
//     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//     //
//     // x------x
//     // |      |
//     // |   ---+--->X
//     // |      |  |
//     // x------x  |
//     //          offset-vector
//     //
//     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//     a_000 = a_000 + xoff * a_100 + yoff * a_010 + zoff * a_001 + xoff_sq * a_200 + yoff_sq * a_020 + zoff_sq * a_002 +
//             xoff * yoff * a_110 + xoff * zoff * a_101 + yoff * zoff * a_011;
//     a_100 = a_100 + c2o1 * xoff * a_200 + yoff * a_110 + zoff * a_101;
//     a_010 = a_010 + c2o1 * yoff * a_020 + xoff * a_110 + zoff * a_011;
//     a_001 = a_001 + c2o1 * zoff * a_002 + xoff * a_101 + yoff * a_011;
//     b_000 = b_000 + xoff * b_100 + yoff * b_010 + zoff * b_001 + xoff_sq * b_200 + yoff_sq * b_020 + zoff_sq * b_002 +
//             xoff * yoff * b_110 + xoff * zoff * b_101 + yoff * zoff * b_011;
//     b_100 = b_100 + c2o1 * xoff * b_200 + yoff * b_110 + zoff * b_101;
//     b_010 = b_010 + c2o1 * yoff * b_020 + xoff * b_110 + zoff * b_011;
//     b_001 = b_001 + c2o1 * zoff * b_002 + xoff * b_101 + yoff * b_011;
//     c_000 = c_000 + xoff * c_100 + yoff * c_010 + zoff * c_001 + xoff_sq * c_200 + yoff_sq * c_020 + zoff_sq * c_002 +
//             xoff * yoff * c_110 + xoff * zoff * c_101 + yoff * zoff * c_011;
//     c_100 = c_100 + c2o1 * xoff * c_200 + yoff * c_110 + zoff * c_101;
//     c_010 = c_010 + c2o1 * yoff * c_020 + xoff * c_110 + zoff * c_011;
//     c_001 = c_001 + c2o1 * zoff * c_002 + xoff * c_101 + yoff * c_011;
//     d_000 = d_000 + xoff * d_100 + yoff * d_010 + zoff * d_001 + 
//             xoff * yoff * d_110 + xoff * zoff * d_101 + yoff * zoff * d_011;
   
   calcMoments(icell.TSW, omega, press_SWT, vx1_SWT, vx2_SWT, vx3_SWT, kxyFromfcNEQ_SWT, kyzFromfcNEQ_SWT, kxzFromfcNEQ_SWT, kxxMyyFromfcNEQ_SWT, kxxMzzFromfcNEQ_SWT);
   calcMoments(icell.TNW, omega, press_NWT, vx1_NWT, vx2_NWT, vx3_NWT, kxyFromfcNEQ_NWT, kyzFromfcNEQ_NWT, kxzFromfcNEQ_NWT, kxxMyyFromfcNEQ_NWT, kxxMzzFromfcNEQ_NWT);
   calcMoments(icell.TNE, omega, press_NET, vx1_NET, vx2_NET, vx3_NET, kxyFromfcNEQ_NET, kyzFromfcNEQ_NET, kxzFromfcNEQ_NET, kxxMyyFromfcNEQ_NET, kxxMzzFromfcNEQ_NET);
   calcMoments(icell.TSE, omega, press_SET, vx1_SET, vx2_SET, vx3_SET, kxyFromfcNEQ_SET, kyzFromfcNEQ_SET, kxzFromfcNEQ_SET, kxxMyyFromfcNEQ_SET, kxxMzzFromfcNEQ_SET);
   calcMoments(icell.BSW, omega, press_SWB, vx1_SWB, vx2_SWB, vx3_SWB, kxyFromfcNEQ_SWB, kyzFromfcNEQ_SWB, kxzFromfcNEQ_SWB, kxxMyyFromfcNEQ_SWB, kxxMzzFromfcNEQ_SWB);
   calcMoments(icell.BNW, omega, press_NWB, vx1_NWB, vx2_NWB, vx3_NWB, kxyFromfcNEQ_NWB, kyzFromfcNEQ_NWB, kxzFromfcNEQ_NWB, kxxMyyFromfcNEQ_NWB, kxxMzzFromfcNEQ_NWB);
   calcMoments(icell.BNE, omega, press_NEB, vx1_NEB, vx2_NEB, vx3_NEB, kxyFromfcNEQ_NEB, kyzFromfcNEQ_NEB, kxzFromfcNEQ_NEB, kxxMyyFromfcNEQ_NEB, kxxMzzFromfcNEQ_NEB);
   calcMoments(icell.BSE, omega, press_SEB, vx1_SEB, vx2_SEB, vx3_SEB, kxyFromfcNEQ_SEB, kyzFromfcNEQ_SEB, kxzFromfcNEQ_SEB, kxxMyyFromfcNEQ_SEB, kxxMzzFromfcNEQ_SEB);

   a0 = (-kxxMyyFromfcNEQ_NEB - kxxMyyFromfcNEQ_NET + kxxMyyFromfcNEQ_NWB + kxxMyyFromfcNEQ_NWT -
      kxxMyyFromfcNEQ_SEB - kxxMyyFromfcNEQ_SET + kxxMyyFromfcNEQ_SWB + kxxMyyFromfcNEQ_SWT -
      kxxMzzFromfcNEQ_NEB - kxxMzzFromfcNEQ_NET + kxxMzzFromfcNEQ_NWB + kxxMzzFromfcNEQ_NWT -
      kxxMzzFromfcNEQ_SEB - kxxMzzFromfcNEQ_SET + kxxMzzFromfcNEQ_SWB + kxxMzzFromfcNEQ_SWT -
      2.*kxyFromfcNEQ_NEB - 2.*kxyFromfcNEQ_NET - 2.*kxyFromfcNEQ_NWB - 2.*kxyFromfcNEQ_NWT +
      2.*kxyFromfcNEQ_SEB + 2.*kxyFromfcNEQ_SET + 2.*kxyFromfcNEQ_SWB + 2.*kxyFromfcNEQ_SWT +
      2.*kxzFromfcNEQ_NEB - 2.*kxzFromfcNEQ_NET + 2.*kxzFromfcNEQ_NWB - 2.*kxzFromfcNEQ_NWT +
      2.*kxzFromfcNEQ_SEB - 2.*kxzFromfcNEQ_SET + 2.*kxzFromfcNEQ_SWB - 2.*kxzFromfcNEQ_SWT +
      8.*vx1_NEB + 8.*vx1_NET + 8.*vx1_NWB + 8.*vx1_NWT + 8.*vx1_SEB +
      8.*vx1_SET + 8.*vx1_SWB + 8.*vx1_SWT + 2.*vx2_NEB + 2.*vx2_NET -
      2.*vx2_NWB - 2.*vx2_NWT - 2.*vx2_SEB - 2.*vx2_SET + 2.*vx2_SWB +
      2.*vx2_SWT - 2.*vx3_NEB + 2.*vx3_NET + 2.*vx3_NWB - 2.*vx3_NWT -
      2.*vx3_SEB + 2.*vx3_SET + 2.*vx3_SWB - 2.*vx3_SWT)/64.;
   b0 = (2.*kxxMyyFromfcNEQ_NEB + 2.*kxxMyyFromfcNEQ_NET + 2.*kxxMyyFromfcNEQ_NWB + 2.*kxxMyyFromfcNEQ_NWT -
      2.*kxxMyyFromfcNEQ_SEB - 2.*kxxMyyFromfcNEQ_SET - 2.*kxxMyyFromfcNEQ_SWB - 2.*kxxMyyFromfcNEQ_SWT -
      kxxMzzFromfcNEQ_NEB - kxxMzzFromfcNEQ_NET - kxxMzzFromfcNEQ_NWB - kxxMzzFromfcNEQ_NWT +
      kxxMzzFromfcNEQ_SEB + kxxMzzFromfcNEQ_SET + kxxMzzFromfcNEQ_SWB + kxxMzzFromfcNEQ_SWT -
      2.*kxyFromfcNEQ_NEB - 2.*kxyFromfcNEQ_NET + 2.*kxyFromfcNEQ_NWB + 2.*kxyFromfcNEQ_NWT -
      2.*kxyFromfcNEQ_SEB - 2.*kxyFromfcNEQ_SET + 2.*kxyFromfcNEQ_SWB + 2.*kxyFromfcNEQ_SWT +
      2.*kyzFromfcNEQ_NEB - 2.*kyzFromfcNEQ_NET + 2.*kyzFromfcNEQ_NWB - 2.*kyzFromfcNEQ_NWT +
      2.*kyzFromfcNEQ_SEB - 2.*kyzFromfcNEQ_SET + 2.*kyzFromfcNEQ_SWB - 2.*kyzFromfcNEQ_SWT +
      2.*vx1_NEB + 2.*vx1_NET - 2.*vx1_NWB - 2.*vx1_NWT -
      2.*vx1_SEB - 2.*vx1_SET + 2.*vx1_SWB + 2.*vx1_SWT +
      8.*vx2_NEB + 8.*vx2_NET + 8.*vx2_NWB + 8.*vx2_NWT +
      8.*vx2_SEB + 8.*vx2_SET + 8.*vx2_SWB + 8.*vx2_SWT -
      2.*vx3_NEB + 2.*vx3_NET - 2.*vx3_NWB + 2.*vx3_NWT +
      2.*vx3_SEB - 2.*vx3_SET + 2.*vx3_SWB - 2.*vx3_SWT)/64.;
   c0 = (kxxMyyFromfcNEQ_NEB - kxxMyyFromfcNEQ_NET + kxxMyyFromfcNEQ_NWB - kxxMyyFromfcNEQ_NWT +
      kxxMyyFromfcNEQ_SEB - kxxMyyFromfcNEQ_SET + kxxMyyFromfcNEQ_SWB - kxxMyyFromfcNEQ_SWT -
      2.*kxxMzzFromfcNEQ_NEB + 2.*kxxMzzFromfcNEQ_NET - 2.*kxxMzzFromfcNEQ_NWB + 2.*kxxMzzFromfcNEQ_NWT -
      2.*kxxMzzFromfcNEQ_SEB + 2.*kxxMzzFromfcNEQ_SET - 2.*kxxMzzFromfcNEQ_SWB + 2.*kxxMzzFromfcNEQ_SWT -
      2.*kxzFromfcNEQ_NEB - 2.*kxzFromfcNEQ_NET + 2.*kxzFromfcNEQ_NWB + 2.*kxzFromfcNEQ_NWT -
      2.*kxzFromfcNEQ_SEB - 2.*kxzFromfcNEQ_SET + 2.*kxzFromfcNEQ_SWB + 2.*kxzFromfcNEQ_SWT -
      2.*kyzFromfcNEQ_NEB - 2.*kyzFromfcNEQ_NET - 2.*kyzFromfcNEQ_NWB - 2.*kyzFromfcNEQ_NWT +
      2.*kyzFromfcNEQ_SEB + 2.*kyzFromfcNEQ_SET + 2.*kyzFromfcNEQ_SWB + 2.*kyzFromfcNEQ_SWT -
      2.*vx1_NEB + 2.*vx1_NET + 2.*vx1_NWB - 2.*vx1_NWT -
      2.*vx1_SEB + 2.*vx1_SET + 2.*vx1_SWB - 2.*vx1_SWT -
      2.*vx2_NEB + 2.*vx2_NET - 2.*vx2_NWB + 2.*vx2_NWT +
      2.*vx2_SEB - 2.*vx2_SET + 2.*vx2_SWB - 2.*vx2_SWT +
      8.*vx3_NEB + 8.*vx3_NET + 8.*vx3_NWB + 8.*vx3_NWT +
      8.*vx3_SEB + 8.*vx3_SET + 8.*vx3_SWB + 8.*vx3_SWT)/64.;

      
   ax = (vx1_NEB + vx1_NET - vx1_NWB - vx1_NWT + vx1_SEB + vx1_SET - vx1_SWB - vx1_SWT)/4.;
   bx = (vx2_NEB + vx2_NET - vx2_NWB - vx2_NWT + vx2_SEB + vx2_SET - vx2_SWB - vx2_SWT)/4.;
   cx = (vx3_NEB + vx3_NET - vx3_NWB - vx3_NWT + vx3_SEB + vx3_SET - vx3_SWB - vx3_SWT)/4.;
   axx= (kxxMyyFromfcNEQ_NEB + kxxMyyFromfcNEQ_NET - kxxMyyFromfcNEQ_NWB - kxxMyyFromfcNEQ_NWT +
      kxxMyyFromfcNEQ_SEB + kxxMyyFromfcNEQ_SET - kxxMyyFromfcNEQ_SWB - kxxMyyFromfcNEQ_SWT +
      kxxMzzFromfcNEQ_NEB + kxxMzzFromfcNEQ_NET - kxxMzzFromfcNEQ_NWB - kxxMzzFromfcNEQ_NWT +
      kxxMzzFromfcNEQ_SEB + kxxMzzFromfcNEQ_SET - kxxMzzFromfcNEQ_SWB - kxxMzzFromfcNEQ_SWT +
      2.*vx2_NEB + 2.*vx2_NET - 2.*vx2_NWB - 2.*vx2_NWT -
      2.*vx2_SEB - 2.*vx2_SET + 2.*vx2_SWB + 2.*vx2_SWT -
      2.*vx3_NEB + 2.*vx3_NET + 2.*vx3_NWB - 2.*vx3_NWT -
      2.*vx3_SEB + 2.*vx3_SET + 2.*vx3_SWB - 2.*vx3_SWT)/16.;
   bxx= (kxyFromfcNEQ_NEB + kxyFromfcNEQ_NET - kxyFromfcNEQ_NWB - kxyFromfcNEQ_NWT +
      kxyFromfcNEQ_SEB + kxyFromfcNEQ_SET - kxyFromfcNEQ_SWB - kxyFromfcNEQ_SWT -
      2.*vx1_NEB - 2.*vx1_NET + 2.*vx1_NWB + 2.*vx1_NWT +
      2.*vx1_SEB + 2.*vx1_SET - 2.*vx1_SWB - 2.*vx1_SWT)/8.;
   cxx= (kxzFromfcNEQ_NEB + kxzFromfcNEQ_NET - kxzFromfcNEQ_NWB - kxzFromfcNEQ_NWT +
      kxzFromfcNEQ_SEB + kxzFromfcNEQ_SET - kxzFromfcNEQ_SWB - kxzFromfcNEQ_SWT +
      2.*vx1_NEB - 2.*vx1_NET - 2.*vx1_NWB + 2.*vx1_NWT +
      2.*vx1_SEB - 2.*vx1_SET - 2.*vx1_SWB + 2.*vx1_SWT)/8.;
   ay = (vx1_NEB + vx1_NET + vx1_NWB + vx1_NWT - vx1_SEB - vx1_SET - vx1_SWB - vx1_SWT)/4.;
   by = (vx2_NEB + vx2_NET + vx2_NWB + vx2_NWT - vx2_SEB - vx2_SET - vx2_SWB - vx2_SWT)/4.;
   cy = (vx3_NEB + vx3_NET + vx3_NWB + vx3_NWT - vx3_SEB - vx3_SET - vx3_SWB - vx3_SWT)/4.;
   ayy= (kxyFromfcNEQ_NEB + kxyFromfcNEQ_NET + kxyFromfcNEQ_NWB + kxyFromfcNEQ_NWT -
      kxyFromfcNEQ_SEB - kxyFromfcNEQ_SET - kxyFromfcNEQ_SWB - kxyFromfcNEQ_SWT -
      2.*vx2_NEB - 2.*vx2_NET + 2.*vx2_NWB + 2.*vx2_NWT +
      2.*vx2_SEB + 2.*vx2_SET - 2.*vx2_SWB - 2.*vx2_SWT)/8.;
   byy= (-2.*kxxMyyFromfcNEQ_NEB - 2.*kxxMyyFromfcNEQ_NET - 2.*kxxMyyFromfcNEQ_NWB - 2.*kxxMyyFromfcNEQ_NWT +
      2.*kxxMyyFromfcNEQ_SEB + 2.*kxxMyyFromfcNEQ_SET + 2.*kxxMyyFromfcNEQ_SWB + 2.*kxxMyyFromfcNEQ_SWT +
      kxxMzzFromfcNEQ_NEB + kxxMzzFromfcNEQ_NET + kxxMzzFromfcNEQ_NWB + kxxMzzFromfcNEQ_NWT -
      kxxMzzFromfcNEQ_SEB - kxxMzzFromfcNEQ_SET - kxxMzzFromfcNEQ_SWB - kxxMzzFromfcNEQ_SWT +
      2.*vx1_NEB + 2.*vx1_NET - 2.*vx1_NWB - 2.*vx1_NWT -
      2.*vx1_SEB - 2.*vx1_SET + 2.*vx1_SWB + 2.*vx1_SWT -
      2.*vx3_NEB + 2.*vx3_NET - 2.*vx3_NWB + 2.*vx3_NWT +
      2.*vx3_SEB - 2.*vx3_SET + 2.*vx3_SWB - 2.*vx3_SWT)/16.;
   cyy= (kyzFromfcNEQ_NEB + kyzFromfcNEQ_NET + kyzFromfcNEQ_NWB + kyzFromfcNEQ_NWT -
      kyzFromfcNEQ_SEB - kyzFromfcNEQ_SET - kyzFromfcNEQ_SWB - kyzFromfcNEQ_SWT +
      2.*vx2_NEB - 2.*vx2_NET + 2.*vx2_NWB - 2.*vx2_NWT -
      2.*vx2_SEB + 2.*vx2_SET - 2.*vx2_SWB + 2.*vx2_SWT)/8.;
   az = (-vx1_NEB + vx1_NET - vx1_NWB + vx1_NWT - vx1_SEB + vx1_SET - vx1_SWB + vx1_SWT)/4.;
   bz = (-vx2_NEB + vx2_NET - vx2_NWB + vx2_NWT - vx2_SEB + vx2_SET - vx2_SWB + vx2_SWT)/4.;
   cz = (-vx3_NEB + vx3_NET - vx3_NWB + vx3_NWT - vx3_SEB + vx3_SET - vx3_SWB + vx3_SWT)/4.;
   azz= (-kxzFromfcNEQ_NEB + kxzFromfcNEQ_NET - kxzFromfcNEQ_NWB + kxzFromfcNEQ_NWT -
      kxzFromfcNEQ_SEB + kxzFromfcNEQ_SET - kxzFromfcNEQ_SWB + kxzFromfcNEQ_SWT +
      2.*vx3_NEB - 2.*vx3_NET - 2.*vx3_NWB + 2.*vx3_NWT +
      2.*vx3_SEB - 2.*vx3_SET - 2.*vx3_SWB + 2.*vx3_SWT)/8.;
   bzz= (-kyzFromfcNEQ_NEB + kyzFromfcNEQ_NET - kyzFromfcNEQ_NWB + kyzFromfcNEQ_NWT -
      kyzFromfcNEQ_SEB + kyzFromfcNEQ_SET - kyzFromfcNEQ_SWB + kyzFromfcNEQ_SWT +
      2.*vx3_NEB - 2.*vx3_NET + 2.*vx3_NWB - 2.*vx3_NWT -
      2.*vx3_SEB + 2.*vx3_SET - 2.*vx3_SWB + 2.*vx3_SWT)/8.;
   czz= (-kxxMyyFromfcNEQ_NEB + kxxMyyFromfcNEQ_NET - kxxMyyFromfcNEQ_NWB + kxxMyyFromfcNEQ_NWT -
      kxxMyyFromfcNEQ_SEB + kxxMyyFromfcNEQ_SET - kxxMyyFromfcNEQ_SWB + kxxMyyFromfcNEQ_SWT +
      2.*kxxMzzFromfcNEQ_NEB - 2.*kxxMzzFromfcNEQ_NET + 2.*kxxMzzFromfcNEQ_NWB - 2.*kxxMzzFromfcNEQ_NWT +
      2.*kxxMzzFromfcNEQ_SEB - 2.*kxxMzzFromfcNEQ_SET + 2.*kxxMzzFromfcNEQ_SWB - 2.*kxxMzzFromfcNEQ_SWT -
      2.*vx1_NEB + 2.*vx1_NET + 2.*vx1_NWB - 2.*vx1_NWT -
      2.*vx1_SEB + 2.*vx1_SET + 2.*vx1_SWB - 2.*vx1_SWT -
      2.*vx2_NEB + 2.*vx2_NET - 2.*vx2_NWB + 2.*vx2_NWT +
      2.*vx2_SEB - 2.*vx2_SET + 2.*vx2_SWB - 2.*vx2_SWT)/16.;
   axy= (vx1_NEB + vx1_NET - vx1_NWB - vx1_NWT - vx1_SEB - vx1_SET + vx1_SWB + vx1_SWT)/2.;
   bxy= (vx2_NEB + vx2_NET - vx2_NWB - vx2_NWT - vx2_SEB - vx2_SET + vx2_SWB + vx2_SWT)/2.;
   cxy= (vx3_NEB + vx3_NET - vx3_NWB - vx3_NWT - vx3_SEB - vx3_SET + vx3_SWB + vx3_SWT)/2.;
   axz= (-vx1_NEB + vx1_NET + vx1_NWB - vx1_NWT - vx1_SEB + vx1_SET + vx1_SWB - vx1_SWT)/2.;
   bxz= (-vx2_NEB + vx2_NET + vx2_NWB - vx2_NWT - vx2_SEB + vx2_SET + vx2_SWB - vx2_SWT)/2.;
   cxz= (-vx3_NEB + vx3_NET + vx3_NWB - vx3_NWT - vx3_SEB + vx3_SET + vx3_SWB - vx3_SWT)/2.;
   ayz= (-vx1_NEB + vx1_NET - vx1_NWB + vx1_NWT + vx1_SEB - vx1_SET + vx1_SWB - vx1_SWT)/2.;
   byz= (-vx2_NEB + vx2_NET - vx2_NWB + vx2_NWT + vx2_SEB - vx2_SET + vx2_SWB - vx2_SWT)/2.;
   cyz= (-vx3_NEB + vx3_NET - vx3_NWB + vx3_NWT + vx3_SEB - vx3_SET + vx3_SWB - vx3_SWT)/2.;
   axyz=-vx1_NEB + vx1_NET + vx1_NWB - vx1_NWT + vx1_SEB - vx1_SET - vx1_SWB + vx1_SWT;
   bxyz=-vx2_NEB + vx2_NET + vx2_NWB - vx2_NWT + vx2_SEB - vx2_SET - vx2_SWB + vx2_SWT;
   cxyz=-vx3_NEB + vx3_NET + vx3_NWB - vx3_NWT + vx3_SEB - vx3_SET - vx3_SWB + vx3_SWT;


//    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    kxyAverage       =0;//(kxyFromfcNEQ_SWB+
//                        //kxyFromfcNEQ_SWT+
//                        //kxyFromfcNEQ_SET+
//                        //kxyFromfcNEQ_SEB+
//                        //kxyFromfcNEQ_NWB+
//                        //kxyFromfcNEQ_NWT+
//                        //kxyFromfcNEQ_NET+
//                        //kxyFromfcNEQ_NEB)*c1o8-(ay+bx);
//    kyzAverage       =0;//(kyzFromfcNEQ_SWB+
//                        //kyzFromfcNEQ_SWT+
//                        //kyzFromfcNEQ_SET+
//                        //kyzFromfcNEQ_SEB+
//                        //kyzFromfcNEQ_NWB+
//                        //kyzFromfcNEQ_NWT+
//                        //kyzFromfcNEQ_NET+
//                        //kyzFromfcNEQ_NEB)*c1o8-(bz+cy);
//    kxzAverage       =0;//(kxzFromfcNEQ_SWB+
//                        //kxzFromfcNEQ_SWT+
//                        //kxzFromfcNEQ_SET+
//                        //kxzFromfcNEQ_SEB+
//                        //kxzFromfcNEQ_NWB+
//                        //kxzFromfcNEQ_NWT+
//                        //kxzFromfcNEQ_NET+
//                        //kxzFromfcNEQ_NEB)*c1o8-(az+cx);
//    kxxMyyAverage    =0;//(kxxMyyFromfcNEQ_SWB+
//                        //kxxMyyFromfcNEQ_SWT+
//                        //kxxMyyFromfcNEQ_SET+
//                        //kxxMyyFromfcNEQ_SEB+
//                        //kxxMyyFromfcNEQ_NWB+
//                        //kxxMyyFromfcNEQ_NWT+
//                        //kxxMyyFromfcNEQ_NET+
//                        //kxxMyyFromfcNEQ_NEB)*c1o8-(ax-by);
//    kxxMzzAverage    =0;//(kxxMzzFromfcNEQ_SWB+
//                        //kxxMzzFromfcNEQ_SWT+
//                        //kxxMzzFromfcNEQ_SET+
//                        //kxxMzzFromfcNEQ_SEB+
//                        //kxxMzzFromfcNEQ_NWB+
//                        //kxxMzzFromfcNEQ_NWT+
//                        //kxxMzzFromfcNEQ_NET+
//                        //kxxMzzFromfcNEQ_NEB)*c1o8-(ax-cz);
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //
   // Bernd das Brot
   //
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   a0 = a0 + xoff * ax + yoff * ay + zoff * az + xoff_sq * axx + yoff_sq * ayy + zoff_sq * azz + xoff*yoff*axy + xoff*zoff*axz + yoff*zoff*ayz + xoff*yoff*zoff*axyz ;
   ax = ax + 2. * xoff * axx + yoff * axy + zoff * axz + yoff*zoff*axyz;
   ay = ay + 2. * yoff * ayy + xoff * axy + zoff * ayz + xoff*zoff*axyz;
   az = az + 2. * zoff * azz + xoff * axz + yoff * ayz + xoff*yoff*axyz;
   b0 = b0 + xoff * bx + yoff * by + zoff * bz + xoff_sq * bxx + yoff_sq * byy + zoff_sq * bzz + xoff*yoff*bxy + xoff*zoff*bxz + yoff*zoff*byz + xoff*yoff*zoff*bxyz;
   bx = bx + 2. * xoff * bxx + yoff * bxy + zoff * bxz + yoff*zoff*bxyz;
   by = by + 2. * yoff * byy + xoff * bxy + zoff * byz + xoff*zoff*bxyz;
   bz = bz + 2. * zoff * bzz + xoff * bxz + yoff * byz + xoff*yoff*bxyz;
   c0 = c0 + xoff * cx + yoff * cy + zoff * cz + xoff_sq * cxx + yoff_sq * cyy + zoff_sq * czz + xoff*yoff*cxy + xoff*zoff*cxz + yoff*zoff*cyz + xoff*yoff*zoff*cxyz;
   cx = cx + 2. * xoff * cxx + yoff * cxy + zoff * cxz + yoff*zoff*cxyz;
   cy = cy + 2. * yoff * cyy + xoff * cxy + zoff * cyz + xoff*zoff*cxyz;
   cz = cz + 2. * zoff * czz + xoff * cxz + yoff * cyz + xoff*yoff*cxyz;

   axy= axy + zoff*axyz;
   axz= axz + yoff*axyz;
   ayz= ayz + xoff*axyz;
   bxy= bxy + zoff*bxyz;
   bxz= bxz + yoff*bxyz;
   byz= byz + xoff*bxyz;
   cxy= cxy + zoff*cxyz;
   cxz= cxz + yoff*cxyz;
   cyz= cyz + xoff*cxyz;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetMomentsInterpolationProcessor::calcInterpolatedNodeCF(real* f, real omega, real x, real y, real z, real press, real xs, real ys, real zs)
{
   using namespace D3Q27System;
   using namespace vf::lbm::dir;
   using namespace vf::basics::constant;

   real eps_new = c1o2;
   real o = omega;
   //bulk viscosity
   real oP = OxxPyyPzzF;

//   LBMReal rho  = press ;//+ (c2o1*axx*x+axy*y+axz*z+axyz*y*z+ax + c2o1*byy*y+bxy*x+byz*z+bxyz*x*z+by + c2o1*czz*z+cxz*x+cyz*y+cxyz*x*y+cz)/c3o1;
   real vx1  = a0 + c1o4*( xs*ax + ys*ay + zs*az) + c1o16*(axx + xs*ys*axy + xs*zs*axz + ayy + ys*zs*ayz + azz) + c1o64*(xs*ys*zs*axyz);
   real vx2  = b0 + c1o4*( xs*bx + ys*by + zs*bz) + c1o16*(bxx + xs*ys*bxy + xs*zs*bxz + byy + ys*zs*byz + bzz) + c1o64*(xs*ys*zs*bxyz);
   real vx3  = c0 + c1o4*( xs*cx + ys*cy + zs*cz) + c1o16*(cxx + xs*ys*cxy + xs*zs*cxz + cyy + ys*zs*cyz + czz) + c1o64*(xs*ys*zs*cxyz);

   real mfcbb = c0o1;
   real mfabb = c0o1;
   real mfbcb = c0o1;
   real mfbab = c0o1;
   real mfbbc = c0o1;
   real mfbba = c0o1;
   real mfccb = c0o1;
   real mfaab = c0o1;
   real mfcab = c0o1;
   real mfacb = c0o1;
   real mfcbc = c0o1;
   real mfaba = c0o1;
   real mfcba = c0o1;
   real mfabc = c0o1;
   real mfbcc = c0o1;
   real mfbaa = c0o1;
   real mfbca = c0o1;
   real mfbac = c0o1;
   real mfbbb = c0o1;
   real mfccc = c0o1;
   real mfaac = c0o1;
   real mfcac = c0o1;
   real mfacc = c0o1;
   real mfcca = c0o1;
   real mfaaa = c0o1;
   real mfcaa = c0o1;
   real mfaca = c0o1;

   mfaaa = press; // if drho is interpolated directly

   real vx1Sq = vx1*vx1;
   real vx2Sq = vx2*vx2;
   real vx3Sq = vx3*vx3;
   real oMdrho = c1o1;

   //c2o1f

   // linear combinations
   real mxxPyyPzz = mfaaa - c2o3*(ax + by + c2o1 *axx*x + bxy*x + axy*y + c2o1 *byy*y + axz*z + byz*z + bxyz*x*z + axyz*y*z + cz - cxz*x + cyz*y + cxyz*x*y + c2o1 *czz*z)*eps_new / oP* (c1o1 + press);
   real mxxMyy    = -c2o3*(ax - by + kxxMyyAverage + c2o1 *axx*x - bxy*x + axy*y - c2o1 *byy*y + axz*z - byz*z - bxyz*x*z + axyz*y*z)*eps_new/o * (c1o1 + press);
   real mxxMzz    = -c2o3*(ax - cz + kxxMzzAverage + c2o1 *axx*x - cxz*x + axy*y - cyz*y - cxyz*x*y + axz*z - c2o1 *czz*z + axyz*y*z)*eps_new/o * (c1o1 + press);

   mfabb     = -c1o3 * (bz + cy + kyzAverage + bxz*x + cxy*x + byz*y + c2o1 *cyy*y + bxyz*x*y + c2o1 *bzz*z + cyz*z + cxyz*x*z)*eps_new/o * (c1o1 + press);
   mfbab     = -c1o3 * (az + cx + kxzAverage + axz*x + c2o1 *cxx*x + ayz*y + cxy*y + axyz*x*y + c2o1 *azz*z + cxz*z + cxyz*y*z)*eps_new/o * (c1o1 + press);
   mfbba     = -c1o3 * (ay + bx + kxyAverage + axy*x + c2o1 *bxx*x + c2o1 *ayy*y + bxy*y + ayz*z + bxz*z + axyz*x*z + bxyz*y*z)*eps_new/o * (c1o1 + press);

   // linear combinations back
   mfcaa = c1o3 * (mxxMyy +       mxxMzz + mxxPyyPzz) ;
   mfaca = c1o3 * (-c2o1 * mxxMyy +       mxxMzz + mxxPyyPzz) ;
   mfaac = c1o3 * (mxxMyy - c2o1 * mxxMzz + mxxPyyPzz) ;

   //three
   mfbbb = c0o1;
   real mxxyPyzz = c0o1;
   real mxxyMyzz = c0o1;
   real mxxzPyyz = c0o1;
   real mxxzMyyz = c0o1;
   real mxyyPxzz = c0o1;
   real mxyyMxzz = c0o1;

   // linear combinations back
   mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
   mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
   mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
   mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
   mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
   mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

   //4.f
   mfacc = mfaaa*c1o9;
   mfcac = mfacc;
   mfcca = mfacc;

   //5.

   //6.

   mfccc = mfaaa*c1o27;
   ////////////////////////////////////////////////////////////////////////////////////
   //back
   ////////////////////////////////////////////////////////////////////////////////////
   //mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
   ////////////////////////////////////////////////////////////////////////////////////
   // Z - Dir
   real m0 =  mfaac * c1o2 +      mfaab * (vx3 - c1o2) + (mfaaa + c1o1 * oMdrho) * (vx3Sq - vx3) * c1o2;
   real m1 = -mfaac        - c2o1 * mfaab *  vx3         +  mfaaa                * (c1o1 - vx3Sq)              - c1o1 * oMdrho * vx3Sq;
   real m2 =  mfaac * c1o2 +      mfaab * (vx3 + c1o2) + (mfaaa + c1o1 * oMdrho) * (vx3Sq + vx3) * c1o2;
   mfaaa = m0;
   mfaab = m1;
   mfaac = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfabc * c1o2 +      mfabb * (vx3 - c1o2) + mfaba * (vx3Sq - vx3) * c1o2;
   m1 = -mfabc        - c2o1 * mfabb *  vx3         + mfaba * (c1o1 - vx3Sq);
   m2 =  mfabc * c1o2 +      mfabb * (vx3 + c1o2) + mfaba * (vx3Sq + vx3) * c1o2;
   mfaba = m0;
   mfabb = m1;
   mfabc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfacc * c1o2 +      mfacb * (vx3 - c1o2) + (mfaca + c1o3 * oMdrho) * (vx3Sq - vx3) * c1o2;
   m1 = -mfacc        - c2o1 * mfacb *  vx3         +  mfaca                  * (c1o1 - vx3Sq)              - c1o3 * oMdrho * vx3Sq;
   m2 =  mfacc * c1o2 +      mfacb * (vx3 + c1o2) + (mfaca + c1o3 * oMdrho) * (vx3Sq + vx3) * c1o2;
   mfaca = m0;
   mfacb = m1;
   mfacc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfbac * c1o2 +      mfbab * (vx3 - c1o2) + mfbaa * (vx3Sq - vx3) * c1o2;
   m1 = -mfbac        - c2o1 * mfbab *  vx3         + mfbaa * (c1o1 - vx3Sq);
   m2 =  mfbac * c1o2 +      mfbab * (vx3 + c1o2) + mfbaa * (vx3Sq + vx3) * c1o2;
   mfbaa = m0;
   mfbab = m1;
   mfbac = m2;
   /////////b//////////////////////////////////////////////////////////////////////////
   m0 =  mfbbc * c1o2 +      mfbbb * (vx3 - c1o2) + mfbba * (vx3Sq - vx3) * c1o2;
   m1 = -mfbbc        - c2o1 * mfbbb *  vx3         + mfbba * (c1o1 - vx3Sq);
   m2 =  mfbbc * c1o2 +      mfbbb * (vx3 + c1o2) + mfbba * (vx3Sq + vx3) * c1o2;
   mfbba = m0;
   mfbbb = m1;
   mfbbc = m2;
   /////////b//////////////////////////////////////////////////////////////////////////
   m0 =  mfbcc * c1o2 +      mfbcb * (vx3 - c1o2) + mfbca * (vx3Sq - vx3) * c1o2;
   m1 = -mfbcc        - c2o1 * mfbcb *  vx3         + mfbca * (c1o1 - vx3Sq);
   m2 =  mfbcc * c1o2 +      mfbcb * (vx3 + c1o2) + mfbca * (vx3Sq + vx3) * c1o2;
   mfbca = m0;
   mfbcb = m1;
   mfbcc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcac * c1o2 +      mfcab * (vx3 - c1o2) + (mfcaa + c1o3 * oMdrho) * (vx3Sq - vx3) * c1o2;
   m1 = -mfcac        - c2o1 * mfcab *  vx3         +  mfcaa                  * (c1o1 - vx3Sq)              - c1o3 * oMdrho * vx3Sq;
   m2 =  mfcac * c1o2 +      mfcab * (vx3 + c1o2) + (mfcaa + c1o3 * oMdrho) * (vx3Sq + vx3) * c1o2;
   mfcaa = m0;
   mfcab = m1;
   mfcac = m2;
   /////////c//////////////////////////////////////////////////////////////////////////
   m0 =  mfcbc * c1o2 +      mfcbb * (vx3 - c1o2) + mfcba * (vx3Sq - vx3) * c1o2;
   m1 = -mfcbc        - c2o1 * mfcbb *  vx3         + mfcba * (c1o1 - vx3Sq);
   m2 =  mfcbc * c1o2 +      mfcbb * (vx3 + c1o2) + mfcba * (vx3Sq + vx3) * c1o2;
   mfcba = m0;
   mfcbb = m1;
   mfcbc = m2;
   /////////c//////////////////////////////////////////////////////////////////////////
   m0 =  mfccc * c1o2 +      mfccb * (vx3 - c1o2) + (mfcca + c1o9 * oMdrho) * (vx3Sq - vx3) * c1o2;
   m1 = -mfccc        - c2o1 * mfccb *  vx3         +  mfcca                  * (c1o1 - vx3Sq)              - c1o9 * oMdrho * vx3Sq;
   m2 =  mfccc * c1o2 +      mfccb * (vx3 + c1o2) + (mfcca + c1o9 * oMdrho) * (vx3Sq + vx3) * c1o2;
   mfcca = m0;
   mfccb = m1;
   mfccc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   //mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
   ////////////////////////////////////////////////////////////////////////////////////
   // Y - Dir
   m0 =  mfaca * c1o2 +      mfaba * (vx2 - c1o2) + (mfaaa + c1o6 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfaca        - c2o1 * mfaba *  vx2         +  mfaaa                  * (c1o1 - vx2Sq)              - c1o6 * oMdrho * vx2Sq;
   m2 =  mfaca * c1o2 +      mfaba * (vx2 + c1o2) + (mfaaa + c1o6 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfaaa = m0;
   mfaba = m1;
   mfaca = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfacb * c1o2 +      mfabb * (vx2 - c1o2) + (mfaab + c2o3 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfacb        - c2o1 * mfabb *  vx2         +  mfaab                  * (c1o1 - vx2Sq)              - c2o3 * oMdrho * vx2Sq;
   m2 =  mfacb * c1o2 +      mfabb * (vx2 + c1o2) + (mfaab + c2o3 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfaab = m0;
   mfabb = m1;
   mfacb = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfacc * c1o2 +      mfabc * (vx2 - c1o2) + (mfaac + c1o6 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfacc        - c2o1 * mfabc *  vx2         +  mfaac                  * (c1o1 - vx2Sq)              - c1o6 * oMdrho * vx2Sq;
   m2 =  mfacc * c1o2 +      mfabc * (vx2 + c1o2) + (mfaac + c1o6 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfaac = m0;
   mfabc = m1;
   mfacc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfbca * c1o2 +      mfbba * (vx2 - c1o2) + mfbaa * (vx2Sq - vx2) * c1o2;
   m1 = -mfbca        - c2o1 * mfbba *  vx2         + mfbaa * (c1o1 - vx2Sq);
   m2 =  mfbca * c1o2 +      mfbba * (vx2 + c1o2) + mfbaa * (vx2Sq + vx2) * c1o2;
   mfbaa = m0;
   mfbba = m1;
   mfbca = m2;
   /////////b//////////////////////////////////////////////////////////////////////////
   m0 =  mfbcb * c1o2 +      mfbbb * (vx2 - c1o2) + mfbab * (vx2Sq - vx2) * c1o2;
   m1 = -mfbcb        - c2o1 * mfbbb *  vx2         + mfbab * (c1o1 - vx2Sq);
   m2 =  mfbcb * c1o2 +      mfbbb * (vx2 + c1o2) + mfbab * (vx2Sq + vx2) * c1o2;
   mfbab = m0;
   mfbbb = m1;
   mfbcb = m2;
   /////////b//////////////////////////////////////////////////////////////////////////
   m0 =  mfbcc * c1o2 +      mfbbc * (vx2 - c1o2) + mfbac * (vx2Sq - vx2) * c1o2;
   m1 = -mfbcc        - c2o1 * mfbbc *  vx2         + mfbac * (c1o1 - vx2Sq);
   m2 =  mfbcc * c1o2 +      mfbbc * (vx2 + c1o2) + mfbac * (vx2Sq + vx2) * c1o2;
   mfbac = m0;
   mfbbc = m1;
   mfbcc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcca * c1o2 +      mfcba * (vx2 - c1o2) + (mfcaa + c1o18 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfcca        - c2o1 * mfcba *  vx2         +  mfcaa                   * (c1o1 - vx2Sq)              - c1o18 * oMdrho * vx2Sq;
   m2 =  mfcca * c1o2 +      mfcba * (vx2 + c1o2) + (mfcaa + c1o18 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfcaa = m0;
   mfcba = m1;
   mfcca = m2;
   /////////c//////////////////////////////////////////////////////////////////////////
   m0 =  mfccb * c1o2 +      mfcbb * (vx2 - c1o2) + (mfcab + c2o9 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfccb        - c2o1 * mfcbb *  vx2         +  mfcab                  * (c1o1 - vx2Sq)              - c2o9 * oMdrho * vx2Sq;
   m2 =  mfccb * c1o2 +      mfcbb * (vx2 + c1o2) + (mfcab + c2o9 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfcab = m0;
   mfcbb = m1;
   mfccb = m2;
   /////////c//////////////////////////////////////////////////////////////////////////
   m0 =  mfccc * c1o2 +      mfcbc * (vx2 - c1o2) + (mfcac + c1o18 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfccc        - c2o1 * mfcbc *  vx2         +  mfcac                   * (c1o1 - vx2Sq)              - c1o18 * oMdrho * vx2Sq;
   m2 =  mfccc * c1o2 +      mfcbc * (vx2 + c1o2) + (mfcac + c1o18 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfcac = m0;
   mfcbc = m1;
   mfccc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   //mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
   ////////////////////////////////////////////////////////////////////////////////////
   // X - Dir
   m0 =  mfcaa * c1o2 +      mfbaa * (vx1 - c1o2) + (mfaaa + c1o36 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcaa        - c2o1 * mfbaa *  vx1         +  mfaaa                   * (c1o1 - vx1Sq)              - c1o36 * oMdrho * vx1Sq;
   m2 =  mfcaa * c1o2 +      mfbaa * (vx1 + c1o2) + (mfaaa + c1o36 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfaaa = m0;
   mfbaa = m1;
   mfcaa = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcba * c1o2 +      mfbba * (vx1 - c1o2) + (mfaba + c1o9 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcba        - c2o1 * mfbba *  vx1         +  mfaba                  * (c1o1 - vx1Sq)              - c1o9 * oMdrho * vx1Sq;
   m2 =  mfcba * c1o2 +      mfbba * (vx1 + c1o2) + (mfaba + c1o9 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfaba = m0;
   mfbba = m1;
   mfcba = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcca * c1o2 +      mfbca * (vx1 - c1o2) + (mfaca + c1o36 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcca        - c2o1 * mfbca *  vx1         +  mfaca                   * (c1o1 - vx1Sq)              - c1o36 * oMdrho * vx1Sq;
   m2 =  mfcca * c1o2 +      mfbca * (vx1 + c1o2) + (mfaca + c1o36 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfaca = m0;
   mfbca = m1;
   mfcca = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcab * c1o2 +      mfbab * (vx1 - c1o2) + (mfaab + c1o9 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcab        - c2o1 * mfbab *  vx1         +  mfaab                  * (c1o1 - vx1Sq)              - c1o9 * oMdrho * vx1Sq;
   m2 =  mfcab * c1o2 +      mfbab * (vx1 + c1o2) + (mfaab + c1o9 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfaab = m0;
   mfbab = m1;
   mfcab = m2;
   ///////////b////////////////////////////////////////////////////////////////////////
   m0 =  mfcbb * c1o2 +      mfbbb * (vx1 - c1o2) + (mfabb + c4o9 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcbb        - c2o1 * mfbbb *  vx1         +  mfabb                  * (c1o1 - vx1Sq)              - c4o9 * oMdrho * vx1Sq;
   m2 =  mfcbb * c1o2 +      mfbbb * (vx1 + c1o2) + (mfabb + c4o9 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfabb = m0;
   mfbbb = m1;
   mfcbb = m2;
   ///////////b////////////////////////////////////////////////////////////////////////
   m0 =  mfccb * c1o2 +      mfbcb * (vx1 - c1o2) + (mfacb + c1o9 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfccb        - c2o1 * mfbcb *  vx1         +  mfacb                  * (c1o1 - vx1Sq)              - c1o9 * oMdrho * vx1Sq;
   m2 =  mfccb * c1o2 +      mfbcb * (vx1 + c1o2) + (mfacb + c1o9 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfacb = m0;
   mfbcb = m1;
   mfccb = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcac * c1o2 +      mfbac * (vx1 - c1o2) + (mfaac + c1o36 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcac        - c2o1 * mfbac *  vx1         +  mfaac                   * (c1o1 - vx1Sq)              - c1o36 * oMdrho * vx1Sq;
   m2 =  mfcac * c1o2 +      mfbac * (vx1 + c1o2) + (mfaac + c1o36 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfaac = m0;
   mfbac = m1;
   mfcac = m2;
   ///////////c////////////////////////////////////////////////////////////////////////
   m0 =  mfcbc * c1o2 +      mfbbc * (vx1 - c1o2) + (mfabc + c1o9 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcbc        - c2o1 * mfbbc *  vx1         +  mfabc                  * (c1o1 - vx1Sq)              - c1o9 * oMdrho * vx1Sq;
   m2 =  mfcbc * c1o2 +      mfbbc * (vx1 + c1o2) + (mfabc + c1o9 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfabc = m0;
   mfbbc = m1;
   mfcbc = m2;
   ///////////c////////////////////////////////////////////////////////////////////////
   m0 =  mfccc * c1o2 +      mfbcc * (vx1 - c1o2) + (mfacc + c1o36 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfccc        - c2o1 * mfbcc *  vx1         +  mfacc                   * (c1o1 - vx1Sq)              - c1o36 * oMdrho * vx1Sq;
   m2 =  mfccc * c1o2 +      mfbcc * (vx1 + c1o2) + (mfacc + c1o36 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfacc = m0;
   mfbcc = m1;
   mfccc = m2;
   ////////////////////////////////////////////////////////////////////////////////////

   f[DIR_P00] = mfcbb;
   f[DIR_M00] = mfabb;
   f[DIR_0P0] = mfbcb;
   f[DIR_0M0] = mfbab;
   f[DIR_00P] = mfbbc;
   f[DIR_00M] = mfbba;
   f[DIR_PP0] = mfccb;
   f[DIR_MM0] = mfaab;
   f[DIR_PM0] = mfcab;
   f[DIR_MP0] = mfacb;
   f[DIR_P0P] = mfcbc;
   f[DIR_M0M] = mfaba;
   f[DIR_P0M] = mfcba;
   f[DIR_M0P] = mfabc;
   f[DIR_0PP] = mfbcc;
   f[DIR_0MM] = mfbaa;
   f[DIR_0PM] = mfbca;
   f[DIR_0MP] = mfbac;
   f[DIR_000] = mfbbb;
   f[DIR_PPP] = mfccc;
   f[DIR_PMP] = mfcac;
   f[DIR_PPM] = mfcca;
   f[DIR_PMM] = mfcaa;
   f[DIR_MPP] = mfacc;
   f[DIR_MMP] = mfaac;
   f[DIR_MPM] = mfaca;
   f[DIR_MMM] = mfaaa;
}
//////////////////////////////////////////////////////////////////////////
//Position SWB -0.25, -0.25, -0.25
real CompressibleOffsetMomentsInterpolationProcessor::calcPressBSW()
{
   return   press_SWT * (c9o64 + c3o16 * xoff + c3o16 * yoff - c9o16 * zoff) +
      press_NWT * (c3o64 + c1o16 * xoff - c3o16 * yoff - c3o16 * zoff) +
      press_SET * (c3o64 - c3o16 * xoff + c1o16 * yoff - c3o16 * zoff) +
      press_NET * (c1o64 - c1o16 * xoff - c1o16 * yoff - c1o16 * zoff) +
      press_NEB * (c3o64 - c3o16 * xoff - c3o16 * yoff + c1o16 * zoff) +
      press_NWB * (c9o64 + c3o16 * xoff - c9o16 * yoff + c3o16 * zoff) +
      press_SEB * (c9o64 - c9o16 * xoff + c3o16 * yoff + c3o16 * zoff) +
      press_SWB * (c27o64 + c9o16 * xoff + c9o16 * yoff + c9o16 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position SWT -0.25, -0.25, 0.25
real CompressibleOffsetMomentsInterpolationProcessor::calcPressTSW()
{
   return   press_SWT * (c27o64 + c9o16 * xoff + c9o16 * yoff - c9o16 * zoff) +
      press_NWT * (c9o64 + c3o16 * xoff - c9o16 * yoff - c3o16 * zoff) +
      press_SET * (c9o64 - c9o16 * xoff + c3o16 * yoff - c3o16 * zoff) +
      press_NET * (c3o64 - c3o16 * xoff - c3o16 * yoff - c1o16 * zoff) +
      press_NEB * (c1o64 - c1o16 * xoff - c1o16 * yoff + c1o16 * zoff) +
      press_NWB * (c3o64 + c1o16 * xoff - c3o16 * yoff + c3o16 * zoff) +
      press_SEB * (c3o64 - c3o16 * xoff + c1o16 * yoff + c3o16 * zoff) +
      press_SWB * (c9o64 + c3o16 * xoff + c3o16 * yoff + c9o16 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position SET 0.25, -0.25, 0.25
real CompressibleOffsetMomentsInterpolationProcessor::calcPressTSE()
{
   return   press_SET * (c27o64 - c9o16 * xoff + c9o16 * yoff - c9o16 * zoff) +
      press_NET * (c9o64 - c3o16 * xoff - c9o16 * yoff - c3o16 * zoff) +
      press_SWT * (c9o64 + c9o16 * xoff + c3o16 * yoff - c3o16 * zoff) +
      press_NWT * (c3o64 + c3o16 * xoff - c3o16 * yoff - c1o16 * zoff) +
      press_NWB * (c1o64 + c1o16 * xoff - c1o16 * yoff + c1o16 * zoff) +
      press_NEB * (c3o64 - c1o16 * xoff - c3o16 * yoff + c3o16 * zoff) +
      press_SWB * (c3o64 + c3o16 * xoff + c1o16 * yoff + c3o16 * zoff) +
      press_SEB * (c9o64 - c3o16 * xoff + c3o16 * yoff + c9o16 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position SEB 0.25, -0.25, -0.25
real CompressibleOffsetMomentsInterpolationProcessor::calcPressBSE()
{
   return   press_SET * (c9o64 - c3o16 * xoff + c3o16 * yoff - c9o16 * zoff) +
      press_NET * (c3o64 - c1o16 * xoff - c3o16 * yoff - c3o16 * zoff) +
      press_SWT * (c3o64 + c3o16 * xoff + c1o16 * yoff - c3o16 * zoff) +
      press_NWT * (c1o64 + c1o16 * xoff - c1o16 * yoff - c1o16 * zoff) +
      press_NWB * (c3o64 + c3o16 * xoff - c3o16 * yoff + c1o16 * zoff) +
      press_NEB * (c9o64 - c3o16 * xoff - c9o16 * yoff + c3o16 * zoff) +
      press_SWB * (c9o64 + c9o16 * xoff + c3o16 * yoff + c3o16 * zoff) +
      press_SEB * (c27o64 - c9o16 * xoff + c9o16 * yoff + c9o16 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position NWB -0.25, 0.25, -0.25
real CompressibleOffsetMomentsInterpolationProcessor::calcPressBNW()
{
   return   press_NWT * (c9o64 + c3o16 * xoff - c3o16 * yoff - c9o16 * zoff) +
      press_NET * (c3o64 - c3o16 * xoff - c1o16 * yoff - c3o16 * zoff) +
      press_SWT * (c3o64 + c1o16 * xoff + c3o16 * yoff - c3o16 * zoff) +
      press_SET * (c1o64 - c1o16 * xoff + c1o16 * yoff - c1o16 * zoff) +
      press_SEB * (c3o64 - c3o16 * xoff + c3o16 * yoff + c1o16 * zoff) +
      press_NEB * (c9o64 - c9o16 * xoff - c3o16 * yoff + c3o16 * zoff) +
      press_SWB * (c9o64 + c3o16 * xoff + c9o16 * yoff + c3o16 * zoff) +
      press_NWB * (c27o64 + c9o16 * xoff - c9o16 * yoff + c9o16 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position NWT -0.25, 0.25, 0.25
real CompressibleOffsetMomentsInterpolationProcessor::calcPressTNW()
{
   return   press_NWT * (c27o64 + c9o16 * xoff - c9o16 * yoff - c9o16 * zoff) +
      press_NET * (c9o64 - c9o16 * xoff - c3o16 * yoff - c3o16 * zoff) +
      press_SWT * (c9o64 + c3o16 * xoff + c9o16 * yoff - c3o16 * zoff) +
      press_SET * (c3o64 - c3o16 * xoff + c3o16 * yoff - c1o16 * zoff) +
      press_SEB * (c1o64 - c1o16 * xoff + c1o16 * yoff + c1o16 * zoff) +
      press_NEB * (c3o64 - c3o16 * xoff - c1o16 * yoff + c3o16 * zoff) +
      press_SWB * (c3o64 + c1o16 * xoff + c3o16 * yoff + c3o16 * zoff) +
      press_NWB * (c9o64 + c3o16 * xoff - c3o16 * yoff + c9o16 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position NET 0.25, 0.25, 0.25
real CompressibleOffsetMomentsInterpolationProcessor::calcPressTNE()
{
   return   press_NET * (c27o64 - c9o16 * xoff - c9o16 * yoff - c9o16 * zoff) +
      press_NWT * (c9o64 + c9o16 * xoff - c3o16 * yoff - c3o16 * zoff) +
      press_SET * (c9o64 - c3o16 * xoff + c9o16 * yoff - c3o16 * zoff) +
      press_SWT * (c3o64 + c3o16 * xoff + c3o16 * yoff - c1o16 * zoff) +
      press_SWB * (c1o64 + c1o16 * xoff + c1o16 * yoff + c1o16 * zoff) +
      press_NWB * (c3o64 + c3o16 * xoff - c1o16 * yoff + c3o16 * zoff) +
      press_SEB * (c3o64 - c1o16 * xoff + c3o16 * yoff + c3o16 * zoff) +
      press_NEB * (c9o64 - c3o16 * xoff - c3o16 * yoff + c9o16 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position NEB 0.25, 0.25, -0.25
real CompressibleOffsetMomentsInterpolationProcessor::calcPressBNE()
{
   return   press_NET * (c9o64 - c3o16 * xoff - c3o16 * yoff - c9o16 * zoff) +
      press_NWT * (c3o64 + c3o16 * xoff - c1o16 * yoff - c3o16 * zoff) +
      press_SET * (c3o64 - c1o16 * xoff + c3o16 * yoff - c3o16 * zoff) +
      press_SWT * (c1o64 + c1o16 * xoff + c1o16 * yoff - c1o16 * zoff) +
      press_SWB * (c3o64 + c3o16 * xoff + c3o16 * yoff + c1o16 * zoff) +
      press_NWB * (c9o64 + c9o16 * xoff - c3o16 * yoff + c3o16 * zoff) +
      press_SEB * (c9o64 - c3o16 * xoff + c9o16 * yoff + c3o16 * zoff) +
      press_NEB * (c27o64 - c9o16 * xoff - c9o16 * yoff + c9o16 * zoff);
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetMomentsInterpolationProcessor::calcInterpolatedVelocity(real x, real y, real z, real& vx1, real& vx2, real& vx3)
{
	vx1  = a0 + ax*x + ay*y + az*z + axx*x*x + ayy*y*y + azz*z*z + axy*x*y + axz*x*z + ayz*y*z+axyz*x*y*z;
	vx2  = b0 + bx*x + by*y + bz*z + bxx*x*x + byy*y*y + bzz*z*z + bxy*x*y + bxz*x*z + byz*y*z+bxyz*x*y*z;
	vx3  = c0 + cx*x + cy*y + cz*z + cxx*x*x + cyy*y*y + czz*z*z + cxy*x*y + cxz*x*z + cyz*y*z+cxyz*x*y*z;
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetMomentsInterpolationProcessor::calcInterpolatedShearStress(real x, real y, real z,real& tauxx, real& tauyy, real& tauzz,real& tauxy, real& tauxz, real& tauyz)
{
	tauxx=ax+c2o1*axx*x+axy*y+axz*z+axyz*y*z;
	tauyy=by+c2o1*byy*y+bxy*x+byz*z+bxyz*x*z;
	tauzz=cz+c2o1*czz*z+cxz*x+cyz*y+cxyz*x*y;
	tauxy= c1o2*((ay+c2o1*ayy*y+axy*x+ayz*z+axyz*x*z)+(bx+c2o1*bxx*x+bxy*y+bxz*z+bxyz*y*z));
	tauxz= c1o2*((az+c2o1*azz*z+axz*x+ayz*y+axyz*x*y)+(cx+c2o1*cxx*x+cxy*y+cxz*z+cxyz*y*z));
	tauyz= c1o2*((bz+c2o1*bzz*z+bxz*x+byz*y+bxyz*x*y)+(cy+c2o1*cyy*y+cxy*x+cyz*z+cxyz*x*z));
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetMomentsInterpolationProcessor::setBulkViscosity(real shearViscosity, real bulkViscosity)
{
   this->shearViscosity = shearViscosity;
   this->bulkViscosity  = bulkViscosity;
}

