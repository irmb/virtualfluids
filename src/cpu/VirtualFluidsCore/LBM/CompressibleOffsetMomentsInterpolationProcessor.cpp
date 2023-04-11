#include "CompressibleOffsetMomentsInterpolationProcessor.h"
#include "D3Q27System.h"

#include <lbm/Scaling.h>
#include <lbm/Interpolation_FC.h>

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
   calcInterpolatedCoefficiets(icellC, omegaC, c1o2);
   calcInterpolatedNodeCF(icellF.BSW, omegaF, -c1o4, -c1o4, -c1o4, calcPressBSW(), -c1o1, -c1o1, -c1o1);
   calcInterpolatedNodeCF(icellF.BNE, omegaF,  c1o4,  c1o4, -c1o4, calcPressBNE(),  c1o1,  c1o1, -c1o1);
   calcInterpolatedNodeCF(icellF.TNW, omegaF, -c1o4,  c1o4,  c1o4, calcPressTNW(), -c1o1,  c1o1,  c1o1);
   calcInterpolatedNodeCF(icellF.TSE, omegaF,  c1o4, -c1o4,  c1o4, calcPressTSE(),  c1o1, -c1o1,  c1o1);
   calcInterpolatedNodeCF(icellF.BNW, omegaF, -c1o4,  c1o4, -c1o4, calcPressBNW(), -c1o1,  c1o1, -c1o1);
   calcInterpolatedNodeCF(icellF.BSE, omegaF,  c1o4, -c1o4, -c1o4, calcPressBSE(),  c1o1, -c1o1, -c1o1);
   calcInterpolatedNodeCF(icellF.TSW, omegaF, -c1o4, -c1o4,  c1o4, calcPressTSW(), -c1o1, -c1o1,  c1o1);
   calcInterpolatedNodeCF(icellF.TNE, omegaF,  c1o4,  c1o4,  c1o4, calcPressTNE(),  c1o1,  c1o1,  c1o1);
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetMomentsInterpolationProcessor::interpolateFineToCoarse(D3Q27ICell& icellF, real* icellC, real xoff, real yoff, real zoff)
{
   setOffsets(xoff, yoff, zoff);
   calcInterpolation(icellF, icellC);
}

void CompressibleOffsetMomentsInterpolationProcessor::calcInterpolation(const D3Q27ICell& icell, real* icellC)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - declare local variables for source nodes
    //!
    const real eps_new = c2o1; // ratio of grid resolutions

    // zeroth and first order moments at the source nodes
    real drho_PPP, vx1_PPP, vx2_PPP, vx3_PPP;
    real drho_MPP, vx1_MPP, vx2_MPP, vx3_MPP;
    real drho_PMP, vx1_PMP, vx2_PMP, vx3_PMP;
    real drho_MMP, vx1_MMP, vx2_MMP, vx3_MMP;
    real drho_PPM, vx1_PPM, vx2_PPM, vx3_PPM;
    real drho_MPM, vx1_MPM, vx2_MPM, vx3_MPM;
    real drho_PMM, vx1_PMM, vx2_PMM, vx3_PMM;
    real drho_MMM, vx1_MMM, vx2_MMM, vx3_MMM;

    // second order moments at the source nodes
    real kxyFromfcNEQ_PPP, kyzFromfcNEQ_PPP, kxzFromfcNEQ_PPP, kxxMyyFromfcNEQ_PPP, kxxMzzFromfcNEQ_PPP;
    real kxyFromfcNEQ_MPP, kyzFromfcNEQ_MPP, kxzFromfcNEQ_MPP, kxxMyyFromfcNEQ_MPP, kxxMzzFromfcNEQ_MPP;
    real kxyFromfcNEQ_PMP, kyzFromfcNEQ_PMP, kxzFromfcNEQ_PMP, kxxMyyFromfcNEQ_PMP, kxxMzzFromfcNEQ_PMP;
    real kxyFromfcNEQ_MMP, kyzFromfcNEQ_MMP, kxzFromfcNEQ_MMP, kxxMyyFromfcNEQ_MMP, kxxMzzFromfcNEQ_MMP;
    real kxyFromfcNEQ_PPM, kyzFromfcNEQ_PPM, kxzFromfcNEQ_PPM, kxxMyyFromfcNEQ_PPM, kxxMzzFromfcNEQ_PPM;
    real kxyFromfcNEQ_MPM, kyzFromfcNEQ_MPM, kxzFromfcNEQ_MPM, kxxMyyFromfcNEQ_MPM, kxxMzzFromfcNEQ_MPM;
    real kxyFromfcNEQ_PMM, kyzFromfcNEQ_PMM, kxzFromfcNEQ_PMM, kxxMyyFromfcNEQ_PMM, kxxMzzFromfcNEQ_PMM;
    real kxyFromfcNEQ_MMM, kyzFromfcNEQ_MMM, kxzFromfcNEQ_MMM, kxxMyyFromfcNEQ_MMM, kxxMzzFromfcNEQ_MMM;

    vf::lbm::calculateMomentsOnSourceNodes(icell.BSW, omegaF, drho_MMM, vx1_MMM, vx2_MMM, vx3_MMM,kxyFromfcNEQ_MMM, kyzFromfcNEQ_MMM, kxzFromfcNEQ_MMM, kxxMyyFromfcNEQ_MMM, kxxMzzFromfcNEQ_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(icell.TSW, omegaF, drho_MMP, vx1_MMP, vx2_MMP, vx3_MMP,kxyFromfcNEQ_MMP, kyzFromfcNEQ_MMP, kxzFromfcNEQ_MMP, kxxMyyFromfcNEQ_MMP, kxxMzzFromfcNEQ_MMP);
    vf::lbm::calculateMomentsOnSourceNodes(icell.TNW, omegaF, drho_MPP, vx1_MPP, vx2_MPP, vx3_MPP,kxyFromfcNEQ_MPP, kyzFromfcNEQ_MPP, kxzFromfcNEQ_MPP, kxxMyyFromfcNEQ_MPP, kxxMzzFromfcNEQ_MPP);
    vf::lbm::calculateMomentsOnSourceNodes(icell.BNW, omegaF, drho_MPM, vx1_MPM, vx2_MPM, vx3_MPM,kxyFromfcNEQ_MPM, kyzFromfcNEQ_MPM, kxzFromfcNEQ_MPM, kxxMyyFromfcNEQ_MPM, kxxMzzFromfcNEQ_MPM);
    vf::lbm::calculateMomentsOnSourceNodes(icell.BSE, omegaF, drho_PMM, vx1_PMM, vx2_PMM, vx3_PMM,kxyFromfcNEQ_PMM, kyzFromfcNEQ_PMM, kxzFromfcNEQ_PMM, kxxMyyFromfcNEQ_PMM, kxxMzzFromfcNEQ_PMM);
    vf::lbm::calculateMomentsOnSourceNodes(icell.TNE, omegaF, drho_PPP, vx1_PPP, vx2_PPP, vx3_PPP,kxyFromfcNEQ_PPP, kyzFromfcNEQ_PPP, kxzFromfcNEQ_PPP, kxxMyyFromfcNEQ_PPP, kxxMzzFromfcNEQ_PPP);
    vf::lbm::calculateMomentsOnSourceNodes(icell.TSE, omegaF, drho_PMP, vx1_PMP, vx2_PMP, vx3_PMP,kxyFromfcNEQ_PMP, kyzFromfcNEQ_PMP, kxzFromfcNEQ_PMP, kxxMyyFromfcNEQ_PMP, kxxMzzFromfcNEQ_PMP);
    vf::lbm::calculateMomentsOnSourceNodes(icell.BNE, omegaF, drho_PPM, vx1_PPM, vx2_PPM, vx3_PPM,kxyFromfcNEQ_PPM, kyzFromfcNEQ_PPM, kxzFromfcNEQ_PPM, kxxMyyFromfcNEQ_PPM, kxxMzzFromfcNEQ_PPM);

    vf::lbm::interpolate_fc(icellC,
        eps_new,
        omegaC,
        xoff,
        yoff,
        zoff,
        xoff_sq,
        yoff_sq,
        zoff_sq,
        drho_PPP, vx1_PPP, vx2_PPP, vx3_PPP,
        drho_MPP, vx1_MPP, vx2_MPP, vx3_MPP,
        drho_PMP, vx1_PMP, vx2_PMP, vx3_PMP,
        drho_MMP, vx1_MMP, vx2_MMP, vx3_MMP,
        drho_PPM, vx1_PPM, vx2_PPM, vx3_PPM,
        drho_MPM, vx1_MPM, vx2_MPM, vx3_MPM,
        drho_PMM, vx1_PMM, vx2_PMM, vx3_PMM,
        drho_MMM, vx1_MMM, vx2_MMM, vx3_MMM,
        kxyFromfcNEQ_PPP, kyzFromfcNEQ_PPP, kxzFromfcNEQ_PPP, kxxMyyFromfcNEQ_PPP, kxxMzzFromfcNEQ_PPP,
        kxyFromfcNEQ_MPP, kyzFromfcNEQ_MPP, kxzFromfcNEQ_MPP, kxxMyyFromfcNEQ_MPP, kxxMzzFromfcNEQ_MPP,
        kxyFromfcNEQ_PMP, kyzFromfcNEQ_PMP, kxzFromfcNEQ_PMP, kxxMyyFromfcNEQ_PMP, kxxMzzFromfcNEQ_PMP,
        kxyFromfcNEQ_MMP, kyzFromfcNEQ_MMP, kxzFromfcNEQ_MMP, kxxMyyFromfcNEQ_MMP, kxxMzzFromfcNEQ_MMP,
        kxyFromfcNEQ_PPM, kyzFromfcNEQ_PPM, kxzFromfcNEQ_PPM, kxxMyyFromfcNEQ_PPM, kxxMzzFromfcNEQ_PPM,
        kxyFromfcNEQ_MPM, kyzFromfcNEQ_MPM, kxzFromfcNEQ_MPM, kxxMyyFromfcNEQ_MPM, kxxMzzFromfcNEQ_MPM,
        kxyFromfcNEQ_PMM, kyzFromfcNEQ_PMM, kxzFromfcNEQ_PMM, kxxMyyFromfcNEQ_PMM, kxxMzzFromfcNEQ_PMM,
        kxyFromfcNEQ_MMM, kyzFromfcNEQ_MMM, kxzFromfcNEQ_MMM, kxxMyyFromfcNEQ_MMM, kxxMzzFromfcNEQ_MMM);
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

