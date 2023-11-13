#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;


__global__ void LB_Init_Comp_AD_27(unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned int* geoD,
    real* Conc,
    real* ux,
    real* uy,
    real* uz,
    unsigned int size_Mat,
    real* DD27,
    bool EvenOrOdd)
{
    ////////////////////////////////////////////////////////////////////////////////
    const unsigned  x = threadIdx.x;  // Globaler x-Index 
    const unsigned  y = blockIdx.x;   // Globaler y-Index 
    const unsigned  z = blockIdx.y;   // Globaler z-Index 

    const unsigned nx = blockDim.x;
    const unsigned ny = gridDim.x;

    const unsigned k = nx*(ny*z + y) + x;
    //////////////////////////////////////////////////////////////////////////

    if (k<size_Mat)
    {
        ////////////////////////////////////////////////////////////////////////////////
        unsigned int BC;
        BC = geoD[k];

        if (BC != GEO_SOLID && BC != GEO_VOID)
        {
            Distributions27 D27;
            if (EvenOrOdd == true)
            {
                D27.f[dP00] = &DD27[dP00   *size_Mat];
                D27.f[dM00] = &DD27[dM00   *size_Mat];
                D27.f[d0P0] = &DD27[d0P0   *size_Mat];
                D27.f[d0M0] = &DD27[d0M0   *size_Mat];
                D27.f[d00P] = &DD27[d00P   *size_Mat];
                D27.f[d00M] = &DD27[d00M   *size_Mat];
                D27.f[dPP0] = &DD27[dPP0  *size_Mat];
                D27.f[dMM0] = &DD27[dMM0  *size_Mat];
                D27.f[dPM0] = &DD27[dPM0  *size_Mat];
                D27.f[dMP0] = &DD27[dMP0  *size_Mat];
                D27.f[dP0P] = &DD27[dP0P  *size_Mat];
                D27.f[dM0M] = &DD27[dM0M  *size_Mat];
                D27.f[dP0M] = &DD27[dP0M  *size_Mat];
                D27.f[dM0P] = &DD27[dM0P  *size_Mat];
                D27.f[d0PP] = &DD27[d0PP  *size_Mat];
                D27.f[d0MM] = &DD27[d0MM  *size_Mat];
                D27.f[d0PM] = &DD27[d0PM  *size_Mat];
                D27.f[d0MP] = &DD27[d0MP  *size_Mat];
                D27.f[d000] = &DD27[d000*size_Mat];
                D27.f[dPPP] = &DD27[dPPP *size_Mat];
                D27.f[dMMP] = &DD27[dMMP *size_Mat];
                D27.f[dPMP] = &DD27[dPMP *size_Mat];
                D27.f[dMPP] = &DD27[dMPP *size_Mat];
                D27.f[dPPM] = &DD27[dPPM *size_Mat];
                D27.f[dMMM] = &DD27[dMMM *size_Mat];
                D27.f[dPMM]= &DD27[dPMM *size_Mat];
                D27.f[dMPM]= &DD27[dMPM *size_Mat];
            }
            else
            {
                D27.f[dM00] = &DD27[dP00   *size_Mat];
                D27.f[dP00] = &DD27[dM00   *size_Mat];
                D27.f[d0M0] = &DD27[d0P0   *size_Mat];
                D27.f[d0P0] = &DD27[d0M0   *size_Mat];
                D27.f[d00M] = &DD27[d00P   *size_Mat];
                D27.f[d00P] = &DD27[d00M   *size_Mat];
                D27.f[dMM0] = &DD27[dPP0  *size_Mat];
                D27.f[dPP0] = &DD27[dMM0  *size_Mat];
                D27.f[dMP0] = &DD27[dPM0  *size_Mat];
                D27.f[dPM0] = &DD27[dMP0  *size_Mat];
                D27.f[dM0M] = &DD27[dP0P  *size_Mat];
                D27.f[dP0P] = &DD27[dM0M  *size_Mat];
                D27.f[dM0P] = &DD27[dP0M  *size_Mat];
                D27.f[dP0M] = &DD27[dM0P  *size_Mat];
                D27.f[d0MM] = &DD27[d0PP  *size_Mat];
                D27.f[d0PP] = &DD27[d0MM  *size_Mat];
                D27.f[d0MP] = &DD27[d0PM  *size_Mat];
                D27.f[d0PM] = &DD27[d0MP  *size_Mat];
                D27.f[d000] = &DD27[d000*size_Mat];
                D27.f[dMMM] = &DD27[dPPP *size_Mat];
                D27.f[dPPM] = &DD27[dMMP *size_Mat];
                D27.f[dMPM]= &DD27[dPMP *size_Mat];
                D27.f[dPMM]= &DD27[dMPP *size_Mat];
                D27.f[dMMP] = &DD27[dPPM *size_Mat];
                D27.f[dPPP] = &DD27[dMMM *size_Mat];
                D27.f[dMPP] = &DD27[dPMM *size_Mat];
                D27.f[dPMP] = &DD27[dMPM *size_Mat];
            }
            //////////////////////////////////////////////////////////////////////////
            real ConcD = Conc[k];
            real   vx1 = ux[k];
            real   vx2 = uy[k];
            real   vx3 = uz[k];
            //real lambdaD     = -three + sqrt(three);
            //real Diffusivity = c1o20;
            //real Lam         = -(c1o2+one/lambdaD);
            //real nue_d       = Lam/three;
            //real ae          = Diffusivity/nue_d - one;
            //real ux_sq       = vx1 * vx1;
            //real uy_sq       = vx2 * vx2;
            //real uz_sq       = vx3 * vx3;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //D3Q7
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //index
            //unsigned int kzero= k;
            //unsigned int ke   = k;
            //unsigned int kw   = neighborX[k];
            //unsigned int kn   = k;
            //unsigned int ks   = neighborY[k];
            //unsigned int kt   = k;
            //unsigned int kb   = neighborZ[k];
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //(D7.f[0])[kzero] = ConcD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
            //(D7.f[1])[ke   ] = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)+vx1*c1o2);
            //(D7.f[2])[kw   ] = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)-vx1*c1o2);
            //(D7.f[3])[kn   ] = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)+vx2*c1o2);
            //(D7.f[4])[ks   ] = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)-vx2*c1o2);
            //(D7.f[5])[kt   ] = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)+vx3*c1o2);
            //(D7.f[6])[kb   ] = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)-vx3*c1o2);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //D3Q27
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //index
            unsigned int kzero = k;
            unsigned int ke = k;
            unsigned int kw = neighborX[k];
            unsigned int kn = k;
            unsigned int ks = neighborY[k];
            unsigned int kt = k;
            unsigned int kb = neighborZ[k];
            unsigned int ksw = neighborY[kw];
            unsigned int kne = k;
            unsigned int kse = ks;
            unsigned int knw = kw;
            unsigned int kbw = neighborZ[kw];
            unsigned int kte = k;
            unsigned int kbe = kb;
            unsigned int ktw = kw;
            unsigned int kbs = neighborZ[ks];
            unsigned int ktn = k;
            unsigned int kbn = kb;
            unsigned int kts = ks;
            unsigned int ktse = ks;
            unsigned int kbnw = kbw;
            unsigned int ktnw = kw;
            unsigned int kbse = kbs;
            unsigned int ktsw = ksw;
            unsigned int kbne = kb;
            unsigned int ktne = k;
            unsigned int kbsw = neighborZ[ksw];
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            real cu_sq = c3o2*(vx1*vx1 + vx2*vx2 + vx3*vx3);

            (D27.f[d000])[kzero] = c8o27* ConcD*(c1o1 - cu_sq);
            (D27.f[dP00])[ke] = c2o27* ConcD*(c1o1 + c3o1*(vx1)+c9o2*(vx1)*(vx1)-cu_sq);
            (D27.f[dM00])[kw] = c2o27* ConcD*(c1o1 + c3o1*(-vx1) + c9o2*(-vx1)*(-vx1) - cu_sq);
            (D27.f[d0P0])[kn] = c2o27* ConcD*(c1o1 + c3o1*(vx2)+c9o2*(vx2)*(vx2)-cu_sq);
            (D27.f[d0M0])[ks] = c2o27* ConcD*(c1o1 + c3o1*(-vx2) + c9o2*(-vx2)*(-vx2) - cu_sq);
            (D27.f[d00P])[kt] = c2o27* ConcD*(c1o1 + c3o1*(vx3)+c9o2*(vx3)*(vx3)-cu_sq);
            (D27.f[d00M])[kb] = c2o27* ConcD*(c1o1 + c3o1*(-vx3) + c9o2*(-vx3)*(-vx3) - cu_sq);
            (D27.f[dPP0])[kne] = c1o54* ConcD*(c1o1 + c3o1*(vx1 + vx2) + c9o2*(vx1 + vx2)*(vx1 + vx2) - cu_sq);
            (D27.f[dMM0])[ksw] = c1o54* ConcD*(c1o1 + c3o1*(-vx1 - vx2) + c9o2*(-vx1 - vx2)*(-vx1 - vx2) - cu_sq);
            (D27.f[dPM0])[kse] = c1o54* ConcD*(c1o1 + c3o1*(vx1 - vx2) + c9o2*(vx1 - vx2)*(vx1 - vx2) - cu_sq);
            (D27.f[dMP0])[knw] = c1o54* ConcD*(c1o1 + c3o1*(-vx1 + vx2) + c9o2*(-vx1 + vx2)*(-vx1 + vx2) - cu_sq);
            (D27.f[dP0P])[kte] = c1o54* ConcD*(c1o1 + c3o1*(vx1 + vx3) + c9o2*(vx1 + vx3)*(vx1 + vx3) - cu_sq);
            (D27.f[dM0M])[kbw] = c1o54* ConcD*(c1o1 + c3o1*(-vx1 - vx3) + c9o2*(-vx1 - vx3)*(-vx1 - vx3) - cu_sq);
            (D27.f[dP0M])[kbe] = c1o54* ConcD*(c1o1 + c3o1*(vx1 - vx3) + c9o2*(vx1 - vx3)*(vx1 - vx3) - cu_sq);
            (D27.f[dM0P])[ktw] = c1o54* ConcD*(c1o1 + c3o1*(-vx1 + vx3) + c9o2*(-vx1 + vx3)*(-vx1 + vx3) - cu_sq);
            (D27.f[d0PP])[ktn] = c1o54* ConcD*(c1o1 + c3o1*(vx2 + vx3) + c9o2*(vx2 + vx3)*(vx2 + vx3) - cu_sq);
            (D27.f[d0MM])[kbs] = c1o54* ConcD*(c1o1 + c3o1*(-vx2 - vx3) + c9o2*(-vx2 - vx3)*(-vx2 - vx3) - cu_sq);
            (D27.f[d0PM])[kbn] = c1o54* ConcD*(c1o1 + c3o1*(vx2 - vx3) + c9o2*(vx2 - vx3)*(vx2 - vx3) - cu_sq);
            (D27.f[d0MP])[kts] = c1o54* ConcD*(c1o1 + c3o1*(-vx2 + vx3) + c9o2*(-vx2 + vx3)*(-vx2 + vx3) - cu_sq);
            (D27.f[dPPP])[ktne] = c1o216*ConcD*(c1o1 + c3o1*(vx1 + vx2 + vx3) + c9o2*(vx1 + vx2 + vx3)*(vx1 + vx2 + vx3) - cu_sq);
            (D27.f[dMMM])[kbsw] = c1o216*ConcD*(c1o1 + c3o1*(-vx1 - vx2 - vx3) + c9o2*(-vx1 - vx2 - vx3)*(-vx1 - vx2 - vx3) - cu_sq);
            (D27.f[dPPM])[kbne] = c1o216*ConcD*(c1o1 + c3o1*(vx1 + vx2 - vx3) + c9o2*(vx1 + vx2 - vx3)*(vx1 + vx2 - vx3) - cu_sq);
            (D27.f[dMMP])[ktsw] = c1o216*ConcD*(c1o1 + c3o1*(-vx1 - vx2 + vx3) + c9o2*(-vx1 - vx2 + vx3)*(-vx1 - vx2 + vx3) - cu_sq);
            (D27.f[dPMP])[ktse] = c1o216*ConcD*(c1o1 + c3o1*(vx1 - vx2 + vx3) + c9o2*(vx1 - vx2 + vx3)*(vx1 - vx2 + vx3) - cu_sq);
            (D27.f[dMPM])[kbnw] = c1o216*ConcD*(c1o1 + c3o1*(-vx1 + vx2 - vx3) + c9o2*(-vx1 + vx2 - vx3)*(-vx1 + vx2 - vx3) - cu_sq);
            (D27.f[dPMM])[kbse] = c1o216*ConcD*(c1o1 + c3o1*(vx1 - vx2 - vx3) + c9o2*(vx1 - vx2 - vx3)*(vx1 - vx2 - vx3) - cu_sq);
            (D27.f[dMPP])[ktnw] = c1o216*ConcD*(c1o1 + c3o1*(-vx1 + vx2 + vx3) + c9o2*(-vx1 + vx2 + vx3)*(-vx1 + vx2 + vx3) - cu_sq);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }
    }
}
