#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;
#include "math.h"

__global__ void LB_Kernel_PM_Cum_One_Comp_SP_27(real omega,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	unsigned long long numberOfLBnodes,
	int level,
	real* forces,
	real porosity,
	real darcy,
	real forchheimer,
	unsigned int sizeOfPorousMedia,
	unsigned int* nodeIdsPorousMedia,
	bool EvenOrOdd)
{
	Distributions27 D;
	if (EvenOrOdd == true)
	{
		D.f[DIR_P00] = &DDStart[DIR_P00   *numberOfLBnodes];
		D.f[DIR_M00] = &DDStart[DIR_M00   *numberOfLBnodes];
		D.f[DIR_0P0] = &DDStart[DIR_0P0   *numberOfLBnodes];
		D.f[DIR_0M0] = &DDStart[DIR_0M0   *numberOfLBnodes];
		D.f[DIR_00P] = &DDStart[DIR_00P   *numberOfLBnodes];
		D.f[DIR_00M] = &DDStart[DIR_00M   *numberOfLBnodes];
		D.f[DIR_PP0] = &DDStart[DIR_PP0  *numberOfLBnodes];
		D.f[DIR_MM0] = &DDStart[DIR_MM0  *numberOfLBnodes];
		D.f[DIR_PM0] = &DDStart[DIR_PM0  *numberOfLBnodes];
		D.f[DIR_MP0] = &DDStart[DIR_MP0  *numberOfLBnodes];
		D.f[DIR_P0P] = &DDStart[DIR_P0P  *numberOfLBnodes];
		D.f[DIR_M0M] = &DDStart[DIR_M0M  *numberOfLBnodes];
		D.f[DIR_P0M] = &DDStart[DIR_P0M  *numberOfLBnodes];
		D.f[DIR_M0P] = &DDStart[DIR_M0P  *numberOfLBnodes];
		D.f[DIR_0PP] = &DDStart[DIR_0PP  *numberOfLBnodes];
		D.f[DIR_0MM] = &DDStart[DIR_0MM  *numberOfLBnodes];
		D.f[DIR_0PM] = &DDStart[DIR_0PM  *numberOfLBnodes];
		D.f[DIR_0MP] = &DDStart[DIR_0MP  *numberOfLBnodes];
		D.f[DIR_000] = &DDStart[DIR_000*numberOfLBnodes];
		D.f[DIR_PPP] = &DDStart[DIR_PPP *numberOfLBnodes];
		D.f[DIR_MMP] = &DDStart[DIR_MMP *numberOfLBnodes];
		D.f[DIR_PMP] = &DDStart[DIR_PMP *numberOfLBnodes];
		D.f[DIR_MPP] = &DDStart[DIR_MPP *numberOfLBnodes];
		D.f[DIR_PPM] = &DDStart[DIR_PPM *numberOfLBnodes];
		D.f[DIR_MMM] = &DDStart[DIR_MMM *numberOfLBnodes];
		D.f[DIR_PMM] = &DDStart[DIR_PMM *numberOfLBnodes];
		D.f[DIR_MPM] = &DDStart[DIR_MPM *numberOfLBnodes];
	}
	else
	{
		D.f[DIR_M00] = &DDStart[DIR_P00   *numberOfLBnodes];
		D.f[DIR_P00] = &DDStart[DIR_M00   *numberOfLBnodes];
		D.f[DIR_0M0] = &DDStart[DIR_0P0   *numberOfLBnodes];
		D.f[DIR_0P0] = &DDStart[DIR_0M0   *numberOfLBnodes];
		D.f[DIR_00M] = &DDStart[DIR_00P   *numberOfLBnodes];
		D.f[DIR_00P] = &DDStart[DIR_00M   *numberOfLBnodes];
		D.f[DIR_MM0] = &DDStart[DIR_PP0  *numberOfLBnodes];
		D.f[DIR_PP0] = &DDStart[DIR_MM0  *numberOfLBnodes];
		D.f[DIR_MP0] = &DDStart[DIR_PM0  *numberOfLBnodes];
		D.f[DIR_PM0] = &DDStart[DIR_MP0  *numberOfLBnodes];
		D.f[DIR_M0M] = &DDStart[DIR_P0P  *numberOfLBnodes];
		D.f[DIR_P0P] = &DDStart[DIR_M0M  *numberOfLBnodes];
		D.f[DIR_M0P] = &DDStart[DIR_P0M  *numberOfLBnodes];
		D.f[DIR_P0M] = &DDStart[DIR_M0P  *numberOfLBnodes];
		D.f[DIR_0MM] = &DDStart[DIR_0PP  *numberOfLBnodes];
		D.f[DIR_0PP] = &DDStart[DIR_0MM  *numberOfLBnodes];
		D.f[DIR_0MP] = &DDStart[DIR_0PM  *numberOfLBnodes];
		D.f[DIR_0PM] = &DDStart[DIR_0MP  *numberOfLBnodes];
		D.f[DIR_000] = &DDStart[DIR_000*numberOfLBnodes];
		D.f[DIR_MMM] = &DDStart[DIR_PPP *numberOfLBnodes];
		D.f[DIR_PPM] = &DDStart[DIR_MMP *numberOfLBnodes];
		D.f[DIR_MPM] = &DDStart[DIR_PMP *numberOfLBnodes];
		D.f[DIR_PMM] = &DDStart[DIR_MPP *numberOfLBnodes];
		D.f[DIR_MMP] = &DDStart[DIR_PPM *numberOfLBnodes];
		D.f[DIR_PPP] = &DDStart[DIR_MMM *numberOfLBnodes];
		D.f[DIR_MPP] = &DDStart[DIR_PMM *numberOfLBnodes];
		D.f[DIR_PMP] = &DDStart[DIR_MPM *numberOfLBnodes];
	}

	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned kt = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if (kt<sizeOfPorousMedia)
	{
		////////////////////////////////////////////////////////////////////////////////
		//index
		unsigned int k = nodeIdsPorousMedia[kt];
		unsigned int kw = neighborX[k];
		unsigned int ks = neighborY[k];
		unsigned int kb = neighborZ[k];
		unsigned int ksw = neighborY[kw];
		unsigned int kbw = neighborZ[kw];
		unsigned int kbs = neighborZ[ks];
		unsigned int kbsw = neighborZ[ksw];
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		real mfcbb = (D.f[DIR_P00])[k];
		real mfabb = (D.f[DIR_M00])[kw];
		real mfbcb = (D.f[DIR_0P0])[k];
		real mfbab = (D.f[DIR_0M0])[ks];
		real mfbbc = (D.f[DIR_00P])[k];
		real mfbba = (D.f[DIR_00M])[kb];
		real mfccb = (D.f[DIR_PP0])[k];
		real mfaab = (D.f[DIR_MM0])[ksw];
		real mfcab = (D.f[DIR_PM0])[ks];
		real mfacb = (D.f[DIR_MP0])[kw];
		real mfcbc = (D.f[DIR_P0P])[k];
		real mfaba = (D.f[DIR_M0M])[kbw];
		real mfcba = (D.f[DIR_P0M])[kb];
		real mfabc = (D.f[DIR_M0P])[kw];
		real mfbcc = (D.f[DIR_0PP])[k];
		real mfbaa = (D.f[DIR_0MM])[kbs];
		real mfbca = (D.f[DIR_0PM])[kb];
		real mfbac = (D.f[DIR_0MP])[ks];
		real mfbbb = (D.f[DIR_000])[k];
		real mfccc = (D.f[DIR_PPP])[k];
		real mfaac = (D.f[DIR_MMP])[ksw];
		real mfcac = (D.f[DIR_PMP])[ks];
		real mfacc = (D.f[DIR_MPP])[kw];
		real mfcca = (D.f[DIR_PPM])[kb];
		real mfaaa = (D.f[DIR_MMM])[kbsw];
		real mfcaa = (D.f[DIR_PMM])[kbs];
		real mfaca = (D.f[DIR_MPM])[kbw];
		////////////////////////////////////////////////////////////////////////////////////
		real drho = ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
			(((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
			((mfabb + mfcbb) + (mfbab + mfbcb)) + (mfbba + mfbbc)) + mfbbb;

		real rho = c1o1 + drho;
		////////////////////////////////////////////////////////////////////////////////////
		real vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
			(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
			(mfcbb - mfabb)) / rho;
		real vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
			(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
			(mfbcb - mfbab)) / rho;
		real vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
			(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
			(mfbbc - mfbba)) / rho;
		real vx2 = vvx*vvx;
		real vy2 = vvy*vvy;
		real vz2 = vvz*vvz;
		////////////////////////////////////////////////////////////////////////////////////
		//porous media
		vvx = -(c2o1 * vvx) / (-c2o1 - darcy - forchheimer * sqrtf(vx2 + vy2 + vz2));
		vvy = -(c2o1 * vvy) / (-c2o1 - darcy - forchheimer * sqrtf(vx2 + vy2 + vz2));
		vvz = -(c2o1 * vvz) / (-c2o1 - darcy - forchheimer * sqrtf(vx2 + vy2 + vz2));
		//vvx = (two * vvx) / (two + 134.4 + 0.0068287 * sqrtf(vx2 + vy2 + vz2));
		//vvy = (two * vvy) / (two + 134.4 + 0.0068287 * sqrtf(vx2 + vy2 + vz2));
		//vvz = (two * vvz) / (two + 134.4 + 0.0068287 * sqrtf(vx2 + vy2 + vz2));
		////////////////////////////////////////////////////////////////////////////////////
		//the force be with you
		real fx = forces[0] / (pow((double)c2o1, (double)level)); //zero;//0.0032653/(pow(two,level)); //0.000000005;//(two/1600000.0) / 120.0; //
		real fy = forces[1] / (pow((double)c2o1, (double)level)); //zero;
		real fz = forces[2] / (pow((double)c2o1, (double)level)); //zero;
		vvx += fx*c1o2;
		vvy += fy*c1o2;
		vvz += fz*c1o2;
		////////////////////////////////////////////////////////////////////////////////////
		//real omega = omega_in;
		////////////////////////////////////////////////////////////////////////////////////
		real oMdrho = c1o1; // comp special
						   ////////////////////////////////////////////////////////////////////////////////////
		real m0, m1, m2;
		//////////////////////////////////////////////////////////////////////////////////////
		real wadjust;
		real qudricLimitP = 0.01f;// * 0.0001f;
		real qudricLimitM = 0.01f;// * 0.0001f;
		real qudricLimitD = 0.01f;// * 0.001f;
								  ////////////////////////////////////////////////////////////////////////////////////
								  //Hin
								  ////////////////////////////////////////////////////////////////////////////////////
								  // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
								  ////////////////////////////////////////////////////////////////////////////////////
								  // Z - Dir
		m2 = mfaaa + mfaac;
		m1 = mfaac - mfaaa;
		m0 = m2 + mfaab;
		mfaaa = m0;
		m0 += c1o36 * oMdrho;
		mfaab = m1 - m0 * vvz;
		mfaac = m2 - c2o1*	m1 * vvz + vz2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfaba + mfabc;
		m1 = mfabc - mfaba;
		m0 = m2 + mfabb;
		mfaba = m0;
		m0 += c1o9 * oMdrho;
		mfabb = m1 - m0 * vvz;
		mfabc = m2 - c2o1*	m1 * vvz + vz2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfaca + mfacc;
		m1 = mfacc - mfaca;
		m0 = m2 + mfacb;
		mfaca = m0;
		m0 += c1o36 * oMdrho;
		mfacb = m1 - m0 * vvz;
		mfacc = m2 - c2o1*	m1 * vvz + vz2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfbaa + mfbac;
		m1 = mfbac - mfbaa;
		m0 = m2 + mfbab;
		mfbaa = m0;
		m0 += c1o9 * oMdrho;
		mfbab = m1 - m0 * vvz;
		mfbac = m2 - c2o1*	m1 * vvz + vz2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfbba + mfbbc;
		m1 = mfbbc - mfbba;
		m0 = m2 + mfbbb;
		mfbba = m0;
		m0 += c4o9 * oMdrho;
		mfbbb = m1 - m0 * vvz;
		mfbbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfbca + mfbcc;
		m1 = mfbcc - mfbca;
		m0 = m2 + mfbcb;
		mfbca = m0;
		m0 += c1o9 * oMdrho;
		mfbcb = m1 - m0 * vvz;
		mfbcc = m2 - c2o1*	m1 * vvz + vz2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfcaa + mfcac;
		m1 = mfcac - mfcaa;
		m0 = m2 + mfcab;
		mfcaa = m0;
		m0 += c1o36 * oMdrho;
		mfcab = m1 - m0 * vvz;
		mfcac = m2 - c2o1*	m1 * vvz + vz2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfcba + mfcbc;
		m1 = mfcbc - mfcba;
		m0 = m2 + mfcbb;
		mfcba = m0;
		m0 += c1o9 * oMdrho;
		mfcbb = m1 - m0 * vvz;
		mfcbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfcca + mfccc;
		m1 = mfccc - mfcca;
		m0 = m2 + mfccb;
		mfcca = m0;
		m0 += c1o36 * oMdrho;
		mfccb = m1 - m0 * vvz;
		mfccc = m2 - c2o1*	m1 * vvz + vz2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
		////////////////////////////////////////////////////////////////////////////////////
		// Y - Dir
		m2 = mfaaa + mfaca;
		m1 = mfaca - mfaaa;
		m0 = m2 + mfaba;
		mfaaa = m0;
		m0 += c1o6 * oMdrho;
		mfaba = m1 - m0 * vvy;
		mfaca = m2 - c2o1*	m1 * vvy + vy2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfaab + mfacb;
		m1 = mfacb - mfaab;
		m0 = m2 + mfabb;
		mfaab = m0;
		mfabb = m1 - m0 * vvy;
		mfacb = m2 - c2o1*	m1 * vvy + vy2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfaac + mfacc;
		m1 = mfacc - mfaac;
		m0 = m2 + mfabc;
		mfaac = m0;
		m0 += c1o18 * oMdrho;
		mfabc = m1 - m0 * vvy;
		mfacc = m2 - c2o1*	m1 * vvy + vy2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfbaa + mfbca;
		m1 = mfbca - mfbaa;
		m0 = m2 + mfbba;
		mfbaa = m0;
		m0 += c2o3 * oMdrho;
		mfbba = m1 - m0 * vvy;
		mfbca = m2 - c2o1*	m1 * vvy + vy2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfbab + mfbcb;
		m1 = mfbcb - mfbab;
		m0 = m2 + mfbbb;
		mfbab = m0;
		mfbbb = m1 - m0 * vvy;
		mfbcb = m2 - c2o1*	m1 * vvy + vy2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfbac + mfbcc;
		m1 = mfbcc - mfbac;
		m0 = m2 + mfbbc;
		mfbac = m0;
		m0 += c2o9 * oMdrho;
		mfbbc = m1 - m0 * vvy;
		mfbcc = m2 - c2o1*	m1 * vvy + vy2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfcaa + mfcca;
		m1 = mfcca - mfcaa;
		m0 = m2 + mfcba;
		mfcaa = m0;
		m0 += c1o6 * oMdrho;
		mfcba = m1 - m0 * vvy;
		mfcca = m2 - c2o1*	m1 * vvy + vy2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfcab + mfccb;
		m1 = mfccb - mfcab;
		m0 = m2 + mfcbb;
		mfcab = m0;
		mfcbb = m1 - m0 * vvy;
		mfccb = m2 - c2o1*	m1 * vvy + vy2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfcac + mfccc;
		m1 = mfccc - mfcac;
		m0 = m2 + mfcbc;
		mfcac = m0;
		m0 += c1o18 * oMdrho;
		mfcbc = m1 - m0 * vvy;
		mfccc = m2 - c2o1*	m1 * vvy + vy2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
		////////////////////////////////////////////////////////////////////////////////////
		// X - Dir
		m2 = mfaaa + mfcaa;
		m1 = mfcaa - mfaaa;
		m0 = m2 + mfbaa;
		mfaaa = m0;
		m0 += c1o1* oMdrho;
		mfbaa = m1 - m0 * vvx;
		mfcaa = m2 - c2o1*	m1 * vvx + vx2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfaba + mfcba;
		m1 = mfcba - mfaba;
		m0 = m2 + mfbba;
		mfaba = m0;
		mfbba = m1 - m0 * vvx;
		mfcba = m2 - c2o1*	m1 * vvx + vx2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfaca + mfcca;
		m1 = mfcca - mfaca;
		m0 = m2 + mfbca;
		mfaca = m0;
		m0 += c1o3 * oMdrho;
		mfbca = m1 - m0 * vvx;
		mfcca = m2 - c2o1*	m1 * vvx + vx2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfaab + mfcab;
		m1 = mfcab - mfaab;
		m0 = m2 + mfbab;
		mfaab = m0;
		mfbab = m1 - m0 * vvx;
		mfcab = m2 - c2o1*	m1 * vvx + vx2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfabb + mfcbb;
		m1 = mfcbb - mfabb;
		m0 = m2 + mfbbb;
		mfabb = m0;
		mfbbb = m1 - m0 * vvx;
		mfcbb = m2 - c2o1*	m1 * vvx + vx2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfacb + mfccb;
		m1 = mfccb - mfacb;
		m0 = m2 + mfbcb;
		mfacb = m0;
		mfbcb = m1 - m0 * vvx;
		mfccb = m2 - c2o1*	m1 * vvx + vx2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfaac + mfcac;
		m1 = mfcac - mfaac;
		m0 = m2 + mfbac;
		mfaac = m0;
		m0 += c1o3 * oMdrho;
		mfbac = m1 - m0 * vvx;
		mfcac = m2 - c2o1*	m1 * vvx + vx2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfabc + mfcbc;
		m1 = mfcbc - mfabc;
		m0 = m2 + mfbbc;
		mfabc = m0;
		mfbbc = m1 - m0 * vvx;
		mfcbc = m2 - c2o1*	m1 * vvx + vx2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		m2 = mfacc + mfccc;
		m1 = mfccc - mfacc;
		m0 = m2 + mfbcc;
		mfacc = m0;
		m0 += c1o9 * oMdrho;
		mfbcc = m1 - m0 * vvx;
		mfccc = m2 - c2o1*	m1 * vvx + vx2 * m0;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////


		////////////////////////////////////////////////////////////////////////////////////
		// Cumulants
		////////////////////////////////////////////////////////////////////////////////////
		real OxxPyyPzz = c1o1;	//set the bulk viscosity one is high / two is very low and zero is (too) high

								////////////////////////////////////////////////////////////
								//3.
								//////////////////////////////
		real OxyyPxzz = c1o1;
		real OxyyMxzz = c1o1;
		real Oxyz = c1o1;
		////////////////////////////////////////////////////////////
		//4.
		//////////////////////////////
		real O4 = c1o1;
		////////////////////////////////////////////////////////////
		//5.
		//////////////////////////////
		real O5 = c1o1;
		////////////////////////////////////////////////////////////
		//6.
		//////////////////////////////
		real O6 = c1o1;
		////////////////////////////////////////////////////////////


		//central moments to cumulants
		//4.
		real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho;
		real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho;
		real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho;

		real CUMcca = mfcca - (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho - c1o9*(drho / rho));
		real CUMcac = mfcac - (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho - c1o9*(drho / rho));
		real CUMacc = mfacc - (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho - c1o9*(drho / rho));

		//5.
		real CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac)) / rho;
		real CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc)) / rho;
		real CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab)) / rho;

		//6.

		real CUMccc = mfccc + ((-c4o1 *  mfbbb * mfbbb
			- (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
			- c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
			- c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
			+ (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
				+ c2o1 * (mfcaa * mfaca * mfaac)
				+ c16o1 *  mfbba * mfbab * mfabb) / (rho * rho)
			- c1o3 * (mfacc + mfcac + mfcca) / rho
			- c1o9 * (mfcaa + mfaca + mfaac) / rho
			+ (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
				+ (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3
			+ c1o27*((drho * drho - drho) / (rho*rho)));
		//+ c1o27*(one -three/rho +two/(rho*rho)));




		//2.
		// linear combinations
		real mxxPyyPzz = mfcaa + mfaca + mfaac;
		real mxxMyy = mfcaa - mfaca;
		real mxxMzz = mfcaa - mfaac;

		//////////////////////////////////////////////////////////////////////////
		//real magicBulk = (CUMacc + CUMcac + CUMcca)*(c1o1 / OxxPyyPzz - c1o2)*c3o2*8.;

		//////////////////////////////////////////////////////////////////////////
		//limiter-Scheise Teil 1
		//real oxxyy,oxxzz,oxy,oxz,oyz;
		//real smag=0.001;
		//oxxyy    = omega+(one-omega)*abs(mxxMyy)/(abs(mxxMyy)+smag);
		//oxxzz    = omega+(one-omega)*abs(mxxMzz)/(abs(mxxMzz)+smag);
		//oxy      = omega+(one-omega)*abs(mfbba)/(abs(mfbba)+smag);
		//oxz      = omega+(one-omega)*abs(mfbab)/(abs(mfbab)+smag);
		//oyz      = omega+(one-omega)*abs(mfabb)/(abs(mfabb)+smag);

		////////////////////////////////////////////////////////////////////////////
		////Teil 1b
		//real constante = 1000.0;
		//real nuEddi = constante * abs(mxxPyyPzz);
		//real omegaLimit = one / (one / omega + three * nuEddi);

		//{
		//	real dxux = c1o2 * (-omegaLimit) *(mxxMyy + mxxMzz) +  OxxPyyPzz * (mfaaa - mxxPyyPzz);
		//	real dyuy = dxux + omegaLimit * c3o2 * mxxMyy;
		//	real dzuz = dxux + omegaLimit * c3o2 * mxxMzz;

		////relax
		//mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
		//mxxMyy    += omegaLimit * (-mxxMyy) - three * (one + c1o2 * (-omegaLimit)) * (vx2 * dxux + vy2 * dyuy);
		//mxxMzz    += omegaLimit * (-mxxMzz) - three * (one + c1o2 * (-omegaLimit)) * (vx2 * dxux + vz2 * dzuz);

		//}
		//mfabb     += omegaLimit * (-mfabb);
		//mfbab     += omegaLimit * (-mfbab);
		//mfbba     += omegaLimit * (-mfbba);
		////////////////////////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)
		{
			real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
			real dyuy = dxux + omega * c3o2 * mxxMyy;
			real dzuz = dxux + omega * c3o2 * mxxMzz;

			////relax original
			//mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
			//mxxMyy    += omega * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
			//mxxMzz    += omega * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);
			//relax porous media Mp
			mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz + rho * (vx2 + vy2 + vz2) * (c1o1 / porosity - c1o1)) - c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
																																													  //relax porous media Ms
			mxxMyy += omega * (rho * (vx2 - vy2) * (c1o1 / porosity - c1o1) - mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
			mxxMzz += omega * (rho * (vx2 - vz2) * (c1o1 / porosity - c1o1) - mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

			//////////////////////////////////////////////////////////////////////////
			//limiter-Scheise Teil 2
			//mxxMyy    += oxxyy * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vy2 * dyuy);
			//mxxMzz    += oxxzz * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vz2 * dzuz);
			//////////////////////////////////////////////////////////////////////////

		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////no correction
		//mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);//-magicBulk*OxxPyyPzz;
		//mxxMyy    += -(-omega) * (-mxxMyy);
		//mxxMzz    += -(-omega) * (-mxxMzz);
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////original
		//mfabb     += omega * (-mfabb);
		//mfbab     += omega * (-mfbab);
		//mfbba     += omega * (-mfbba);
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//porous media M11
		mfabb += omega * (rho * vvy * vvz * (c1o1 / porosity - c1o1) - mfabb);
		mfbab += omega * (rho * vvx * vvz * (c1o1 / porosity - c1o1) - mfbab);
		mfbba += omega * (rho * vvx * vvy * (c1o1 / porosity - c1o1) - mfbba);

		//////////////////////////////////////////////////////////////////////////
		//limiter-Scheise Teil 3
		//mfabb     += oyz * (-mfabb);
		//mfbab     += oxz * (-mfbab);
		//mfbba     += oxy * (-mfbba);
		//////////////////////////////////////////////////////////////////////////

		// linear combinations back
		mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
		mfaca = c1o3 * (-c2o1*  mxxMyy + mxxMzz + mxxPyyPzz);
		mfaac = c1o3 * (mxxMyy - c2o1* mxxMzz + mxxPyyPzz);

		//3.
		// linear combinations

		real mxxyPyzz = mfcba + mfabc;
		real mxxyMyzz = mfcba - mfabc;

		real mxxzPyyz = mfcab + mfacb;
		real mxxzMyyz = mfcab - mfacb;

		real mxyyPxzz = mfbca + mfbac;
		real mxyyMxzz = mfbca - mfbac;

		//relax
		//////////////////////////////////////////////////////////////////////////
		//das ist der limiter
		wadjust = Oxyz + (c1o1 - Oxyz)*abs(mfbbb) / (abs(mfbbb) + qudricLimitD);
		mfbbb += wadjust * (-mfbbb);
		wadjust = OxyyPxzz + (c1o1 - OxyyPxzz)*abs(mxxyPyzz) / (abs(mxxyPyzz) + qudricLimitP);
		mxxyPyzz += wadjust * (-mxxyPyzz);
		wadjust = OxyyMxzz + (c1o1 - OxyyMxzz)*abs(mxxyMyzz) / (abs(mxxyMyzz) + qudricLimitM);
		mxxyMyzz += wadjust * (-mxxyMyzz);
		wadjust = OxyyPxzz + (c1o1 - OxyyPxzz)*abs(mxxzPyyz) / (abs(mxxzPyyz) + qudricLimitP);
		mxxzPyyz += wadjust * (-mxxzPyyz);
		wadjust = OxyyMxzz + (c1o1 - OxyyMxzz)*abs(mxxzMyyz) / (abs(mxxzMyyz) + qudricLimitM);
		mxxzMyyz += wadjust * (-mxxzMyyz);
		wadjust = OxyyPxzz + (c1o1 - OxyyPxzz)*abs(mxyyPxzz) / (abs(mxyyPxzz) + qudricLimitP);
		mxyyPxzz += wadjust * (-mxyyPxzz);
		wadjust = OxyyMxzz + (c1o1 - OxyyMxzz)*abs(mxyyMxzz) / (abs(mxyyMxzz) + qudricLimitM);
		mxyyMxzz += wadjust * (-mxyyMxzz);
		//////////////////////////////////////////////////////////////////////////
		//ohne limiter
		//mfbbb     += OxyyMxzz * (-mfbbb);
		//mxxyPyzz  += OxyyPxzz * (-mxxyPyzz);
		//mxxyMyzz  += OxyyMxzz * (-mxxyMyzz);
		//mxxzPyyz  += OxyyPxzz * (-mxxzPyyz);
		//mxxzMyyz  += OxyyMxzz * (-mxxzMyyz);
		//mxyyPxzz  += OxyyPxzz * (-mxyyPxzz);
		//mxyyMxzz  += OxyyMxzz * (-mxyyMxzz);
		//////////////////////////////////////////////////////////////////////////

		//// linear combinations back

		mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
		mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
		mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
		mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
		mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
		mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

		//4.
		//////////////////////////////////////////////////////////////////////////
		//mit limiter
		//	wadjust    = O4+(one-O4)*abs(CUMacc)/(abs(CUMacc)+qudricLimit);
		//CUMacc    += wadjust * (-CUMacc);
		//	wadjust    = O4+(one-O4)*abs(CUMcac)/(abs(CUMcac)+qudricLimit);
		//CUMcac    += wadjust * (-CUMcac); 
		//	wadjust    = O4+(one-O4)*abs(CUMcca)/(abs(CUMcca)+qudricLimit);
		//CUMcca    += wadjust * (-CUMcca); 

		//	wadjust    = O4+(one-O4)*abs(CUMbbc)/(abs(CUMbbc)+qudricLimit);
		//CUMbbc    += wadjust * (-CUMbbc); 
		//	wadjust    = O4+(one-O4)*abs(CUMbcb)/(abs(CUMbcb)+qudricLimit);
		//CUMbcb    += wadjust * (-CUMbcb); 
		//	wadjust    = O4+(one-O4)*abs(CUMcbb)/(abs(CUMcbb)+qudricLimit);
		//CUMcbb    += wadjust * (-CUMcbb); 
		//////////////////////////////////////////////////////////////////////////
		//ohne limiter
		CUMacc += O4 * (-CUMacc);
		CUMcac += O4 * (-CUMcac);
		CUMcca += O4 * (-CUMcca);

		CUMbbc += O4 * (-CUMbbc);
		CUMbcb += O4 * (-CUMbcb);
		CUMcbb += O4 * (-CUMcbb);
		//////////////////////////////////////////////////////////////////////////


		//5.
		CUMbcc += O5 * (-CUMbcc);
		CUMcbc += O5 * (-CUMcbc);
		CUMccb += O5 * (-CUMccb);

		//6.
		CUMccc += O6 * (-CUMccc);



		//back cumulants to central moments
		//4.
		mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho;
		mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho;
		mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho;

		mfcca = CUMcca + (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho - c1o9*(drho / rho));
		mfcac = CUMcac + (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho - c1o9*(drho / rho));
		mfacc = CUMacc + (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho - c1o9*(drho / rho));

		//5.
		mfbcc = CUMbcc + ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac)) / rho;
		mfcbc = CUMcbc + ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc)) / rho;
		mfccb = CUMccb + ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab)) / rho;

		//6.

		mfccc = CUMccc - ((-c4o1 *  mfbbb * mfbbb
			- (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
			- c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
			- c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
			+ (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
				+ c2o1 * (mfcaa * mfaca * mfaac)
				+ c16o1 *  mfbba * mfbab * mfabb) / (rho * rho)
			- c1o3 * (mfacc + mfcac + mfcca) / rho
			- c1o9 * (mfcaa + mfaca + mfaac) / rho
			+ (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
				+ (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3
			+ c1o27*((drho * drho - drho) / (rho*rho)));
		//+ c1o27*(one -three/rho +two/(rho*rho)));
		////////////////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////////
		//the force be with you
		mfbaa = -mfbaa;
		mfaba = -mfaba;
		mfaab = -mfaab;
		////////////////////////////////////////////////////////////////////////////////////


		////////////////////////////////////////////////////////////////////////////////////
		//back
		////////////////////////////////////////////////////////////////////////////////////
		//mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
		////////////////////////////////////////////////////////////////////////////////////
		// Z - Dir
		m0 = mfaac * c1o2 + mfaab * (vvz - c1o2) + (mfaaa + c1o1* oMdrho) * (vz2 - vvz) * c1o2;
		m1 = -mfaac - c2o1* mfaab *  vvz + mfaaa                * (c1o1 - vz2) - c1o1* oMdrho * vz2;
		m2 = mfaac * c1o2 + mfaab * (vvz + c1o2) + (mfaaa + c1o1* oMdrho) * (vz2 + vvz) * c1o2;
		mfaaa = m0;
		mfaab = m1;
		mfaac = m2;
		////////////////////////////////////////////////////////////////////////////////////
		m0 = mfabc * c1o2 + mfabb * (vvz - c1o2) + mfaba * (vz2 - vvz) * c1o2;
		m1 = -mfabc - c2o1* mfabb *  vvz + mfaba * (c1o1 - vz2);
		m2 = mfabc * c1o2 + mfabb * (vvz + c1o2) + mfaba * (vz2 + vvz) * c1o2;
		mfaba = m0;
		mfabb = m1;
		mfabc = m2;
		////////////////////////////////////////////////////////////////////////////////////
		m0 = mfacc * c1o2 + mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (vz2 - vvz) * c1o2;
		m1 = -mfacc - c2o1* mfacb *  vvz + mfaca                  * (c1o1 - vz2) - c1o3 * oMdrho * vz2;
		m2 = mfacc * c1o2 + mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (vz2 + vvz) * c1o2;
		mfaca = m0;
		mfacb = m1;
		mfacc = m2;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		m0 = mfbac * c1o2 + mfbab * (vvz - c1o2) + mfbaa * (vz2 - vvz) * c1o2;
		m1 = -mfbac - c2o1* mfbab *  vvz + mfbaa * (c1o1 - vz2);
		m2 = mfbac * c1o2 + mfbab * (vvz + c1o2) + mfbaa * (vz2 + vvz) * c1o2;
		mfbaa = m0;
		mfbab = m1;
		mfbac = m2;
		/////////b//////////////////////////////////////////////////////////////////////////
		m0 = mfbbc * c1o2 + mfbbb * (vvz - c1o2) + mfbba * (vz2 - vvz) * c1o2;
		m1 = -mfbbc - c2o1* mfbbb *  vvz + mfbba * (c1o1 - vz2);
		m2 = mfbbc * c1o2 + mfbbb * (vvz + c1o2) + mfbba * (vz2 + vvz) * c1o2;
		mfbba = m0;
		mfbbb = m1;
		mfbbc = m2;
		/////////b//////////////////////////////////////////////////////////////////////////
		m0 = mfbcc * c1o2 + mfbcb * (vvz - c1o2) + mfbca * (vz2 - vvz) * c1o2;
		m1 = -mfbcc - c2o1* mfbcb *  vvz + mfbca * (c1o1 - vz2);
		m2 = mfbcc * c1o2 + mfbcb * (vvz + c1o2) + mfbca * (vz2 + vvz) * c1o2;
		mfbca = m0;
		mfbcb = m1;
		mfbcc = m2;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		m0 = mfcac * c1o2 + mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (vz2 - vvz) * c1o2;
		m1 = -mfcac - c2o1* mfcab *  vvz + mfcaa                  * (c1o1 - vz2) - c1o3 * oMdrho * vz2;
		m2 = mfcac * c1o2 + mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (vz2 + vvz) * c1o2;
		mfcaa = m0;
		mfcab = m1;
		mfcac = m2;
		/////////c//////////////////////////////////////////////////////////////////////////
		m0 = mfcbc * c1o2 + mfcbb * (vvz - c1o2) + mfcba * (vz2 - vvz) * c1o2;
		m1 = -mfcbc - c2o1* mfcbb *  vvz + mfcba * (c1o1 - vz2);
		m2 = mfcbc * c1o2 + mfcbb * (vvz + c1o2) + mfcba * (vz2 + vvz) * c1o2;
		mfcba = m0;
		mfcbb = m1;
		mfcbc = m2;
		/////////c//////////////////////////////////////////////////////////////////////////
		m0 = mfccc * c1o2 + mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (vz2 - vvz) * c1o2;
		m1 = -mfccc - c2o1* mfccb *  vvz + mfcca                  * (c1o1 - vz2) - c1o9 * oMdrho * vz2;
		m2 = mfccc * c1o2 + mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (vz2 + vvz) * c1o2;
		mfcca = m0;
		mfccb = m1;
		mfccc = m2;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		//mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
		////////////////////////////////////////////////////////////////////////////////////
		// Y - Dir
		m0 = mfaca * c1o2 + mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (vy2 - vvy) * c1o2;
		m1 = -mfaca - c2o1* mfaba *  vvy + mfaaa                  * (c1o1 - vy2) - c1o6 * oMdrho * vy2;
		m2 = mfaca * c1o2 + mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (vy2 + vvy) * c1o2;
		mfaaa = m0;
		mfaba = m1;
		mfaca = m2;
		////////////////////////////////////////////////////////////////////////////////////
		m0 = mfacb * c1o2 + mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (vy2 - vvy) * c1o2;
		m1 = -mfacb - c2o1* mfabb *  vvy + mfaab                  * (c1o1 - vy2) - c2o3 * oMdrho * vy2;
		m2 = mfacb * c1o2 + mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (vy2 + vvy) * c1o2;
		mfaab = m0;
		mfabb = m1;
		mfacb = m2;
		////////////////////////////////////////////////////////////////////////////////////
		m0 = mfacc * c1o2 + mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (vy2 - vvy) * c1o2;
		m1 = -mfacc - c2o1* mfabc *  vvy + mfaac                  * (c1o1 - vy2) - c1o6 * oMdrho * vy2;
		m2 = mfacc * c1o2 + mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (vy2 + vvy) * c1o2;
		mfaac = m0;
		mfabc = m1;
		mfacc = m2;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		m0 = mfbca * c1o2 + mfbba * (vvy - c1o2) + mfbaa * (vy2 - vvy) * c1o2;
		m1 = -mfbca - c2o1* mfbba *  vvy + mfbaa * (c1o1 - vy2);
		m2 = mfbca * c1o2 + mfbba * (vvy + c1o2) + mfbaa * (vy2 + vvy) * c1o2;
		mfbaa = m0;
		mfbba = m1;
		mfbca = m2;
		/////////b//////////////////////////////////////////////////////////////////////////
		m0 = mfbcb * c1o2 + mfbbb * (vvy - c1o2) + mfbab * (vy2 - vvy) * c1o2;
		m1 = -mfbcb - c2o1* mfbbb *  vvy + mfbab * (c1o1 - vy2);
		m2 = mfbcb * c1o2 + mfbbb * (vvy + c1o2) + mfbab * (vy2 + vvy) * c1o2;
		mfbab = m0;
		mfbbb = m1;
		mfbcb = m2;
		/////////b//////////////////////////////////////////////////////////////////////////
		m0 = mfbcc * c1o2 + mfbbc * (vvy - c1o2) + mfbac * (vy2 - vvy) * c1o2;
		m1 = -mfbcc - c2o1* mfbbc *  vvy + mfbac * (c1o1 - vy2);
		m2 = mfbcc * c1o2 + mfbbc * (vvy + c1o2) + mfbac * (vy2 + vvy) * c1o2;
		mfbac = m0;
		mfbbc = m1;
		mfbcc = m2;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		m0 = mfcca * c1o2 + mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (vy2 - vvy) * c1o2;
		m1 = -mfcca - c2o1* mfcba *  vvy + mfcaa                   * (c1o1 - vy2) - c1o18 * oMdrho * vy2;
		m2 = mfcca * c1o2 + mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (vy2 + vvy) * c1o2;
		mfcaa = m0;
		mfcba = m1;
		mfcca = m2;
		/////////c//////////////////////////////////////////////////////////////////////////
		m0 = mfccb * c1o2 + mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (vy2 - vvy) * c1o2;
		m1 = -mfccb - c2o1* mfcbb *  vvy + mfcab                  * (c1o1 - vy2) - c2o9 * oMdrho * vy2;
		m2 = mfccb * c1o2 + mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (vy2 + vvy) * c1o2;
		mfcab = m0;
		mfcbb = m1;
		mfccb = m2;
		/////////c//////////////////////////////////////////////////////////////////////////
		m0 = mfccc * c1o2 + mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (vy2 - vvy) * c1o2;
		m1 = -mfccc - c2o1* mfcbc *  vvy + mfcac                   * (c1o1 - vy2) - c1o18 * oMdrho * vy2;
		m2 = mfccc * c1o2 + mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (vy2 + vvy) * c1o2;
		mfcac = m0;
		mfcbc = m1;
		mfccc = m2;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		//mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
		////////////////////////////////////////////////////////////////////////////////////
		// X - Dir
		m0 = mfcaa * c1o2 + mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
		m1 = -mfcaa - c2o1* mfbaa *  vvx + mfaaa                   * (c1o1 - vx2) - c1o36 * oMdrho * vx2;
		m2 = mfcaa * c1o2 + mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
		mfaaa = m0;
		mfbaa = m1;
		mfcaa = m2;
		////////////////////////////////////////////////////////////////////////////////////
		m0 = mfcba * c1o2 + mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
		m1 = -mfcba - c2o1* mfbba *  vvx + mfaba                  * (c1o1 - vx2) - c1o9 * oMdrho * vx2;
		m2 = mfcba * c1o2 + mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
		mfaba = m0;
		mfbba = m1;
		mfcba = m2;
		////////////////////////////////////////////////////////////////////////////////////
		m0 = mfcca * c1o2 + mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
		m1 = -mfcca - c2o1* mfbca *  vvx + mfaca                   * (c1o1 - vx2) - c1o36 * oMdrho * vx2;
		m2 = mfcca * c1o2 + mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
		mfaca = m0;
		mfbca = m1;
		mfcca = m2;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		m0 = mfcab * c1o2 + mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
		m1 = -mfcab - c2o1* mfbab *  vvx + mfaab                  * (c1o1 - vx2) - c1o9 * oMdrho * vx2;
		m2 = mfcab * c1o2 + mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
		mfaab = m0;
		mfbab = m1;
		mfcab = m2;
		///////////b////////////////////////////////////////////////////////////////////////
		m0 = mfcbb * c1o2 + mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (vx2 - vvx) * c1o2;
		m1 = -mfcbb - c2o1* mfbbb *  vvx + mfabb                  * (c1o1 - vx2) - c4o9 * oMdrho * vx2;
		m2 = mfcbb * c1o2 + mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (vx2 + vvx) * c1o2;
		mfabb = m0;
		mfbbb = m1;
		mfcbb = m2;
		///////////b////////////////////////////////////////////////////////////////////////
		m0 = mfccb * c1o2 + mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
		m1 = -mfccb - c2o1* mfbcb *  vvx + mfacb                  * (c1o1 - vx2) - c1o9 * oMdrho * vx2;
		m2 = mfccb * c1o2 + mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
		mfacb = m0;
		mfbcb = m1;
		mfccb = m2;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		m0 = mfcac * c1o2 + mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
		m1 = -mfcac - c2o1* mfbac *  vvx + mfaac                   * (c1o1 - vx2) - c1o36 * oMdrho * vx2;
		m2 = mfcac * c1o2 + mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
		mfaac = m0;
		mfbac = m1;
		mfcac = m2;
		///////////c////////////////////////////////////////////////////////////////////////
		m0 = mfcbc * c1o2 + mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
		m1 = -mfcbc - c2o1* mfbbc *  vvx + mfabc                  * (c1o1 - vx2) - c1o9 * oMdrho * vx2;
		m2 = mfcbc * c1o2 + mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
		mfabc = m0;
		mfbbc = m1;
		mfcbc = m2;
		///////////c////////////////////////////////////////////////////////////////////////
		m0 = mfccc * c1o2 + mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
		m1 = -mfccc - c2o1* mfbcc *  vvx + mfacc                   * (c1o1 - vx2) - c1o36 * oMdrho * vx2;
		m2 = mfccc * c1o2 + mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
		mfacc = m0;
		mfbcc = m1;
		mfccc = m2;
		////////////////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////////
		(D.f[DIR_P00])[k] = mfabb;
		(D.f[DIR_M00])[kw] = mfcbb;
		(D.f[DIR_0P0])[k] = mfbab;
		(D.f[DIR_0M0])[ks] = mfbcb;
		(D.f[DIR_00P])[k] = mfbba;
		(D.f[DIR_00M])[kb] = mfbbc;
		(D.f[DIR_PP0])[k] = mfaab;
		(D.f[DIR_MM0])[ksw] = mfccb;
		(D.f[DIR_PM0])[ks] = mfacb;
		(D.f[DIR_MP0])[kw] = mfcab;
		(D.f[DIR_P0P])[k] = mfaba;
		(D.f[DIR_M0M])[kbw] = mfcbc;
		(D.f[DIR_P0M])[kb] = mfabc;
		(D.f[DIR_M0P])[kw] = mfcba;
		(D.f[DIR_0PP])[k] = mfbaa;
		(D.f[DIR_0MM])[kbs] = mfbcc;
		(D.f[DIR_0PM])[kb] = mfbac;
		(D.f[DIR_0MP])[ks] = mfbca;
		(D.f[DIR_000])[k] = mfbbb;
		(D.f[DIR_PPP])[k] = mfaaa;
		(D.f[DIR_PMP])[ks] = mfaca;
		(D.f[DIR_PPM])[kb] = mfaac;
		(D.f[DIR_PMM])[kbs] = mfacc;
		(D.f[DIR_MPP])[kw] = mfcaa;
		(D.f[DIR_MMP])[ksw] = mfcca;
		(D.f[DIR_MPM])[kbw] = mfcac;
		(D.f[DIR_MMM])[kbsw] = mfccc;
		////////////////////////////////////////////////////////////////////////////////////
	}
}