#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
#include "math.h"


extern "C" __global__ void LB_Init_F3(unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int* geoD,
	real* rho,
	real* ux,
	real* uy,
	real* uz,
	unsigned int size_Mat,
	real* G6,
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

	if (k < size_Mat)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = geoD[k];

		if (BC != GEO_SOLID &&  BC != GEO_VOID)
		{
			Distributions6 D;
			if (EvenOrOdd == true)
			{
				D.g[dirE] = &G6[dirE   *size_Mat];
				D.g[dirW] = &G6[dirW   *size_Mat];
				D.g[dirN] = &G6[dirN   *size_Mat];
				D.g[dirS] = &G6[dirS   *size_Mat];
				D.g[dirT] = &G6[dirT   *size_Mat];
				D.g[dirB] = &G6[dirB   *size_Mat];
			}
			else
			{
				D.g[dirW] = &G6[dirE   *size_Mat];
				D.g[dirE] = &G6[dirW   *size_Mat];
				D.g[dirS] = &G6[dirN   *size_Mat];
				D.g[dirN] = &G6[dirS   *size_Mat];
				D.g[dirB] = &G6[dirT   *size_Mat];
				D.g[dirT] = &G6[dirB   *size_Mat];
			}
			//////////////////////////////////////////////////////////////////////////
			//index
			//////////////////////////////////////////////////////////////////////////
			// unsigned int kzero = k;
			unsigned int ke = k;
			unsigned int kw = neighborX[k];
			unsigned int kn = k;
			unsigned int ks = neighborY[k];
			unsigned int kt = k;
			unsigned int kb = neighborZ[k];
			//////////////////////////////////////////////////////////////////////////

			(D.g[dirE])[ke] = 0.0f;
			(D.g[dirW])[kw] = 0.0f;
			(D.g[dirN])[kn] = 0.0f;
			(D.g[dirS])[ks] = 0.0f;
			(D.g[dirT])[kt] = 0.0f;
			(D.g[dirB])[kb] = 0.0f;
		}
	}
}