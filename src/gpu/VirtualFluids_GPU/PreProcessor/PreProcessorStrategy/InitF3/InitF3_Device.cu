#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;
#include "math.h"


__global__ void LB_Init_F3(unsigned int* neighborX,
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
				D.g[DIR_P00] = &G6[DIR_P00   *size_Mat];
				D.g[DIR_M00] = &G6[DIR_M00   *size_Mat];
				D.g[DIR_0P0] = &G6[DIR_0P0   *size_Mat];
				D.g[DIR_0M0] = &G6[DIR_0M0   *size_Mat];
				D.g[DIR_00P] = &G6[DIR_00P   *size_Mat];
				D.g[DIR_00M] = &G6[DIR_00M   *size_Mat];
			}
			else
			{
				D.g[DIR_M00] = &G6[DIR_P00   *size_Mat];
				D.g[DIR_P00] = &G6[DIR_M00   *size_Mat];
				D.g[DIR_0M0] = &G6[DIR_0P0   *size_Mat];
				D.g[DIR_0P0] = &G6[DIR_0M0   *size_Mat];
				D.g[DIR_00M] = &G6[DIR_00P   *size_Mat];
				D.g[DIR_00P] = &G6[DIR_00M   *size_Mat];
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

			(D.g[DIR_P00])[ke] = 0.0f;
			(D.g[DIR_M00])[kw] = 0.0f;
			(D.g[DIR_0P0])[kn] = 0.0f;
			(D.g[DIR_0M0])[ks] = 0.0f;
			(D.g[DIR_00P])[kt] = 0.0f;
			(D.g[DIR_00M])[kb] = 0.0f;
		}
	}
}