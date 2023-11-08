#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
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
				D.g[dP00] = &G6[dP00   *size_Mat];
				D.g[dM00] = &G6[dM00   *size_Mat];
				D.g[d0P0] = &G6[d0P0   *size_Mat];
				D.g[d0M0] = &G6[d0M0   *size_Mat];
				D.g[d00P] = &G6[d00P   *size_Mat];
				D.g[d00M] = &G6[d00M   *size_Mat];
			}
			else
			{
				D.g[dM00] = &G6[dP00   *size_Mat];
				D.g[dP00] = &G6[dM00   *size_Mat];
				D.g[d0M0] = &G6[d0P0   *size_Mat];
				D.g[d0P0] = &G6[d0M0   *size_Mat];
				D.g[d00M] = &G6[d00P   *size_Mat];
				D.g[d00P] = &G6[d00M   *size_Mat];
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

			(D.g[dP00])[ke] = 0.0f;
			(D.g[dM00])[kw] = 0.0f;
			(D.g[d0P0])[kn] = 0.0f;
			(D.g[d0M0])[ks] = 0.0f;
			(D.g[d00P])[kt] = 0.0f;
			(D.g[d00M])[kb] = 0.0f;
		}
	}
}