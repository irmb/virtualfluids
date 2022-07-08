#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;
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
				D.g[E] = &G6[E   *size_Mat];
				D.g[W] = &G6[W   *size_Mat];
				D.g[N] = &G6[N   *size_Mat];
				D.g[S] = &G6[S   *size_Mat];
				D.g[T] = &G6[T   *size_Mat];
				D.g[B] = &G6[B   *size_Mat];
			}
			else
			{
				D.g[W] = &G6[E   *size_Mat];
				D.g[E] = &G6[W   *size_Mat];
				D.g[S] = &G6[N   *size_Mat];
				D.g[N] = &G6[S   *size_Mat];
				D.g[B] = &G6[T   *size_Mat];
				D.g[T] = &G6[B   *size_Mat];
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

			(D.g[E])[ke] = 0.0f;
			(D.g[W])[kw] = 0.0f;
			(D.g[N])[kn] = 0.0f;
			(D.g[S])[ks] = 0.0f;
			(D.g[T])[kt] = 0.0f;
			(D.g[B])[kb] = 0.0f;
		}
	}
}