#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include "Core/RealConstants.h"
#include "math.h"

extern "C" __global__ void LB_Kernel_WaleCumulantK17Comp(
	real omega_in,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int* neighborWSB,
	real* veloX,
	real* veloY,
	real* veloZ,
	real* DDStart,
	real* turbulentViscosity,
	int size_Mat,
	int level,
	unsigned int timestep, 
	real* forces,
	real* quadricLimiters,
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

	if(k<size_Mat)
	{
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = bcMatD[k];

		if( BC == GEO_FLUID /*(BC != GEO_SOLID) && (BC != GEO_VOID)*/ )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
			{
				D.f[dirE   ] = &DDStart[dirE   *size_Mat];
				D.f[dirW   ] = &DDStart[dirW   *size_Mat];
				D.f[dirN   ] = &DDStart[dirN   *size_Mat];
				D.f[dirS   ] = &DDStart[dirS   *size_Mat];
				D.f[dirT   ] = &DDStart[dirT   *size_Mat];
				D.f[dirB   ] = &DDStart[dirB   *size_Mat];
				D.f[dirNE  ] = &DDStart[dirNE  *size_Mat];
				D.f[dirSW  ] = &DDStart[dirSW  *size_Mat];
				D.f[dirSE  ] = &DDStart[dirSE  *size_Mat];
				D.f[dirNW  ] = &DDStart[dirNW  *size_Mat];
				D.f[dirTE  ] = &DDStart[dirTE  *size_Mat];
				D.f[dirBW  ] = &DDStart[dirBW  *size_Mat];
				D.f[dirBE  ] = &DDStart[dirBE  *size_Mat];
				D.f[dirTW  ] = &DDStart[dirTW  *size_Mat];
				D.f[dirTN  ] = &DDStart[dirTN  *size_Mat];
				D.f[dirBS  ] = &DDStart[dirBS  *size_Mat];
				D.f[dirBN  ] = &DDStart[dirBN  *size_Mat];
				D.f[dirTS  ] = &DDStart[dirTS  *size_Mat];
				D.f[dirZERO] = &DDStart[dirZERO*size_Mat];
				D.f[dirTNE ] = &DDStart[dirTNE *size_Mat];
				D.f[dirTSW ] = &DDStart[dirTSW *size_Mat];
				D.f[dirTSE ] = &DDStart[dirTSE *size_Mat];
				D.f[dirTNW ] = &DDStart[dirTNW *size_Mat];
				D.f[dirBNE ] = &DDStart[dirBNE *size_Mat];
				D.f[dirBSW ] = &DDStart[dirBSW *size_Mat];
				D.f[dirBSE ] = &DDStart[dirBSE *size_Mat];
				D.f[dirBNW ] = &DDStart[dirBNW *size_Mat];
			}
			else
			{
				D.f[dirW   ] = &DDStart[dirE   *size_Mat];
				D.f[dirE   ] = &DDStart[dirW   *size_Mat];
				D.f[dirS   ] = &DDStart[dirN   *size_Mat];
				D.f[dirN   ] = &DDStart[dirS   *size_Mat];
				D.f[dirB   ] = &DDStart[dirT   *size_Mat];
				D.f[dirT   ] = &DDStart[dirB   *size_Mat];
				D.f[dirSW  ] = &DDStart[dirNE  *size_Mat];
				D.f[dirNE  ] = &DDStart[dirSW  *size_Mat];
				D.f[dirNW  ] = &DDStart[dirSE  *size_Mat];
				D.f[dirSE  ] = &DDStart[dirNW  *size_Mat];
				D.f[dirBW  ] = &DDStart[dirTE  *size_Mat];
				D.f[dirTE  ] = &DDStart[dirBW  *size_Mat];
				D.f[dirTW  ] = &DDStart[dirBE  *size_Mat];
				D.f[dirBE  ] = &DDStart[dirTW  *size_Mat];
				D.f[dirBS  ] = &DDStart[dirTN  *size_Mat];
				D.f[dirTN  ] = &DDStart[dirBS  *size_Mat];
				D.f[dirTS  ] = &DDStart[dirBN  *size_Mat];
				D.f[dirBN  ] = &DDStart[dirTS  *size_Mat];
				D.f[dirZERO] = &DDStart[dirZERO*size_Mat];
				D.f[dirBSW ] = &DDStart[dirTNE *size_Mat];
				D.f[dirBNE ] = &DDStart[dirTSW *size_Mat];
				D.f[dirBNW ] = &DDStart[dirTSE *size_Mat];
				D.f[dirBSE ] = &DDStart[dirTNW *size_Mat];
				D.f[dirTSW ] = &DDStart[dirBNE *size_Mat];
				D.f[dirTNE ] = &DDStart[dirBSW *size_Mat];
				D.f[dirTNW ] = &DDStart[dirBSE *size_Mat];
				D.f[dirTSE ] = &DDStart[dirBNW *size_Mat];
			}

			////////////////////////////////////////////////////////////////////////////////
			//index
			//unsigned int kzero= k;
			//unsigned int ke   = k;
			unsigned int kw   = neighborX[k];
			//unsigned int kn   = k;
			unsigned int ks   = neighborY[k];
			//unsigned int kt   = k;
			unsigned int kb   = neighborZ[k];
			unsigned int ksw  = neighborY[kw];
			//unsigned int kne  = k;
			//unsigned int kse  = ks;
			//unsigned int knw  = kw;
			unsigned int kbw  = neighborZ[kw];
			//unsigned int kte  = k;
			//unsigned int kbe  = kb;
			//unsigned int ktw  = kw;
			unsigned int kbs  = neighborZ[ks];
			//unsigned int ktn  = k;
			//unsigned int kbn  = kb;
			//unsigned int kts  = ks;
			//unsigned int ktse = ks;
			//unsigned int kbnw = kbw;
			//unsigned int ktnw = kw;
			//unsigned int kbse = kbs;
			//unsigned int ktsw = ksw;
			//unsigned int kbne = kb;
			//unsigned int ktne = k;
			unsigned int kbsw = neighborZ[ksw];

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real mfcbb = (D.f[dirE   ])[k  ];
			real mfabb = (D.f[dirW   ])[kw ];
			real mfbcb = (D.f[dirN   ])[k  ];
			real mfbab = (D.f[dirS   ])[ks ];
			real mfbbc = (D.f[dirT   ])[k  ];
			real mfbba = (D.f[dirB   ])[kb ];
			real mfccb = (D.f[dirNE  ])[k  ];
			real mfaab = (D.f[dirSW  ])[ksw];
			real mfcab = (D.f[dirSE  ])[ks ];
			real mfacb = (D.f[dirNW  ])[kw ];
			real mfcbc = (D.f[dirTE  ])[k  ];
			real mfaba = (D.f[dirBW  ])[kbw];
			real mfcba = (D.f[dirBE  ])[kb ];
			real mfabc = (D.f[dirTW  ])[kw ];
			real mfbcc = (D.f[dirTN  ])[k  ];
			real mfbaa = (D.f[dirBS  ])[kbs];
			real mfbca = (D.f[dirBN  ])[kb ];
			real mfbac = (D.f[dirTS  ])[ks ];
			real mfbbb = (D.f[dirZERO])[k  ];
			real mfccc = (D.f[dirTNE ])[k  ];
			real mfaac = (D.f[dirTSW ])[ksw];
			real mfcac = (D.f[dirTSE ])[ks ];
			real mfacc = (D.f[dirTNW ])[kw ];
			real mfcca = (D.f[dirBNE ])[kb ];
			real mfaaa = (D.f[dirBSW])[kbsw];
			real mfcaa = (D.f[dirBSE ])[kbs];
			real mfaca = (D.f[dirBNW ])[kbw];
			////////////////////////////////////////////////////////////////////////////////////
			real drho = ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
							(((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
							((mfabb+mfcbb) + (mfbab+mfbcb) + (mfbba+mfbbc))) + mfbbb;

			real rho = c1o1+drho;
			////////////////////////////////////////////////////////////////////////////////////
			//slow
			real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb)) / rho;
			real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab)) / rho;
			real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				               (mfbbc-mfbba)) / rho;
			////////////////////////////////////////////////////////////////////////////////////
			real nuTurb = c0o1;
			{
				/////////////      Wale Model     ///////////////
				//neighbor index
				unsigned int kPx = neighborX[k];
				unsigned int kPy = neighborY[k];
				unsigned int kPz = neighborZ[k];
				unsigned int kMxyz = neighborWSB[k];
				unsigned int kMx = neighborZ[neighborY[kMxyz]];
				unsigned int kMy = neighborZ[neighborX[kMxyz]];
				unsigned int kMz = neighborY[neighborX[kMxyz]];
				//getVeloX//
				real veloXNeighborPx = veloX[kPx];
				real veloXNeighborMx = veloX[kMx];
				real veloXNeighborPy = veloX[kPy];
				real veloXNeighborMy = veloX[kMy];
				real veloXNeighborPz = veloX[kPz];
				real veloXNeighborMz = veloX[kMz];
				//getVeloY//
				real veloYNeighborPx = veloY[kPx];
				real veloYNeighborMx = veloY[kMx];
				real veloYNeighborPy = veloY[kPy];
				real veloYNeighborMy = veloY[kMy];
				real veloYNeighborPz = veloY[kPz];
				real veloYNeighborMz = veloY[kMz];
				//getVeloZ//
				real veloZNeighborPx = veloZ[kPx];
				real veloZNeighborMx = veloZ[kMx];
				real veloZNeighborPy = veloZ[kPy];
				real veloZNeighborMy = veloZ[kMy];
				real veloZNeighborPz = veloZ[kPz];
				real veloZNeighborMz = veloZ[kMz];
				//getVeloLocal//
				real veloLocalX = veloX[k];
				real veloLocalY = veloY[k];
				real veloLocalZ = veloZ[k];

				//////////////////////////////////////////////////////////////////////////////
				real dxvx = c0o1;
				real dyvx = c0o1;
				real dzvx = c0o1;
				real dxvy = c0o1;
				real dyvy = c0o1;
				real dzvy = c0o1;
				real dxvz = c0o1;
				real dyvz = c0o1;
				real dzvz = c0o1;
				real SumSd = c0o1;
				real SumS = c0o1;
				//////////////////////////////////////////////////////////////////////////////
				//no central differences at the boundary nodes
				//Dx
				if (bcMatD[kPx] != GEO_FLUID)
				{
					dxvx = (veloLocalX - veloXNeighborMx);
					dxvy = (veloLocalY - veloYNeighborMx);
					dxvz = (veloLocalZ - veloZNeighborMx);
				}
				else if (bcMatD[kMx] != GEO_FLUID)
				{
					dxvx = (veloXNeighborPx - veloLocalX);
					dxvy = (veloYNeighborPx - veloLocalY);
					dxvz = (veloZNeighborPx - veloLocalZ);
				}
				else
				{
					dxvx = (veloXNeighborPx - veloXNeighborMx) / c2o1;
					dxvy = (veloYNeighborPx - veloYNeighborMx) / c2o1;
					dxvz = (veloZNeighborPx - veloZNeighborMx) / c2o1;
				}
				//Dy
				if (bcMatD[kPy] != GEO_FLUID)
				{
					dyvx = (veloLocalX - veloXNeighborMy);
					dyvy = (veloLocalY - veloYNeighborMy);
					dyvz = (veloLocalZ - veloZNeighborMy);
				}
				else if (bcMatD[kMy] != GEO_FLUID)
				{
					dyvx = (veloXNeighborPy - veloLocalX);
					dyvy = (veloYNeighborPy - veloLocalY);
					dyvz = (veloZNeighborPy - veloLocalZ);
				}
				else
				{
					dyvx = (veloXNeighborPy - veloXNeighborMy) / c2o1;
					dyvy = (veloYNeighborPy - veloYNeighborMy) / c2o1;
					dyvz = (veloZNeighborPy - veloZNeighborMy) / c2o1;
				}
				//Dz
				if (bcMatD[kPz] != GEO_FLUID)
				{
					dzvx = (veloLocalX - veloXNeighborMz);
					dzvy = (veloLocalY - veloYNeighborMz);
					dzvz = (veloLocalZ - veloZNeighborMz);
				}
				else if (bcMatD[kMz] != GEO_FLUID)
				{
					dzvx = (veloXNeighborPz - veloLocalX);
					dzvy = (veloYNeighborPz - veloLocalY);
					dzvz = (veloZNeighborPz - veloLocalZ);
				}
				else
				{
					dzvx = (veloXNeighborPz - veloXNeighborMz) / c2o1;
					dzvy = (veloYNeighborPz - veloYNeighborMz) / c2o1;
					dzvz = (veloZNeighborPz - veloZNeighborMz) / c2o1;
				}
				//////////////////////////////////////////////////////////////////////////////
				////partial Div vx in x, y, z//
				//dxvx = (veloXNeighborPx - veloXNeighborMx) / two; //deltaX * two??
				//dyvx = (veloXNeighborPy - veloXNeighborMy) / two; //deltaX * two??
				//dzvx = (veloXNeighborPz - veloXNeighborMz) / two; //deltaX * two??
				//												  //partial Div vy in x, y, z//
				//dxvy = (veloYNeighborPx - veloYNeighborMx) / two; //deltaX * two??
				//dyvy = (veloYNeighborPy - veloYNeighborMy) / two; //deltaX * two??
				//dzvy = (veloYNeighborPz - veloYNeighborMz) / two; //deltaX * two??
				//												  //partial Div vz in x, y, z//
				//dxvz = (veloZNeighborPx - veloZNeighborMx) / two; //deltaX * two??
				//dyvz = (veloZNeighborPy - veloZNeighborMy) / two; //deltaX * two??
				//dzvz = (veloZNeighborPz - veloZNeighborMz) / two; //deltaX * two??

				real g11, g12, g13, g21, g22, g23, g31, g32, g33;
				real g11sq, g12sq, g13sq, g21sq, g22sq, g23sq, g31sq, g32sq, g33sq;

				g11 = dxvx;
				g12 = dyvx;
				g13 = dzvx;
				g21 = dxvy;
				g22 = dyvy;
				g23 = dzvy;
				g31 = dxvz;
				g32 = dyvz;
				g33 = dzvz;

				g11sq = g11 * g11 + g12 * g21 + g13 * g31;
				g12sq = g11 * g12 + g12 * g22 + g13 * g32;
				g13sq = g11 * g13 + g12 * g23 + g13 * g33;
				g21sq = g21 * g11 + g22 * g21 + g23 * g31;
				g22sq = g21 * g12 + g22 * g22 + g23 * g32;
				g23sq = g21 * g13 + g22 * g23 + g23 * g33;
				g31sq = g31 * g11 + g32 * g21 + g33 * g31;
				g32sq = g31 * g12 + g32 * g22 + g33 * g32;
				g33sq = g31 * g13 * g32 * g23 + g33 * g33;

				real gkk = g11sq + g22sq + g33sq;

				real Sd11 = c1o2 * (g11sq + g11sq) - c1o3 * gkk;
				real Sd12 = c1o2 * (g12sq + g21sq);
				real Sd13 = c1o2 * (g13sq + g31sq);
				real Sd21 = c1o2 * (g21sq + g12sq); // ==Sd12
				real Sd22 = c1o2 * (g22sq + g22sq) - c1o3 * gkk;
				real Sd23 = c1o2 * (g23sq + g32sq);
				real Sd31 = c1o2 * (g31sq + g13sq); // ==Sd13
				real Sd32 = c1o2 * (g32sq + g23sq); // ==Sd23
				real Sd33 = c1o2 * (g33sq + g33sq) - c1o3 * gkk;

				SumSd = Sd11*Sd11 + Sd12*Sd12 + Sd13*Sd13 + Sd21*Sd21 + Sd22*Sd22 + Sd23*Sd23 + Sd31*Sd31 + Sd32*Sd32 + Sd33*Sd33;



				//real SumSd =
				//	c1o2 * powf(dzvx, four) +
				//	c1o2 * powf(dzvy, four) +
				//	c2o3 * powf(dzvz, four) +
				//	c2o3 * powf(dyvy, four) +
				//	c1o2 * powf(dyvz, four) +
				//	c2o3 * powf(dxvx, four) +
				//	c1o2 * powf(dxvy, four) +
				//	c1o2 * powf(dxvz, four) +
				//	c1o2 * powf(dyvx, four) +
				//	powf(dyvx, two) * powf(dxvy, two) +
				//	powf(dzvx, two) * powf(dxvz, two) +
				//	powf(dzvy, two) * powf(dyvz, two) -
				//	c2o3 * powf(dzvz, two) * powf(dyvy, two) -
				//	c2o3 * powf(dzvz, two) * powf(dxvx, two) -
				//	c2o3 * powf(dyvy, two) * powf(dxvx, two);   //powf
				//real SumSd =
				//	((c1o2 * (((dzvx*dzvx)*(dzvx*dzvx)) + 
				//		((dzvy*dzvy)*(dzvy*dzvy)) + 
				//		((dyvz*dyvz)*(dyvz*dyvz)) + 
				//		((dyvx*dyvx)*(dyvx*dyvx)) + 
				//		((dxvy*dxvy)*(dxvy*dxvy)) + 
				//		((dxvz*dxvz)*(dxvz*dxvz)))) +
				//	(c2o3 * (((dxvx*dxvx)*(dxvx*dxvx)) + 
				//		((dyvy*dyvy)*(dyvy*dyvy)) + 
				//		((dzvz*dzvz)*(dzvz*dzvz)))) +
				//	((dyvx * dyvx) * (dxvy * dxvy)) +
				//	((dzvx * dzvx) * (dxvz * dxvz)) +
				//	((dzvy * dzvy) * (dyvz * dyvz))) -
				//	(c2o3 * ((dzvz * dzvz) * (dyvy * dyvy)) + 
				//		((dzvz * dzvz) * (dxvx * dxvx)) + 
				//		((dyvy * dyvy) * (dxvx * dxvx)));   //explicit
				//real SumSd =
				//	c1o2 * pow(dzvx, four) +
				//	c1o2 * pow(dzvy, four) +
				//	c2o3 * pow(dzvz, four) +
				//	c2o3 * pow(dyvy, four) +
				//	c1o2 * pow(dyvz, four) +
				//	c2o3 * pow(dxvx, four) +
				//	c1o2 * pow(dxvy, four) +
				//	c1o2 * pow(dxvz, four) +
				//	c1o2 * pow(dyvx, four) +
				//	pow(dyvx, two) * pow(dxvy, two) +
				//	pow(dzvx, two) * pow(dxvz, two) +
				//	pow(dzvy, two) * pow(dyvz, two) -
				//	c2o3 * pow(dzvz, two) * pow(dyvy, two) -
				//	c2o3 * pow(dzvz, two) * pow(dxvx, two) -
				//	c2o3 * pow(dyvy, two) * pow(dxvx, two);    //pow
				//SumS
				//real SumS =
				//	powf(dxvx, two) +
				//	powf(dyvy, two) +
				//	powf(dzvz, two) +
				//	c1o2 * powf(dyvx + dxvy, two) +
				//	c1o2 * powf(dzvx + dxvz, two) +
				//	c1o2 * powf(dyvz + dzvy, two);   //powf
				SumS =
					((dxvx * dxvx) +
					(dyvy * dyvy) +
						(dzvz * dzvz)) +
						(c1o2 * (((dyvx + dxvy) * (dyvx + dxvy)) +
					((dzvx + dxvz) * (dzvx + dxvz)) +
							((dyvz + dzvy) * (dyvz + dzvy))));   //explicit
																 //real SumS =
																 //	pow(dxvx, two) +
																 //	pow(dyvy, two) +
																 //	pow(dzvz, two) +
																 //	c1o2 * pow(dyvx + dxvy, two) +
																 //	c1o2 * pow(dzvx + dxvz, two) +
																 //	c1o2 * pow(dyvz + dzvy, two);   //pow
																 //nu turbulent
				real coefficient = 0.5; //0.325; //
				real delta = coefficient * c1o1;
				/////////////////////////////////
				//real SumSsq = SumS * SumS;
				//real SumSDsq = SumSd * SumSd;
				real SumSsq = SumS;
				real SumSDsq = SumSd;
				//nuTurb = powf(delta, two) * powf(SumSDsq, c3o2) / (powf(SumSsq, c5o2) + powf(SumSDsq, c5o4) + smallSingle); //powf
				//nuTurb = pow(delta, two) * pow(SumSDsq, c3o2) / (pow(SumSsq, c5o2) + pow(SumSDsq, c5o4) + smallSingle);     //pow
				//nuTurb = (delta * delta) * (real)pow((double)SumSDsq, 1.5) / ((real)pow((double)SumSsq, 2.5) + (real)pow((double)SumSDsq, 1.25) + smallSingle); //SMversion//
				nuTurb = (delta * delta) * (real)pow((double)SumSDsq, 0.25) /
					((real)pow(((real)pow((double)SumSsq, 0.5) / ((real)pow((double)SumSDsq, 0.25) + c10eM10)), 5.0) + c1o1); //SMversion2//
																																	  /////////////////////////////////
																																	  //nuTurb = rho * powf(delta, two) * powf(SumSd*SumSd, c3o2) / (powf(SumS*SumS, c5o2) + powf(SumSd*SumSd, c5o4));
																																	  //nuTurb = powf(delta, two) * powf(SumSd*SumSd, c3o2) / (powf(SumS*SumS, c5o2) + powf(SumSd*SumSd, c5o4));
			}
			//Test debug
			//if (nuTurb > c1o100) nuTurb = c1o100;
			//if (nuTurb < -c1o100) nuTurb = -c1o100;
			//nuTurb = c1o100;
			turbulentViscosity[k] = nuTurb;
			////////////////////////////////////////////////////////////////////////////////////
			//the force be with you
			{
				real fx = forces[0] / (pow((double)c2o1, (double)level)); //zero;//0.0032653/(pow(two,level)); //0.000000005;//(two/1600000.0) / 120.0; //
				real fy = forces[1] / (pow((double)c2o1, (double)level)); //zero;
				real fz = forces[2] / (pow((double)c2o1, (double)level)); //zero;
				vvx += fx*c1o2;
				vvy += fy*c1o2;
				vvz += fz*c1o2;
			}
			////////////////////////////////////////////////////////////////////////////////////
			real omega = omega_in;
			if (timestep >= 500000)
			{
				real nuOld = c1o3 * (c1o1 / omega_in - c1o2);
				omega = c1o1 / (c3o1 * (nuOld + nuTurb) + c1o2);
			}
			////////////////////////////////////////////////////////////////////////////////////
			//fast
			real oMdrho = c1o1; // comp special
			real m0, m1, m2;	
			real vx2;
			real vy2;
			real vz2;
			vx2=vvx*vvx;
			vy2=vvy*vvy;
			vz2=vvz*vvz;
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
			m2    = mfaaa	+ mfaac;
			m1    = mfaac	- mfaaa;
			m0    = m2		+ mfaab;
			mfaaa = m0;
			m0   += c1o36 * oMdrho;	
			mfaab = m1 -		m0 * vvz;
			mfaac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m2    = mfaaa	+ mfaca;
			m1    = mfaca	- mfaaa;
			m0    = m2		+ mfaba;
			mfaaa = m0;
			m0   += c1o6 * oMdrho;
			mfaba = m1 -		m0 * vvy;
			mfaca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += c1o1* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			real OxxPyyPzz = c1o1;	//set the bulk viscosity one is high / two is very low and zero is (too) high ... (also called omega 2)

			////////////////////////////////////////////////////////////
			//3.
			//////////////////////////////
			real OxyyPxzz  = c8o1*(-c2o1+omega)*(c1o1+c2o1*omega)/(-c8o1-c14o1*omega+c7o1*omega*omega);//one;
			real OxyyMxzz  = c8o1*(-c2o1+omega)*(-c7o1+c4o1*omega)/(c56o1-c50o1*omega+c9o1*omega*omega);//one;
			real Oxyz      = c24o1*(-c2o1+omega)*(-c2o1-c7o1*omega+c3o1*omega*omega)/(c48o1+c152o1*omega-c130o1*omega*omega+c29o1*omega*omega*omega);//one;
			////////////////////////////////////////////////////////////
			//4.
			//////////////////////////////
			real O4        = c1o1;
			//////////////////////////////
			//real O4        = omega;//TRT
			////////////////////////////////////////////////////////////
			//5.
			//////////////////////////////
			real O5        = c1o1;
			////////////////////////////////////////////////////////////
			//6.
			//////////////////////////////
			real O6        = c1o1;
			////////////////////////////////////////////////////////////


			//central moments to cumulants
			//4.
			real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho;
			real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho;
			real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho;
				  	 		
			real CUMcca = mfcca - (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho  - c1o9*(drho/rho));
			real CUMcac = mfcac - (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho  - c1o9*(drho/rho));
			real CUMacc = mfacc - (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho  - c1o9*(drho/rho));

			//5.
			real CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) ) / rho ;
			real CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) ) / rho ;
			real CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) ) / rho ;
			
			//6.

			real CUMccc = mfccc + ((-c4o1 *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     c2o1 * (mfcaa * mfaca * mfaac)
							+ c16o1 *  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3 * (mfacc + mfcac + mfcca) /rho 
							-    c1o9 * (mfcaa + mfaca + mfaac) /rho 
							+(    c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) 
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3 
							+ c1o27*((drho * drho - drho)/(rho*rho)));

			//2.
			// linear combinations
			real mxxPyyPzz = mfcaa + mfaca + mfaac;
			real mxxMyy    = mfcaa - mfaca;
			real mxxMzz	   = mfcaa - mfaac;
			
			////////////////////////////////////////////////////////////////////////////
            real Dxy =-c3o1*omega*mfbba;
            real Dxz =-c3o1*omega*mfbab;
            real Dyz =-c3o1*omega*mfabb;

			//3.
			// linear combinations

			real mxxyPyzz = mfcba + mfabc;
			real mxxyMyzz = mfcba - mfabc;

			real mxxzPyyz = mfcab + mfacb;
			real mxxzMyyz = mfcab - mfacb;

			real mxyyPxzz = mfbca + mfbac;
			real mxyyMxzz = mfbca - mfbac;

 			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 			//incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)
 			
 				real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
 				real dyuy = dxux + omega * c3o2 * mxxMyy;
 				real dzuz = dxux + omega * c3o2 * mxxMzz;
 
 				//relax
				mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
 				mxxMyy    += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
 				mxxMzz    += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);
 
 			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 			////no correction
 			//mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);//-magicBulk*OxxPyyPzz;
 			//mxxMyy    += -(-omega) * (-mxxMyy);
 			//mxxMzz    += -(-omega) * (-mxxMzz);
 			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			mfabb     += omega * (-mfabb);
			mfbab     += omega * (-mfbab);
			mfbba     += omega * (-mfbba);
			//////////////////////////////////////////////////////////////////////////

			// linear combinations back
			mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaca = c1o3 * (-c2o1*  mxxMyy +      mxxMzz + mxxPyyPzz);
			mfaac = c1o3 * (       mxxMyy - c2o1* mxxMzz + mxxPyyPzz);


			//relax
			//////////////////////////////////////////////////////////////////////////
			//das ist der limiter
 			wadjust    = Oxyz+(c1o1-Oxyz)*abs(mfbbb)/(abs(mfbbb)+qudricLimitD);
 			mfbbb     += wadjust * (-mfbbb);
 			wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimitP);
 			mxxyPyzz  += wadjust * (-mxxyPyzz);
 			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimitM);
 			mxxyMyzz  += wadjust * (-mxxyMyzz);
 			wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimitP);
 			mxxzPyyz  += wadjust * (-mxxzPyyz);
 			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimitM);
 			mxxzMyyz  += wadjust * (-mxxzMyyz);
 			wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimitP);
 			mxyyPxzz  += wadjust * (-mxyyPxzz);
 			wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimitM);
 			mxyyMxzz  += wadjust * (-mxyyMxzz);
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

			// linear combinations back
			mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
			mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
			mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
			mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
			mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
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
			real A = (c4o1 + c2o1*omega - c3o1*omega*omega) / (c2o1 - c7o1*omega + c5o1*omega*omega);
			real B = (c4o1 + c28o1*omega - c14o1*omega*omega) / (c6o1 - c21o1*omega + c15o1*omega*omega);
			//////////////////////////////////////////////////////////////////////////
			//ohne limiter
			//CUMacc += O4 * (-CUMacc); 
			//CUMcac += O4 * (-CUMcac); 
			//CUMcca += O4 * (-CUMcca); 
			//CUMbbc += O4 * (-CUMbbc); 
			//CUMbcb += O4 * (-CUMbcb); 
			//CUMcbb += O4 * (-CUMcbb); 
			CUMacc = -O4*(c1o1 / omega - c1o2) * (dyuy + dzuz) * c2o3 * A + (c1o1 - O4) * (CUMacc);
			CUMcac = -O4*(c1o1 / omega - c1o2) * (dxux + dzuz) * c2o3 * A + (c1o1 - O4) * (CUMcac);
			CUMcca = -O4*(c1o1 / omega - c1o2) * (dyuy + dxux) * c2o3 * A + (c1o1 - O4) * (CUMcca);
			CUMbbc = -O4*(c1o1 / omega - c1o2) * Dxy           * c1o3 * B + (c1o1 - O4) * (CUMbbc);
			CUMbcb = -O4*(c1o1 / omega - c1o2) * Dxz           * c1o3 * B + (c1o1 - O4) * (CUMbcb);
			CUMcbb = -O4*(c1o1 / omega - c1o2) * Dyz           * c1o3 * B + (c1o1 - O4) * (CUMcbb);
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
						   
			mfcca = CUMcca + (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho  - c1o9*(drho/rho));
			mfcac = CUMcac + (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho  - c1o9*(drho/rho));
			mfacc = CUMacc + (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho  - c1o9*(drho/rho));

			//5.
			mfbcc = CUMbcc + ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) ) / rho ;
			mfcbc = CUMcbc + ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) ) / rho ;
			mfccb = CUMccb + ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) ) / rho ;
			
			//6.

			mfccc = CUMccc - ((-c4o1 *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     c2o1 * (mfcaa * mfaca * mfaac)
							+ c16o1 *  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3 * (mfacc + mfcac + mfcca) /rho 
							-    c1o9 * (mfcaa + mfaca + mfaac) /rho 
							+(    c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) 
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3 
							+ c1o27*((drho * drho - drho)/(rho*rho)));
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
			m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + c1o1* oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfaac        - c2o1* mfaab *  vvz         +  mfaaa                * (c1o1- vz2)              - c1o1* oMdrho * vz2; 
			m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + c1o1* oMdrho) * (     vz2 + vvz) * c1o2;
			mfaaa = m0;
			mfaab = m1;
			mfaac = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
			m1 = -mfabc        - c2o1* mfabb *  vvz         + mfaba * (c1o1- vz2); 
			m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
			mfaba = m0;
			mfabb = m1;
			mfabc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfacc        - c2o1* mfacb *  vvz         +  mfaca                  * (c1o1- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfaca = m0;
			mfacb = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
			m1 = -mfbac        - c2o1* mfbab *  vvz         + mfbaa * (c1o1- vz2); 
			m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
			mfbaa = m0;
			mfbab = m1;
			mfbac = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
			m1 = -mfbbc        - c2o1* mfbbb *  vvz         + mfbba * (c1o1- vz2); 
			m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
			mfbba = m0;
			mfbbb = m1;
			mfbbc = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
			m1 = -mfbcc        - c2o1* mfbcb *  vvz         + mfbca * (c1o1- vz2); 
			m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
			mfbca = m0;
			mfbcb = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfcac        - c2o1* mfcab *  vvz         +  mfcaa                  * (c1o1- vz2)              - c1o3 * oMdrho * vz2; 
			m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcaa = m0;
			mfcab = m1;
			mfcac = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
			m1 = -mfcbc        - c2o1* mfcbb *  vvz         + mfcba * (c1o1- vz2); 
			m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
			mfcba = m0;
			mfcbb = m1;
			mfcbc = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2; 
			m1 = -mfccc        - c2o1* mfccb *  vvz         +  mfcca                  * (c1o1- vz2)              - c1o9 * oMdrho * vz2; 
			m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 + vvz) * c1o2;
			mfcca = m0;
			mfccb = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfaca        - c2o1* mfaba *  vvy         +  mfaaa                  * (c1o1- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaaa = m0;
			mfaba = m1;
			mfaca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacb        - c2o1* mfabb *  vvy         +  mfaab                  * (c1o1- vy2)              - c2o3 * oMdrho * vy2; 
			m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaab = m0;
			mfabb = m1;
			mfacb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfacc        - c2o1* mfabc *  vvy         +  mfaac                  * (c1o1- vy2)              - c1o6 * oMdrho * vy2; 
			m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfaac = m0;
			mfabc = m1;
			mfacc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
			m1 = -mfbca        - c2o1* mfbba *  vvy         + mfbaa * (c1o1- vy2); 
			m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
			mfbaa = m0;
			mfbba = m1;
			mfbca = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcb        - c2o1* mfbbb *  vvy         + mfbab * (c1o1- vy2); 
			m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
			mfbab = m0;
			mfbbb = m1;
			mfbcb = m2;
			/////////b//////////////////////////////////////////////////////////////////////////
			m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
			m1 = -mfbcc        - c2o1* mfbbc *  vvy         + mfbac * (c1o1- vy2); 
			m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
			mfbac = m0;
			mfbbc = m1;
			mfbcc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfcca        - c2o1* mfcba *  vvy         +  mfcaa                   * (c1o1- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcaa = m0;
			mfcba = m1;
			mfcca = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccb        - c2o1* mfcbb *  vvy         +  mfcab                  * (c1o1- vy2)              - c2o9 * oMdrho * vy2; 
			m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcab = m0;
			mfcbb = m1;
			mfccb = m2;
			/////////c//////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
			m1 = -mfccc        - c2o1* mfcbc *  vvy         +  mfcac                   * (c1o1- vy2)              - c1o18 * oMdrho * vy2; 
			m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
			mfcac = m0;
			mfcbc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			//mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcaa        - c2o1* mfbaa *  vvx         +  mfaaa                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaaa = m0;
			mfbaa = m1;
			mfcaa = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcba        - c2o1* mfbba *  vvx         +  mfaba                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaba = m0;
			mfbba = m1;
			mfcba = m2;
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcca        - c2o1* mfbca *  vvx         +  mfaca                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaca = m0;
			mfbca = m1;
			mfcca = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcab        - c2o1* mfbab *  vvx         +  mfaab                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaab = m0;
			mfbab = m1;
			mfcab = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbb        - c2o1* mfbbb *  vvx         +  mfabb                  * (c1o1- vx2)              - c4o9 * oMdrho * vx2; 
			m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabb = m0;
			mfbbb = m1;
			mfcbb = m2;
			///////////b////////////////////////////////////////////////////////////////////////
			m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccb        - c2o1* mfbcb *  vvx         +  mfacb                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacb = m0;
			mfbcb = m1;
			mfccb = m2;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcac        - c2o1* mfbac *  vvx         +  mfaac                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfaac = m0;
			mfbac = m1;
			mfcac = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfcbc        - c2o1* mfbbc *  vvx         +  mfabc                  * (c1o1- vx2)              - c1o9 * oMdrho * vx2; 
			m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfabc = m0;
			mfbbc = m1;
			mfcbc = m2;
			///////////c////////////////////////////////////////////////////////////////////////
			m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
			m1 = -mfccc        - c2o1* mfbcc *  vvx         +  mfacc                   * (c1o1- vx2)              - c1o36 * oMdrho * vx2; 
			m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
			mfacc = m0;
			mfbcc = m1;
			mfccc = m2;
			////////////////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////////////////
			(D.f[ dirE   ])[k   ] = mfabb;                                                                    
			(D.f[ dirW   ])[kw  ] = mfcbb;                                                                  
			(D.f[ dirN   ])[k   ] = mfbab;
			(D.f[ dirS   ])[ks  ] = mfbcb;
			(D.f[ dirT   ])[k   ] = mfbba;
			(D.f[ dirB   ])[kb  ] = mfbbc;
			(D.f[ dirNE  ])[k   ] = mfaab;
			(D.f[ dirSW  ])[ksw ] = mfccb;
			(D.f[ dirSE  ])[ks  ] = mfacb;
			(D.f[ dirNW  ])[kw  ] = mfcab;
			(D.f[ dirTE  ])[k   ] = mfaba;
			(D.f[ dirBW  ])[kbw ] = mfcbc;
			(D.f[ dirBE  ])[kb  ] = mfabc;
			(D.f[ dirTW  ])[kw  ] = mfcba;
			(D.f[ dirTN  ])[k   ] = mfbaa;
			(D.f[ dirBS  ])[kbs ] = mfbcc;
			(D.f[ dirBN  ])[kb  ] = mfbac;
			(D.f[ dirTS  ])[ks  ] = mfbca;
			(D.f[ dirZERO])[k   ] = mfbbb;
			(D.f[ dirTNE ])[k   ] = mfaaa;
			(D.f[ dirTSE ])[ks  ] = mfaca;
			(D.f[ dirBNE ])[kb  ] = mfaac;
			(D.f[ dirBSE ])[kbs ] = mfacc;
			(D.f[ dirTNW ])[kw  ] = mfcaa;
			(D.f[ dirTSW ])[ksw ] = mfcca;
			(D.f[ dirBNW ])[kbw ] = mfcac;
			(D.f[ dirBSW ])[kbsw] = mfccc;
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}