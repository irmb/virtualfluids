//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ / 
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////
#ifndef _LB_H_
#define _LB_H_

//////////////////////////////////////////////////////////////////////////
#define GEO_FLUID_OLD    1
#define GEO_VELO         2
#define GEO_PRESS        4

//////////////////////////
//porous media
#define GEO_PM_0		 5
#define GEO_PM_1		 6
#define GEO_PM_2		 7
//////////////////////////

#define GEO_SOLID       15
#define GEO_VOID        16

#define GEO_FLUID       19
#define OFFSET_BCsInGeo 20
//////////////////////////////////////////////////////////////////////////

#define LES false // LES Simulation

#define STARTOFFX 16
#define STARTOFFY 16
#define STARTOFFZ 16

#define X1PERIODIC true
#define X2PERIODIC true
#define X3PERIODIC true

#define INTERFACE_E 0
#define INTERFACE_W 1
#define INTERFACE_N 2
#define INTERFACE_S 3
#define INTERFACE_T 4
#define INTERFACE_B 5

//////////////////////////////////////////////////////////////////////////
// precision                             (change between double and float)
//
#include "Core/DataTypes.h"
//////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
//#include <boost/shared_ptr.hpp>
//#define  BSP boost::shared_ptr

// Initial condition
typedef struct InitCond{
   real Re;
   real factorPressBC;
   real Diffusivity, Temp, TempBC;
   real RealX, RealY;
   int numprocs, myid, maxdev;
   unsigned int tend, tout, tStartOut, tCalcMedStart, tCalcMedEnd, tDoCheckPoint, tDoRestart;
   unsigned int PressInID, PressOutID;
   unsigned int PressInZ, PressOutZ;
   std::vector<uint> devices;
   std::vector<int> GridX, GridY, GridZ, DistX, DistY, DistZ;
   std::vector<real> scaleLBMtoSI, translateLBMtoSI;
   std::vector<real> minCoordX, minCoordY, minCoordZ, maxCoordX, maxCoordY, maxCoordZ;
   std::vector<bool> NeedInterface;
   std::string fname, oPath, oPrefix;
   std::string geometryFileC, geometryFileM, geometryFileF;
   std::string kFull, geoFull, geoVec, coordX, coordY, coordZ, neighborX, neighborY, neighborZ, neighborWSB, scaleCFC, scaleCFF, scaleFCC, scaleFCF, scaleOffsetCF, scaleOffsetFC;
   std::string noSlipBcPos, noSlipBcQs, noSlipBcValue;
   std::string slipBcPos, slipBcQs, slipBcValue;
   std::string pressBcPos, pressBcQs, pressBcValue;
   std::string geomBoundaryBcQs,velBcQs;
   std::string geomBoundaryBcValues,velBcValues,pressBcValues,noSlipBcValues;
   std::string propellerCylinder, propellerValues, propellerQs, measurePoints;
   std::string inletBcQs, inletBcValues;
   std::string outletBcQs, outletBcValues;
   std::string topBcQs, topBcValues;
   std::string bottomBcQs, bottomBcValues;
   std::string frontBcQs, frontBcValues;
   std::string backBcQs, backBcValues;
   std::string wallBcQs, wallBcValues;
   std::string periodicBcQs, periodicBcValues;
   std::string numberNodes, LBMvsSI;
   std::string cpTop, cpBottom, cpBottom2;
   std::string concentration, streetVelocity;
   std::string geomNormalX, geomNormalY, geomNormalZ, inflowNormalX, inflowNormalY, inflowNormalZ, outflowNormalX, outflowNormalY, outflowNormalZ;
   unsigned int timeStepForMP;
   real clockCycleForMP;
   real vis, vis_ratio;
   real u0, u0_ratio;
   real delta_rho, delta_press;
   bool  printFiles, readGeo, doRestart, doCheckPoint, isGeo, isProp, isCp, calcMedian, GeometryValues, isConc, is2ndOrderMoments, is3rdOrderMoments, isHighOrderMoments, isWale, isMeasurePoints, isInitNeq;
   bool isGeoNormal, isInflowNormal, isOutflowNormal;
   bool simulatePorousMedia;
   bool streetVelocityFile;
} InitCondition;

//Interface Cells
typedef struct ICellCF{
   uint* ICellCFF;
   uint* ICellCFC;
   uint kCF;
} InterpolationCellCF;

typedef struct ICellFC{
   uint* ICellFCF;
   uint* ICellFCC;
   uint kFC;
} InterpolationCellFC;

//Offset of the interface cells at the wall
typedef struct OffCF{
   real* xOffCF;
   real* yOffCF;
   real* zOffCF;
} OffsetCF;

typedef struct OffFC{
   real* xOffFC;
   real* yOffFC;
   real* zOffFC;
} OffsetFC;

// Distribution functions g 6
typedef struct  Distri6 {
	real* g[6];
} Distributions6;

// Distribution functions f 7
typedef struct  Distri7{
   real* f[7];
} Distributions7;

// Distribution functions f 19
typedef struct  Distri19{
   real* f[19];
} Distributions19;

// Distribution functions f 27
typedef struct  Distri27{
   real* f[27];
   ////////////////////////////////////////////////////////////////////////////
   ////Restart
   //friend class boost::serialization::access;
   //template<class Archive>
   //void serialize(Archive & ar, const unsigned int version)
   //{
	  // ar & f[0];
   //}
   ////////////////////////////////////////////////////////////////////////////
} Distributions27;

//Q for second order BCs
typedef struct QforBC{
	//boost::shared_ptr<int> k;
   int* k;
   int* kN;
   long long* valueQ;
   real* qread;
   real* q27[27];
   real* q19[19];
   int kQ;
   int kArray;
   real *Vx, *Vy, *Vz, *deltaVz, *RhoBC;
}QforBoundaryConditions;

//BCTemp
typedef struct TempforBC{
   int* k;
   real* temp;
   int kTemp;
}TempforBoundaryConditions;

//BCTempVel
typedef struct TempVelforBC{
   int* k;
   real* temp;
   real* tempPulse;
   real* velo;
   int kTemp;
}TempVelforBoundaryConditions;

//BCTempPress
typedef struct TempPressforBC{
   int* k;
   real* temp;
   real* velo;
   int kTemp;
}TempPressforBoundaryConditions;

//measurePoints
typedef struct MeasP{
	std::string name;
	uint k;
	std::vector<real> Vx;
	std::vector<real> Vy;
	std::vector<real> Vz;
	std::vector<real> Rho;
	//real* Vx;
	//real* Vy;
	//real* Vz;
	//real* Rho;
}MeasurePoints;

//Process Neighbors
typedef struct PN27{
	real* f[27];
	uint memsizeFs;
	int* index;
	uint memsizeIndex;
	uint rankNeighbor;
	int numberOfNodes;
	int numberOfFs;
}ProcessNeighbor27;

//path line particles
typedef struct PLP{
	bool *stuck, *hot;
	real *coordXabsolut, *coordYabsolut, *coordZabsolut;
	real *coordXlocal,   *coordYlocal,   *coordZlocal;
	real *veloX,         *veloY,         *veloZ;
	real *randomLocationInit;
	uint *timestep;
	uint *ID;
	uint *cellBaseID;
	uint numberOfParticles, numberOfTimestepsParticles;
	uint memSizeID, memSizeTimestep, memSizerealAll, memSizereal, memSizeBool, memSizeBoolBC;
}PathLineParticles;

//////////////////////////////////////////////////////////////////////////
inline int vectorPosition(int i, int j, int k, int Lx, int Ly )
{
   //return((j+15)*(Lx+2*16)+(i+15));
   return((Lx+2*STARTOFFX)*((Ly+2*STARTOFFY)*(k+STARTOFFZ)+(j+STARTOFFY))+(i+STARTOFFX));
}
//////////////////////////////////////////////////////////////////////////


#endif


