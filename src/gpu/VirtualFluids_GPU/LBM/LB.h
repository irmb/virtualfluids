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


#include "Core/DataTypes.h"

#include <string>
#include <vector>


struct InitCondition
{
   real Re;
   real factorPressBC {1.0};
   real Diffusivity {0.001};
   real Temp {0.0};
   real TempBC {1.0};
   real RealX {1.0};
   real RealY {1.0};
   int numprocs {1};
   int myid {0};
   int maxdev {1};
   uint tDoCheckPoint {0};
   uint tDoRestart {0};
   uint tCalcMedStart {0};
   uint tCalcMedEnd {10};
   uint tend {10};
   uint tout {1};
   uint tStartOut {0};
   uint PressInID {0};
   uint PressOutID {0};
   uint PressInZ {1};
   uint PressOutZ {2};
   std::vector<uint> devices {1, 0}; // one device with ID = 0
   std::vector<int> GridX, GridY, GridZ, DistX, DistY, DistZ;
   std::vector<real> scaleLBMtoSI, translateLBMtoSI;
   std::vector<real> minCoordX, minCoordY, minCoordZ, maxCoordX, maxCoordY, maxCoordZ;
   std::string fname, oPath;
   std::string oPrefix {"MyFile"};
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
   uint timeStepForMP {10};
   real clockCycleForMP {1.0};
   real vis {0.001};
   real vis_ratio {1.0};
   real u0 {0.01};
   real u0_ratio {1.0};
   real delta_rho {0.0};
   real delta_press {1.0};
   bool printFiles {false};
   bool doRestart {false};
   bool doCheckPoint {false};
   bool readGeo {false};
   bool isGeo, isProp, isCp;
   bool GeometryValues {false};
   bool is2ndOrderMoments {false};
   bool is3rdOrderMoments {false};
   bool isHighOrderMoments {false};
   bool calcMedian {false};
   bool isConc {false};
   bool isWale {false};
   bool isTurbulentViscosity {false};
   bool isAMD {false};
   real SGSConstant {0.0};
   bool isMeasurePoints {false};
   bool isInitNeq {false};
   bool isGeoNormal, isInflowNormal, isOutflowNormal;
   bool hasWallModelMonitor {false};
   bool simulatePorousMedia {false};
   bool streetVelocityFile {false};
};

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
} Distributions27, DistributionReferences27;

//Q for second order BCs
//! \struct to manage sub-grid-distances (q) for second order Boundary Conditions (BCs)
typedef struct QforBC{
   int* k;
   int* kN;
   long long* valueQ;
   real* qread;
   real* q27[27];
   real* q19[19];
   int numberOfBCnodes=0;
   int kArray;
   real *Vx,      *Vy,      *Vz;
   real *Vx1,     *Vy1,     *Vz1;
   real *deltaVz, *RhoBC;
   real *normalX, *normalY, *normalZ;
}QforBoundaryConditions;

//BCTemp
typedef struct TempforBC{
   int* k;
   real* temp;
   int kTemp=0;
}TempforBoundaryConditions;

//BCTempVel
typedef struct TempVelforBC{
   int* k;
   real* temp;
   real* tempPulse;
   real* velo;
   int kTemp=0;
}TempVelforBoundaryConditions;

//BCTempPress
typedef struct TempPressforBC{
   int* k;
   real* temp;
   real* velo;
   int kTemp=0;
}TempPressforBoundaryConditions;

// Settings for wall model used in StressBC
typedef struct WMparas{
   real* z0;
   int* samplingOffset;
   bool hasMonitor;
   real* u_star;
   real* Fx;
   real* Fy;
   real* Fz;
}WallModelParameters;

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

typedef struct PN_F3 {
	real* g[6];
	uint memsizeGs;
	int* index;
	uint memsizeIndex;
	uint rankNeighbor;
	int numberOfNodes;
	int numberOfGs;
}ProcessNeighborF3;

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


