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
#define ISFLOAT
//#define ISDOUBLE
///////////////////////////
#ifdef ISDOUBLE
typedef double doubflo;
#define DOUBFLO double
#endif
///////////////////////////
#ifdef ISFLOAT
typedef float doubflo;
#define DOUBFLO float
#endif
//////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
//#include <boost/shared_ptr.hpp>
//#define  BSP boost::shared_ptr

// Initial condition
typedef struct InitCond{
   doubflo Re;
   doubflo factorPressBC;
   doubflo Diffusivity, Temp, TempBC;
   doubflo RealX, RealY;
   int numprocs, myid, maxdev;
   unsigned int tend, tout, tStartOut, tCalcMedStart, tCalcMedEnd, tDoCheckPoint, tDoRestart;
   unsigned int PressInID, PressOutID;
   unsigned int PressInZ, PressOutZ;
   std::vector<int> devices;
   std::vector<int> GridX, GridY, GridZ, DistX, DistY, DistZ;
   std::vector<doubflo> scaleLBMtoSI, translateLBMtoSI;
   std::vector<doubflo> minCoordX, minCoordY, minCoordZ, maxCoordX, maxCoordY, maxCoordZ;
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
   std::string concentration;
   std::string geomNormalX, geomNormalY, geomNormalZ, inflowNormalX, inflowNormalY, inflowNormalZ, outflowNormalX, outflowNormalY, outflowNormalZ;
   unsigned int timeStepForMP;
   doubflo clockCycleForMP;
   doubflo vis, vis_ratio;
   doubflo u0, u0_ratio;
   doubflo delta_rho, delta_press;
   bool  printFiles, readGeo, doRestart, doCheckPoint, isGeo, isProp, isCp, calcMedian, GeometryValues, isConc, is2ndOrderMoments, is3rdOrderMoments, isHighOrderMoments, isWale, isMeasurePoints;
   bool isGeoNormal, isInflowNormal, isOutflowNormal;
} InitCondition;

//Interface Cells
typedef struct ICellCF{
   unsigned int* ICellCFF;
   unsigned int* ICellCFC;
   unsigned int kCF;
} InterpolationCellCF;

typedef struct ICellFC{
   unsigned int* ICellFCF;
   unsigned int* ICellFCC;
   unsigned int kFC;
} InterpolationCellFC;

//Offset of the interface cells at the wall
typedef struct OffCF{
   doubflo* xOffCF;
   doubflo* yOffCF;
   doubflo* zOffCF;
} OffsetCF;

typedef struct OffFC{
   doubflo* xOffFC;
   doubflo* yOffFC;
   doubflo* zOffFC;
} OffsetFC;

// Distribution functions g 6
typedef struct  Distri6 {
	doubflo* g[6];
} Distributions6;

// Distribution functions f 7
typedef struct  Distri7{
   doubflo* f[7];
} Distributions7;

// Distribution functions f 19
typedef struct  Distri19{
   doubflo* f[19];
} Distributions19;

// Distribution functions f 27
typedef struct  Distri27{
   doubflo* f[27];
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
   doubflo* qread;
   doubflo* q27[27];
   doubflo* q19[19];
   int kQ;
   int kArray;
   doubflo *Vx, *Vy, *Vz, *deltaVz, *RhoBC;
}QforBoundaryConditions;

//BCTemp
typedef struct TempforBC{
   int* k;
   doubflo* temp;
   int kTemp;
}TempforBoundaryConditions;

//BCTempVel
typedef struct TempVelforBC{
   int* k;
   doubflo* temp;
   doubflo* tempPulse;
   doubflo* velo;
   int kTemp;
}TempVelforBoundaryConditions;

//BCTempPress
typedef struct TempPressforBC{
   int* k;
   doubflo* temp;
   doubflo* velo;
   int kTemp;
}TempPressforBoundaryConditions;

//measurePoints
typedef struct MeasP{
	std::string name;
	unsigned int k;
	std::vector<doubflo> Vx;
	std::vector<doubflo> Vy;
	std::vector<doubflo> Vz;
	std::vector<doubflo> Rho;
	//doubflo* Vx;
	//doubflo* Vy;
	//doubflo* Vz;
	//doubflo* Rho;
}MeasurePoints;

//Process Neighbors
typedef struct PN27{
	doubflo* f[27];
	unsigned int memsizeFs;
	int* index;
	unsigned int memsizeIndex;
	unsigned int rankNeighbor;
	int numberOfNodes;
	int numberOfFs;
}ProcessNeighbor27;

//path line particles
typedef struct PLP{
	bool    *stuck, *hot;
	doubflo *coordXabsolut, *coordYabsolut, *coordZabsolut;
	doubflo *coordXlocal,   *coordYlocal,   *coordZlocal;
	doubflo *veloX,         *veloY,         *veloZ;
	doubflo *randomLocationInit;
	unsigned int *timestep;
	unsigned int *ID;
	unsigned int *cellBaseID;
	unsigned int numberOfParticles, numberOfTimestepsParticles;
	unsigned int memSizeID, memSizeTimestep, memSizeDoubfloAll, memSizeDoubflo, memSizeBool, memSizeBoolBC;
}PathLineParticles;

//////////////////////////////////////////////////////////////////////////
inline int vectorPosition(int i, int j, int k, int Lx, int Ly )
{
   //return((j+15)*(Lx+2*16)+(i+15));
   return((Lx+2*STARTOFFX)*((Ly+2*STARTOFFY)*(k+STARTOFFZ)+(j+STARTOFFY))+(i+STARTOFFX));
}
//////////////////////////////////////////////////////////////////////////


#endif


