#ifndef _LB_H_
#define _LB_H_

#define GEO_FLUID        1

#define GEO_SOLID       15
#define GEO_VOID        16

#define STARTOFFX 1
#define STARTOFFY 1
#define STARTOFFZ 1

#define X1PERIODIC false
#define X2PERIODIC false
#define X3PERIODIC false

//////////////////////////////////////////////////////////////////////////
// precision    (change between double and float)                       //
#define ISFLOAT                                                         //
//#define ISDOUBLE                                                      //
//////////////////////////////////////////////////////////////////////////

#ifdef ISDOUBLE
typedef double doubflo;
#define DOUBFLO double
#endif

#ifdef ISFLOAT
typedef float doubflo;
#define DOUBFLO float
#endif
//////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>


struct InitCondition
{
   doubflo Re;
   doubflo factorPressBC;
   doubflo RealX, RealY;
   int numprocs, myid, maxdev;
   unsigned int tend,tout;
   std::vector<int> devices;
   std::vector<int> GridX, GridY, GridZ, DistX, DistY, DistZ;
   std::vector<doubflo> minCoordX, minCoordY, minCoordZ, maxCoordX, maxCoordY, maxCoordZ;
   std::string fname;
   std::string kFull, geoFull, geoVec, coordX, coordY, coordZ, neighborX, neighborY, neighborZ, neighborWSB;
   std::string noSlipBcPos, noSlipBcQs, noSlipBcValue;
   std::string pressBcPos, pressBcQs, pressBcValue;
   std::string geomBoundaryBcQs,velBcQs;
   std::string geomBoundaryBcValues,velBcValues,pressBcValues,noSlipBcValues;
   std::string inletBcQs, inletBcValues;
   std::string outletBcQs, outletBcValues;
   std::string topBcQs, topBcValues;
   std::string bottomBcQs, bottomBcValues;
   std::string frontBcQs, frontBcValues;
   std::string backBcQs, backBcValues;
   std::string wallBcQs, wallBcValues;
   std::string periodicBcQs, periodicBcValues;
   std::string numberNodes, LBMvsSI;
   doubflo vis, vis_ratio;
   doubflo u0, u0_ratio;
   doubflo delta_rho, delta_press;
   bool  printFiles, readGeo, isGeo, calcMedian, GeometryValues;
};

struct Distributions27
{
   doubflo* f[27];
};

struct QforBoundaryConditions 
{
   int* k;
   int* kN;
   long long* valueQ;
   doubflo* qread;
   doubflo* q27[27];
   doubflo* q19[19];
   int kQ;
   doubflo *Vx, *Vy, *Vz, *deltaVz, *RhoBC;
};

struct ProcessNeighbor27 
{
	doubflo* f[27];
	unsigned int memsizeFs;
	int* index;
	unsigned int memsizeIndex;
	unsigned int rankNeighbor;
	int numberOfNodes;
	int numberOfFs;
};


#endif
