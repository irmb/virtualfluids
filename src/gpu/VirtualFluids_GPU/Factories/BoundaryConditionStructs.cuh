#ifndef BoundaryConditionStructs_H
#define BoundaryConditionStructs_H
#include "Core/DataTypes.h"

//Q for second order BCs
//! \struct to manage sub-grid-distances (q) for second order Boundary Conditions (BCs)
typedef struct QforBC{
   int* k;
   int* kN;
   long long* valueQ;
   real* qread;
   real* q27[27];
   real* q19[19];
   unsigned int numberOfBCnodes=0;
   int kArray;
   real *Vx,      *Vy,      *Vz;
   real *Vx1,     *Vy1,     *Vz1;
   real *deltaVz, *RhoBC;
   real *normalX, *normalY, *normalZ;
}QforBoundaryConditions;

typedef struct QforPrecursorBC{
   int* k;
   int numberOfBCnodes=0;
   int sizeQ;
   int numberOfPrecursorNodes=0;
   uint nPrecursorReads=0;
   uint nTRead;
   size_t numberOfQuantities;
   real* q27[27];
   uint* planeNeighborNT, *planeNeighborNB, *planeNeighborST, *planeNeighborSB;
   real* weightsNT, *weightsNB, *weightsST,  *weightsSB;
   real* last, *current, *next;
   real velocityX, velocityY, velocityZ;
}QforPrecursorBoundaryConditions;

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

#endif