//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file lb.h
//! \ingroup LBM
//! \author Martin Schoenherr
//=======================================================================================#ifndef _LB_H_
#ifndef _LB_H_
#define _LB_H_

//////////////////////////////////////////////////////////////////////////
#define GEO_FLUID_OLD    1
#define GEO_VELO         2
#define GEO_PRESS        4

//////////////////////////
//porous media
#define GEO_PM_0         5
#define GEO_PM_1         6
#define GEO_PM_2         7
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


#include <basics/DataTypes.h>

#include <string>
#include <vector>


//! \brief An enumeration for selecting a template of the collision kernel (CumulantK17)
enum class CollisionTemplate {
   //! - Default: plain collision without additional read/write
   Default,
   //!  - WriteMacroVars: collision \w write out macroscopic variables
   WriteMacroVars,
   //! - ApplyBodyForce: collision \w read and apply body force in the collision kernel
   ApplyBodyForce,
   //! - AllFeatures: collision \w write out macroscopic variables AND read and apply body force
   AllFeatures,
   //! - Border: collision on border nodes
   SubDomainBorder
};
constexpr std::initializer_list<CollisionTemplate> all_CollisionTemplate  = { CollisionTemplate::Default, CollisionTemplate::WriteMacroVars, CollisionTemplate::ApplyBodyForce, CollisionTemplate::AllFeatures, CollisionTemplate::SubDomainBorder};
constexpr std::initializer_list<CollisionTemplate> bulk_CollisionTemplate = { CollisionTemplate::Default, CollisionTemplate::WriteMacroVars, CollisionTemplate::ApplyBodyForce, CollisionTemplate::AllFeatures};

//Interface Cells
// example of old names (pre 2023) ICellCFC: interpolation from Coarse (C) to Fine (F), indices of the Coarse cells (C)
struct ICells
{
   uint* fineCellIndices;
   uint* coarseCellIndices;
   uint numberOfCells;
};

using InterpolationCells = ICells;

//! \brief stores location of neighboring cell (necessary for refinement into the wall)
struct ICellNeigh
{
   real* x;
   real* y;
   real* z;
};

using InterpolationCellNeighbor = ICellNeigh;

// ADD IN FUTURE RELEASE
struct Distributions6 
{
   real* g[6];
};

// ADD IN FUTURE RELEASE
struct  Distributions7
{
   real* f[7];
};

// DEPRECATED
struct  Distributions19
{
   real* f[19];
};

struct  Distributions27
{
   real* f[27];
};

using DistributionReferences27 = Distributions27;

struct SubgridDistances27
{
   real* q[27];
};

//Q for second order BCs
//! \struct to manage sub-grid-distances (q) for second order Boundary Conditions (BCs)
struct QforBoundaryConditions
{
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
};

struct QforPrecursorBoundaryConditions
{
   int* k;
   int numberOfBCnodes=0;
   int sizeQ;
   int numberOfPrecursorNodes=0;
   uint streamIndex=0;
   uint nPrecursorReads=0;
   uint timeStepsBetweenReads;
   size_t numberOfQuantities;
   real* q27[27];
   uint* planeNeighbor0PP, *planeNeighbor0PM, *planeNeighbor0MP, *planeNeighbor0MM;
   real* weights0PP, *weights0PM, *weights0MP,  *weights0MM;
   real* last, *current, *next;
   real velocityX, velocityY, velocityZ;
};

// ADD IN FUTURE RELEASE
struct TempforBoundaryConditions
{
   int* k;
   real* temp;
   int kTemp=0;
};

// ADD IN FUTURE RELEASE
struct TempVelforBoundaryConditions
{
   int* k;
   real* temp;
   real* tempPulse;
   real* velo;
   int kTemp=0;
};

// ADD IN FUTURE RELEASE
struct TempPressforBoundaryConditions
{
   int* k;
   real* temp;
   real* velo;
   int kTemp=0;
};

// Settings for wall model used in StressBC
struct WallModelParameters
{
   real* z0;
   int* samplingOffset;
   bool hasMonitor;
   real* u_star;
   real* Fx;
   real* Fy;
   real* Fz;
};


// ADD IN FUTURE RELEASE
struct MeasurePoints
{
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
};


struct ProcessNeighbor27
{
   real* f[27];
   uint memsizeFs;
   int* index;
   uint memsizeIndex;
   uint rankNeighbor;
   int numberOfNodes;
   int numberOfFs;
};

// ADD IN FUTURE RELEASE
struct ProcessNeighborF3
{
   real* g[6];
   uint memsizeGs;
   int* index;
   uint memsizeIndex;
   uint rankNeighbor;
   int numberOfNodes;
   int numberOfGs;
};

//////////////////////////////////////////////////////////////////////////
// DEPRECATED
inline int vectorPosition(int i, int j, int k, int Lx, int Ly)
{
   return((Lx+2*STARTOFFX)*((Ly+2*STARTOFFY)*(k+STARTOFFZ)+(j+STARTOFFY))+(i+STARTOFFX));
}
//////////////////////////////////////////////////////////////////////////

#endif
