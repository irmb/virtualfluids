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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_Calculation Calculation
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#ifndef CALCULATION_H
#define CALCULATION_H

//////////////////////////
//porous media
#define GEO_PM_0         5
#define GEO_PM_1         6
#define GEO_PM_2         7
//////////////////////////
#define GEO_SOLID       15
#define GEO_VOID        16
#define GEO_FLUID       19
//////////////////////////

#include <basics/DataTypes.h>
#include <lbm/constants/D3Q27.h>

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

struct Distributions27
{
   real* f[27] = { nullptr };
   constexpr Distributions27() = default;
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
   real* q27[27];
   real* q19[19];
   unsigned int numberOfBCnodes=0;
   int kArray;
   real *Vx,      *Vy,      *Vz;
   real *deltaVz, *RhoBC;
   real *normalX, *normalY, *normalZ;
};

struct QforDirectionalBoundaryCondition
{
   int* k;
   int* kN;
   real* q27[27];
   unsigned int numberOfBCnodes = 0;
   real *RhoBC;
   size_t direction;
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

struct AdvectionDiffusionNoFluxBoundaryConditions
{
   int* BCNodeIndices;
   real* q27[27];
   int numberOfBCnodes=0;
};

struct AdvectionDiffusionDirichletBoundaryConditions
{
   int* BCNodeIndices;
   real* q27[27];
   real* concentration;
   real *vx, *vy, *vz;
   int numberOfBCnodes=0;
};

struct AdvectionDiffusionNeumannBoundaryConditions
{
   int* BCNodeIndices;
   real* gradient;
   real* q27[27];
   real *vx, *vy, *vz;
   int numberOfBCnodes=0;
};

struct AdvectionDiffusionFluxBoundaryConditions
{
   int* BCNodeIndices;
   real* q27[27];
   real* normalX, *normalY, *normalZ;
   real* gradient;
   int numberOfBCnodes=0;
};

// Settings for wall model used in StressBC
struct WallModelParameters
{
   real* roughnessLength;
   real* samplingDistance;
   real* vonKarmanConstant;
   uint* samplingIndices;
   bool hasMonitor;
   real* frictionVelocity;
   real* velocityNodeX;
   real* velocityNodeY;
   real* velocityNodeZ;
   real* velocityExchangeLocationX;
   real* velocityExchangeLocationY;
   real* velocityExchangeLocationZ;
   real* forceX;
   real* forceY;
   real* forceZ;
};


struct MeasurePoints
{
   std::string name;
   uint k;
   std::vector<real> Vx;
   std::vector<real> Vy;
   std::vector<real> Vz;
   std::vector<real> Rho;
};


struct ProcessNeighbor27
{
   real* f[27];
   real* fAD[27];
   size_t memsizeFs;
   uint* index;
   size_t memsizeIndex;
   uint rankNeighbor;
   uint numberOfNodes;
   uint numberOfFs;
   ProcessNeighbor27()= default;
   ProcessNeighbor27(uint numberOfNodes, uint rankNeighbor) : 
      memsizeFs(numberOfNodes*sizeof(real)*27), memsizeIndex(numberOfNodes*sizeof(real)), rankNeighbor(rankNeighbor), numberOfNodes(numberOfNodes), numberOfFs(numberOfNodes*27)
   {}
};


#endif

//! \}
