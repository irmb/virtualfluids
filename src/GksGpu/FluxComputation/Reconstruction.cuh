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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
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
//! \file Reconstruction.cuh
//! \ingroup FluxComputation
//! \author Stephan Lenz
//=======================================================================================
#ifndef Reconstruction_CUH
#define Reconstruction_CUH

#include "VirtualFluidsDefinitions.h"

#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"
#include "Parameters/Parameters.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! \brief Reads the cell indices of the cells on positive and negative side of the face
//!
//! \param[in]  faceIndex       index of the current face
//! \param[in]  dataBase        \ref DataBaseStruct that holds the memory pointers
//! \param[out] posCellIndexN   index of the cell on the positive side of the face
//! \param[out] negCellIndexN   index of the cell on the negative side of the face
__host__ __device__ inline void getCellIndicesN ( const uint faceIndex,
                                                  const DataBaseStruct& dataBase,
                                                  uint& posCellIndexN,
                                                  uint& negCellIndexN )
{
    posCellIndexN = dataBase.faceToCell[ POS_CELL( faceIndex, dataBase.numberOfFaces ) ];
    negCellIndexN = dataBase.faceToCell[ NEG_CELL( faceIndex, dataBase.numberOfFaces ) ];
}

//! \brief Reads the cell indices of the tangential cells in x direction on positive and
//! negative side of the face
//!
//! The indices are obtained by pointer chasing. First one goes to the positive/negative 
//! cell. Then from there to the adjacent cells in x direction.
//! 
//! \param[in]  faceIndex       index of the current face
//! \param[in]  dataBase        \ref DataBaseStruct that holds the memory pointers
//! \param[in]  posCellIndexN   index of the cell on the positive side of the face
//! \param[in]  negCellIndexN   index of the cell on the negative side of the face
//! \param[out] posCellIndexTX  two cell indices of the cells that are neighbors in positive x direction of the positive and negative cells
//! \param[out] negCellIndexTX  two cell indices of the cells that are neighbors in negative x direction of the positive and negative cells
__host__ __device__ inline void getCellIndicesTX( const uint faceIndex,
                                                  const DataBaseStruct& dataBase,
                                                  const uint posCellIndexN,
                                                  const uint negCellIndexN,
                                                  uint* posCellIndexTX,
                                                  uint* negCellIndexTX )
{
    posCellIndexTX[0] = dataBase.cellToCell[ CELL_TO_CELL( posCellIndexN, 0, dataBase.numberOfCells ) ];
    posCellIndexTX[1] = dataBase.cellToCell[ CELL_TO_CELL( negCellIndexN, 0, dataBase.numberOfCells ) ];

    negCellIndexTX[0] = dataBase.cellToCell[ CELL_TO_CELL( posCellIndexN, 1, dataBase.numberOfCells ) ];
    negCellIndexTX[1] = dataBase.cellToCell[ CELL_TO_CELL( negCellIndexN, 1, dataBase.numberOfCells ) ];
}

//! \brief Reads the cell indices of the tangential cells in y direction on positive and
//! negative side of the face
//!
//! The indices are obtained by pointer chasing. First one goes to the positive/negative 
//! cell. Then from there to the adjacent cells in y direction.
//! 
//! \param[in]  faceIndex       index of the current face
//! \param[in]  dataBase        \ref DataBaseStruct that holds the memory pointers
//! \param[in]  posCellIndexN   index of the cell on the positive side of the face
//! \param[in]  negCellIndexN   index of the cell on the negative side of the face
//! \param[out] posCellIndexTY  two cell indices of the cells that are neighbors in positive y direction of the positive and negative cells
//! \param[out] negCellIndexTY  two cell indices of the cells that are neighbors in negative y direction of the positive and negative cells
__host__ __device__ inline void getCellIndicesTY( const uint faceIndex,
                                                  const DataBaseStruct& dataBase,
                                                  const uint posCellIndexN,
                                                  const uint negCellIndexN,
                                                  uint* posCellIndexTY,
                                                  uint* negCellIndexTY )
{
    posCellIndexTY[0] = dataBase.cellToCell[ CELL_TO_CELL( posCellIndexN, 2, dataBase.numberOfCells ) ];
    posCellIndexTY[1] = dataBase.cellToCell[ CELL_TO_CELL( negCellIndexN, 2, dataBase.numberOfCells ) ];

    negCellIndexTY[0] = dataBase.cellToCell[ CELL_TO_CELL( posCellIndexN, 3, dataBase.numberOfCells ) ];
    negCellIndexTY[1] = dataBase.cellToCell[ CELL_TO_CELL( negCellIndexN, 3, dataBase.numberOfCells ) ];
}

//! \brief Reads the cell indices of the tangential cells in z direction on positive and
//! negative side of the face
//!
//! The indices are obtained by pointer chasing. First one goes to the positive/negative 
//! cell. Then from there to the adjacent cells in z direction.
//! 
//! \param[in]  faceIndex       index of the current face
//! \param[in]  dataBase        \ref DataBaseStruct that holds the memory pointers
//! \param[in]  posCellIndexN   index of the cell on the positive side of the face
//! \param[in]  negCellIndexN   index of the cell on the negative side of the face
//! \param[out] posCellIndexTZ  two cell indices of the cells that are neighbors in positive z direction of the positive and negative cells
//! \param[out] negCellIndexTZ  two cell indices of the cells that are neighbors in negative z direction of the positive and negative cells
__host__ __device__ inline void getCellIndicesTZ( const uint faceIndex,
                                                  const DataBaseStruct& dataBase,
                                                  const uint posCellIndexN,
                                                  const uint negCellIndexN,
                                                  uint* posCellIndexTZ,
                                                  uint* negCellIndexTZ )
{
    posCellIndexTZ[0] = dataBase.cellToCell[ CELL_TO_CELL( posCellIndexN, 4, dataBase.numberOfCells ) ];
    posCellIndexTZ[1] = dataBase.cellToCell[ CELL_TO_CELL( negCellIndexN, 4, dataBase.numberOfCells ) ];

    negCellIndexTZ[0] = dataBase.cellToCell[ CELL_TO_CELL( posCellIndexN, 5, dataBase.numberOfCells ) ];
    negCellIndexTZ[1] = dataBase.cellToCell[ CELL_TO_CELL( negCellIndexN, 5, dataBase.numberOfCells ) ];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! \brief Computes the flow state on the cell face
//! 
//! Linear interpolation of the conserved variables is used. Due to the uniform grid
//! this is a simple average of the two adjacent cell states.
//!
//! \param[in]  posCons    \ref ConservedVariables object with flow state on positive side of the face
//! \param[in]  negCons    \ref ConservedVariables object with flow state on negative side of the face
//! \param[out] faceCons   \ref ConservedVariables object with flow state on the face
__host__ __device__ inline void computeFaceCons( const ConservedVariables& posCons,
                                                 const ConservedVariables& negCons,
                                                 ConservedVariables& faceCons )
{
    faceCons.rho  = c1o2 * ( negCons.rho  + posCons.rho  );
    faceCons.rhoU = c1o2 * ( negCons.rhoU + posCons.rhoU );
    faceCons.rhoV = c1o2 * ( negCons.rhoV + posCons.rhoV );
    faceCons.rhoW = c1o2 * ( negCons.rhoW + posCons.rhoW );
    faceCons.rhoE = c1o2 * ( negCons.rhoE + posCons.rhoE );
#ifdef USE_PASSIVE_SCALAR
	faceCons.rhoS_1 = c1o2 * ( negCons.rhoS_1 + posCons.rhoS_1 );
	faceCons.rhoS_2 = c1o2 * ( negCons.rhoS_2 + posCons.rhoS_2 );
#endif // USE_PASSIVE_SCALAR
}

//! \brief Computes the normal derivative of the conserved variables
//! on the face by central finite difference
//! 
//! \param[in]  parameters   \ref Parameters struct
//! \param[in]  posCons      \ref ConservedVariables object with flow state on positive side of the face
//! \param[in]  negCons      \ref ConservedVariables object with flow state on negative side of the face
//! \param[in]  facePrim     \ref PrimitiveVariables object with flow state on the face
//! \param[out] gradN        \ref ConservedVariables object with normal derivative of the flow state in terms of conserved variables and divided by density
__host__ __device__ inline void computeGradN( const Parameters& parameters,
                                              const ConservedVariables& posCons,
                                              const ConservedVariables& negCons,
                                              const PrimitiveVariables& facePrim,
                                              ConservedVariables& gradN )
{
    gradN.rho  = ( posCons.rho  - negCons.rho  ) / ( parameters.dx * facePrim.rho );
    gradN.rhoU = ( posCons.rhoU - negCons.rhoU ) / ( parameters.dx * facePrim.rho );
    gradN.rhoV = ( posCons.rhoV - negCons.rhoV ) / ( parameters.dx * facePrim.rho );
    gradN.rhoW = ( posCons.rhoW - negCons.rhoW ) / ( parameters.dx * facePrim.rho );
    gradN.rhoE = ( posCons.rhoE - negCons.rhoE ) / ( parameters.dx * facePrim.rho );
#ifdef USE_PASSIVE_SCALAR
	gradN.rhoS_1 = ( posCons.rhoS_1 - negCons.rhoS_1 ) / ( parameters.dx * facePrim.rho );
	gradN.rhoS_2 = ( posCons.rhoS_2 - negCons.rhoS_2 ) / ( parameters.dx * facePrim.rho );
#endif // USE_PASSIVE_SCALAR
}

//! \brief Computes the tangential derivative of the conserved variables
//! on the face by central finite difference
//! 
//! This function can be used for all three directions, by providing suitable 
//! cell indices.
//! 
//! \param[in]  dataBase          \ref DataBaseStruct that holds the memory pointers
//! \param[in]  parameters        \ref Parameters struct
//! \param[in]  posCellIndexT     two ConservedVariables object with flow state in positive direction
//! \param[in]  negCellIndexT     two ConservedVariables object with flow state in negative direction
//! \param[in]  facePrim          \ref PrimitiveVariables object with flow state on the face
//! \param[out] gradN             \ref ConservedVariables object with derivative of the flow state in terms of conserved variables and divided by density
__host__ __device__ inline void computeGradT( const DataBaseStruct& dataBase,
                                              const Parameters& parameters,
                                              const uint posCellIndexT[2],
                                              const uint negCellIndexT[2],
                                              const PrimitiveVariables& facePrim,
                                              ConservedVariables& gradN )
{
    ConservedVariables cons;

    //////////////////////////////////////////////////////////////////////////
    {
        readCellData(posCellIndexT[0], dataBase, cons);

        gradN.rho  += c1o2 * cons.rho;
        gradN.rhoU += c1o2 * cons.rhoU;
        gradN.rhoV += c1o2 * cons.rhoV;
        gradN.rhoW += c1o2 * cons.rhoW;
        gradN.rhoE += c1o2 * cons.rhoE;
    #ifdef USE_PASSIVE_SCALAR
        gradN.rhoS_1 += c1o2 * cons.rhoS_1;
        gradN.rhoS_2 += c1o2 * cons.rhoS_2;
    #endif // USE_PASSIVE_SCALAR
    }
    {
        readCellData(posCellIndexT[1], dataBase, cons);

        gradN.rho  += c1o2 * cons.rho;
        gradN.rhoU += c1o2 * cons.rhoU;
        gradN.rhoV += c1o2 * cons.rhoV;
        gradN.rhoW += c1o2 * cons.rhoW;
        gradN.rhoE += c1o2 * cons.rhoE;
    #ifdef USE_PASSIVE_SCALAR
        gradN.rhoS_1 += c1o2 * cons.rhoS_1;
        gradN.rhoS_2 += c1o2 * cons.rhoS_2;
    #endif // USE_PASSIVE_SCALAR
    }
    //////////////////////////////////////////////////////////////////////////
    {
        readCellData(negCellIndexT[0], dataBase, cons);

        gradN.rho  -= c1o2 * cons.rho;
        gradN.rhoU -= c1o2 * cons.rhoU;
        gradN.rhoV -= c1o2 * cons.rhoV;
        gradN.rhoW -= c1o2 * cons.rhoW;
        gradN.rhoE -= c1o2 * cons.rhoE;
    #ifdef USE_PASSIVE_SCALAR
        gradN.rhoS_1 -= c1o2 * cons.rhoS_1;
        gradN.rhoS_2 -= c1o2 * cons.rhoS_2;
    #endif // USE_PASSIVE_SCALAR
    }
    {
        readCellData(negCellIndexT[1], dataBase, cons);

        gradN.rho  -= c1o2 * cons.rho;
        gradN.rhoU -= c1o2 * cons.rhoU;
        gradN.rhoV -= c1o2 * cons.rhoV;
        gradN.rhoW -= c1o2 * cons.rhoW;
        gradN.rhoE -= c1o2 * cons.rhoE;
    #ifdef USE_PASSIVE_SCALAR
        gradN.rhoS_1 -= c1o2 * cons.rhoS_1;
        gradN.rhoS_2 -= c1o2 * cons.rhoS_2;
    #endif // USE_PASSIVE_SCALAR
    }
    //////////////////////////////////////////////////////////////////////////
    {
        gradN.rho  /= c2o1 * parameters.dx * facePrim.rho;
        gradN.rhoU /= c2o1 * parameters.dx * facePrim.rho;
        gradN.rhoV /= c2o1 * parameters.dx * facePrim.rho;
        gradN.rhoW /= c2o1 * parameters.dx * facePrim.rho;
        gradN.rhoE /= c2o1 * parameters.dx * facePrim.rho;
    #ifdef USE_PASSIVE_SCALAR
        gradN.rhoS_1 /= c2o1 * parameters.dx * facePrim.rho;
        gradN.rhoS_2 /= c2o1 * parameters.dx * facePrim.rho;
    #endif // USE_PASSIVE_SCALAR
    }
    //////////////////////////////////////////////////////////////////////////
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! \brief Computes the linear reconstruction of the flow field around the face
//!
//! \param[in]  faceIndex        index of the current face
//! \param[in]  dataBase         \ref DataBaseStruct that holds the memory pointers
//! \param[in]  parameters       \ref Parameters struct
//! \param[in]  direction        char with 'x', 'y' or 'z', used for deciding on tangential derivative directions
//! \param[out] gradN            \ref ConservedVariables object with normal derivative of the flow state in terms of conserved variables and divided by density
//! \param[out] gradT1           \ref ConservedVariables object with first  tangential derivative of the flow state in terms of conserved variables and divided by density
//! \param[out] gradT2           \ref ConservedVariables object with second tangential derivative of the flow state in terms of conserved variables and divided by density
//! \param[out] facePrim         \ref PrimitiveVariables object with flow state on the face
//! \param[out] K                deprecated
__host__ __device__ inline void reconstructFiniteDifferences( const uint faceIndex,
                                                              const DataBaseStruct& dataBase,
                                                              const Parameters& parameters,
                                                              const char direction,
                                                              ConservedVariables& gradN,
                                                              ConservedVariables& gradT1,
                                                              ConservedVariables& gradT2,
                                                              PrimitiveVariables& facePrim,
                                                              real& K )
{
    uint posCellIndexN, negCellIndexN;

    getCellIndicesN( faceIndex, dataBase, posCellIndexN, negCellIndexN );
    
    {
        ConservedVariables posCons, negCons, faceCons;

        readCellData(posCellIndexN, dataBase, posCons);
        readCellData(negCellIndexN, dataBase, negCons);
        
        computeFaceCons(posCons, negCons, faceCons);

        facePrim = toPrimitiveVariables( faceCons, K, false );

        computeGradN( parameters, posCons, negCons, facePrim, gradN );
    }

    {
        uint posCellIndexT1[2];
        uint negCellIndexT1[2];
    
        if( direction == 'x' ) getCellIndicesTY(faceIndex, dataBase, posCellIndexN, negCellIndexN, posCellIndexT1, negCellIndexT1);
        if( direction == 'y' ) getCellIndicesTZ(faceIndex, dataBase, posCellIndexN, negCellIndexN, posCellIndexT1, negCellIndexT1);
        if( direction == 'z' ) getCellIndicesTX(faceIndex, dataBase, posCellIndexN, negCellIndexN, posCellIndexT1, negCellIndexT1);

        computeGradT( dataBase, parameters, posCellIndexT1, negCellIndexT1, facePrim, gradT1 );
    }

    {
        uint posCellIndexT2[2];
        uint negCellIndexT2[2];
    
        if( direction == 'x' ) getCellIndicesTZ(faceIndex, dataBase, posCellIndexN, negCellIndexN, posCellIndexT2, negCellIndexT2);
        if( direction == 'y' ) getCellIndicesTX(faceIndex, dataBase, posCellIndexN, negCellIndexN, posCellIndexT2, negCellIndexT2);
        if( direction == 'z' ) getCellIndicesTY(faceIndex, dataBase, posCellIndexN, negCellIndexN, posCellIndexT2, negCellIndexT2);

        computeGradT( dataBase, parameters, posCellIndexT2, negCellIndexT2, facePrim, gradT2 );
    }
}






#endif