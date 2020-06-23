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
//! \file MemoryAccessPattern.h
//! \ingroup Definitions
//! \author Stephan Lenz
//=======================================================================================
#ifndef DataDefinitions_H
#define DataDefinitions_H

#include "PassiveScalar.h"

//! \file MemoryAccessPattern.h 
//! The Macros defined in this file are used to switch between AOS and SOA data structures.
//! The usage is:
//! \code value = array[ MACRO(index, arrayLength) ] \endcode
//! or 
//! \code value = array[ MACRO(index, subIndex, arrayLength) ] \endcode

//! this defines the memory access pattern and can either be AOS or SOA
#define SOA

//////////////////////////////////////////////////////////////////////////

#define LENGTH_VECTOR       3           //!< number of components of a physical vector

#ifdef USE_PASSIVE_SCALAR
    #define LENGTH_CELL_DATA    7       //!< number of components of the flow state. <br>This accounts for density, three momentum components, energy and two passive scalars
#else
    #define LENGTH_CELL_DATA    5       //!< number of components of the flow state. <br>This accounts for density, three momentum components and energy
#endif

#define LENGTH_CELL_TO_CELL 6           //!< number of neighbor cells

#define LENGTH_FACE_TO_CELL 2           //!< number of cells connected by a face

#define LENGTH_FINE_TO_COARSE 9         //!< number of cells in the fine to coarse interpolation

#define LENGTH_COARSE_TO_FINE 9         //!< number of cells in the coarse to fine interpolation

//////////////////////////////////////////////////////////////////////////

#ifdef SOA

#define VEC_X(vecIdx, numberOfVectors)  ( 0 * numberOfVectors + vecIdx )        //!< computes position of x component in vector array
#define VEC_Y(vecIdx, numberOfVectors)  ( 1 * numberOfVectors + vecIdx )        //!< computes position of y component in vector array
#define VEC_Z(vecIdx, numberOfVectors)  ( 2 * numberOfVectors + vecIdx )        //!< computes position of z component in vector array
                                                           
#define RHO__( cellIdx, numberOfCells ) ( 0 * numberOfCells   + cellIdx )       //!< computes position of rho  component in DataBase::data array, or similar arrays
#define RHO_U( cellIdx, numberOfCells ) ( 1 * numberOfCells   + cellIdx )       //!< computes position of rhoU component in DataBase::data array, or similar arrays
#define RHO_V( cellIdx, numberOfCells ) ( 2 * numberOfCells   + cellIdx )       //!< computes position of rhoV component in DataBase::data array, or similar arrays
#define RHO_W( cellIdx, numberOfCells ) ( 3 * numberOfCells   + cellIdx )       //!< computes position of rhoW component in DataBase::data array, or similar arrays
#define RHO_E( cellIdx, numberOfCells ) ( 4 * numberOfCells   + cellIdx )       //!< computes position of rhoE component in DataBase::data array, or similar arrays

#ifdef USE_PASSIVE_SCALAR
    #define RHO_S_1( cellIdx, numberOfCells ) ( 5 * numberOfCells   + cellIdx ) //!< computes position of rho_S_1 component in data array
    #define RHO_S_2( cellIdx, numberOfCells ) ( 6 * numberOfCells   + cellIdx ) //!< computes position of rho_S_2 component in data array
#endif // USE_PASSIVE_SCALAR

#define CELL_TO_CELL( cellIdx, neighborIdx, numberOfCells ) ( neighborIdx * numberOfCells + cellIdx )   //!< computes position of index of neighbor cell in  DataBase::cellToCell array

#define NEG_CELL( faceIdx, numberOfFaces ) (                 faceIdx )          //!< computes position of positive cell index in faceToCell array
#define POS_CELL( faceIdx, numberOfFaces ) ( numberOfFaces + faceIdx )          //!< computes position of negative cell index in faceToCell array

#define FINE_TO_COARSE( idx, cellIdx, number ) ( cellIdx * number + idx )       //!< computes position of fine to coarse interpolation cell index in DataBase::fineToCoarse
#define COARSE_TO_FINE( idx, cellIdx, number ) ( cellIdx * number + idx )       //!< computes position of coarse to fine interpolation cell index in DataBase::coarseToFine

#endif

//////////////////////////////////////////////////////////////////////////

#ifdef AOS

#define VEC_X(vecIdx, numberOfVectors)  ( vecIdx * LENGTH_VECTOR     )          //!< computes position of x component in vector array
#define VEC_Y(vecIdx, numberOfVectors)  ( vecIdx * LENGTH_VECTOR + 1 )          //!< computes position of y component in vector array
#define VEC_Z(vecIdx, numberOfVectors)  ( vecIdx * LENGTH_VECTOR + 2 )          //!< computes position of z component in vector array

#define RHO__( cellIdx, numberOfCells ) ( cellIdx * LENGTH_CELL_DATA     )      //!< computes position of rho  component in DataBase::data array, or similar arrays
#define RHO_U( cellIdx, numberOfCells ) ( cellIdx * LENGTH_CELL_DATA + 1 )      //!< computes position of rhoU component in DataBase::data array, or similar arrays
#define RHO_V( cellIdx, numberOfCells ) ( cellIdx * LENGTH_CELL_DATA + 2 )      //!< computes position of rhoV component in DataBase::data array, or similar arrays
#define RHO_W( cellIdx, numberOfCells ) ( cellIdx * LENGTH_CELL_DATA + 3 )      //!< computes position of rhoW component in DataBase::data array, or similar arrays
#define RHO_E( cellIdx, numberOfCells ) ( cellIdx * LENGTH_CELL_DATA + 4 )      //!< computes position of rhoE component in DataBase::data array, or similar arrays

#ifdef USE_PASSIVE_SCALAR
    #define RHO_S_1( cellIdx, numberOfCells ) ( cellIdx * LENGTH_CELL_DATA + 5 )//!< computes position of rho_S_1 component in data array
    #define RHO_S_2( cellIdx, numberOfCells ) ( cellIdx * LENGTH_CELL_DATA + 6 )//!< computes position of rho_S_2 component in data array
#endif // USE_PASSIVE_SCALAR
                                                                         
#define CELL_TO_CELL( cellIdx, neighborIdx, numberOfCells ) ( cellIdx * LENGTH_CELL_TO_CELL + neighborIdx )//!< computes position of index of neighbor cell in  DataBase::cellToCell array

#define NEG_CELL( faceIdx, numberOfFaces ) ( faceIdx * LENGTH_FACE_TO_CELL     )            //!< computes position of positive cell index in faceToCell array
#define POS_CELL( faceIdx, numberOfFaces ) ( faceIdx * LENGTH_FACE_TO_CELL + 1 )            //!< computes position of negative cell index in faceToCell array

#define FINE_TO_COARSE( idx, cellIdx, number ) ( cellIdx * LENGTH_FINE_TO_COARSE + idx )    //!< Computes position of fine to coarse interpolation cell index in DataBase::fineToCoarse
#define COARSE_TO_FINE( idx, cellIdx, number ) ( cellIdx * LENGTH_COARSE_TO_FINE + idx )    //!< Computes position of coarse to fine interpolation cell index in DataBase::coarseToFine

#endif

//////////////////////////////////////////////////////////////////////////

#endif
