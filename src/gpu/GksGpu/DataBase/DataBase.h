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
//! \file DataBase.h
//! \ingroup DataBase
//! \author Stephan Lenz
//=======================================================================================
#ifndef DataBase_H
#define DataBase_H

#include <memory>
#include <string>
#include <vector>
#include <array>

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"
#include "Core/ArrayTypes.h"

#include "VirtualFluidsDefinitions.h"

#include "Definitions/AccumulatorDataType.h"

#include "CellProperties/CellProperties.cuh"

#include "GksGpu_export.h"

class  GksMeshAdapter;

struct BoundaryCondition;
class  DataBaseAllocator;
struct DataBase;
struct PerLevelCounts;
struct DataBaseStruct;

//! \brief Class for memory management.
//!
//! This class holds and manages the memory of the simulation
//! For this purpose it holds vectors and pointer for host and device/host data, respectively.
struct GKSGPU_EXPORT DataBase : public std::enable_shared_from_this<DataBase>
{
    //////////////////////////////////////////////////////////////////////////
    // Management
    //////////////////////////////////////////////////////////////////////////

    //! shared pointer to an \ref DataBaseAllocator that is used for all memory operations
    SPtr<DataBaseAllocator> myAllocator;

    //! vector of shared pointers to boundary conditions that are used in the simulation
    std::vector< SPtr<BoundaryCondition> > boundaryConditions;

    //////////////////////////////////////////////////////////////////////////
    // Sizes
    //////////////////////////////////////////////////////////////////////////

    uint numberOfNodes;

    uint numberOfCells;

    uint numberOfFaces;

    uint numberOfLevels;

    //! stores number of cells and faces and their start indices per level
    std::vector<PerLevelCounts> perLevelCount;

    //////////////////////////////////////////////////////////////////////////
    // Host only geometry and connectivity
    //////////////////////////////////////////////////////////////////////////

    std::vector<Vec3>   nodeCoordinates;            //!< length = DataBase::numberOfNodes

    std::vector<uint_8> cellToNode;                 //!< Stores node indices. <br>length = DataBase::numberOfCells
    std::vector<uint_4> faceToNode;                 //!< Stores node indices. <br>length = DataBase::numberOfFaces

    std::vector<CellProperties> cellPropertiesHost; //!< length = DataBase::numberOfCells

    //////////////////////////////////////////////////////////////////////////
    // Host/Device geometry and connectivity - READ ONLY
    //////////////////////////////////////////////////////////////////////////

    uint* cellToCell;               //!< Stores cell indices. <br>length = 6 * DataBase::numberOfCells

    uint* faceToCell;               //!< Stores cell indices. <br>length = 2 * DataBase::numberOfFaces

    uint* parentCell;               //!< Stores cell indices. <br>length =     DataBase::numberOfCells

    real* faceCenter;               //!< length = 3 * DataBase::numberOfFaces
    real* cellCenter;               //!< length = 3 * DataBase::numberOfCells

    CellProperties* cellProperties; //!< length =     DataBase::numberOfCells

    char* faceOrientation;          //!< Can be 'x', 'y' or 'z'. <br>length =     DataBase::numberOfFaces; 

    //////////////////////////////////////////////////////////////////////////
    // Host/Device data - READ MODIFY
    //////////////////////////////////////////////////////////////////////////

    real*            data;          //!< Cell averaged flow state data in terms of conserved variables. <br>length = LENGTH_CELL_DATA * DataBase::numberOfCells
    realAccumulator* dataUpdate;    //!< Flux accumulator in terms of conserved variables.              <br>length = LENGTH_CELL_DATA * DataBase::numberOfCells

    real* massFlux;                 //!< Mass flux accumulator.  <br>length = 3 * DataBase::numberOfCells

    real* diffusivity;              //!< Turbulent diffusivity accumulator. <br>length = DataBase::numberOfCells

    int* crashCellIndex;            //!< Index of the crashed cell. It is negative if no cell crashed. <br>length = 1;

    //////////////////////////////////////////////////////////////////////////
    // Host only data
    //////////////////////////////////////////////////////////////////////////

    std::vector<real> dataHost;         //!< Copy of data for initialization and post processing. <br>length = LENGTH_CELL_DATA * DataBase::numberOfCells

    std::vector<real> diffusivityHost;  //!< Copy of diffusivity for post processing. <br>length = DataBase::numberOfCells

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    // Methods
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    //! constructor
    //! initializes the pointers and creates an DataBaseAllocator
    //! \param type string that can be either "GPU" or "CPU", depending on where the computations should be executed
    DataBase( std::string type );

    //! destructor
    //! initiates memory deallocation
    ~DataBase();

    //! copies a mesh from a <b>GksMeshAdapter</b> to the DataStructure and allocates device memory for the simulation
    //! \param adapter   <b>GksMeshAdapter</b> that contains a mesh
    void setMesh( GksMeshAdapter& adapter );

    
    void copyDataHostToDevice();                    //!< copies the flow state data from DataBase::dataHost to DataBase::data
    
    void copyDataDeviceToHost();                    //!< copies the flow state data from DataBase::data to DataBase::dataHost
    
    void copyDataDeviceToHost( real* dataHost );    //!< copies the flow state data from DataBase::data to dataHost, where dataHost may be any array of suitable size
    
    int getCrashCellIndex();                        //!< downloads the crash cell index from the device

    DataBaseStruct toStruct();                      //!< exports all memory pointers and sizes for kernel usage

    //////////////////////////////////////////////////////////////////////////

    uint getCellLevel( uint cellIdx );              //!< \return grid level of cell with index cellIdx
    uint getFaceLevel( uint faceIdx );              //!< \return grid level of face with index faceIdx

    Vec3 getCellCenter( uint cellIdx );             //!< \return cell of cell with index cellIdx

    bool isGhostCell( uint cellIdx );               //!< \return true if cell is ghost cell

    std::string getDeviceType();                    //!< \return "GPU" or "CPU" depending on which DataBaseAllocator is used
};

//! \brief Stores number of cells and faces and their start indices per level
//!
//! Additionally information for refinement and 
struct GKSGPU_EXPORT PerLevelCounts
{
    uint numberOfCells;
    uint startOfCells;

    uint numberOfBulkCells;

    uint numberOfFaces;

    uint numberOfInnerFaces;

    uint numberOfFacesX;
    uint startOfFacesX;

    uint numberOfFacesY;
    uint startOfFacesY;

    uint numberOfFacesZ;
    uint startOfFacesZ;

    uint numberOfCoarseToFine;
    uint startOfCoarseToFine;

    uint numberOfFineToCoarse;
    uint startOfFineToCoarse;
};

#endif