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
//! \file DataBaseAllocator.h
//! \ingroup DataBase
//! \author Stephan Lenz
//=======================================================================================
#ifndef DataBaseAllocator_H
#define DataBaseAllocator_H

#include <string>
#include <vector>

#include "Core/DataTypes.h"
#include "Core/PointerDefinitions.h"

#include "VirtualFluidsDefinitions.h"

#include "GksGpu_export.h"

class  GksMeshAdapter;
struct DataBase;
struct BoundaryCondition;

//! Virtual base class for GPU and CPU memory management
class GKSGPU_EXPORT DataBaseAllocator {

public:

    //! \return shared pointer to DataBaseAllocatorGPU or DataBaseAllocatorCPU based on type
    static std::shared_ptr<DataBaseAllocator> create( std::string type );

    //////////////////////////////////////////////////////////////////////////

    //! frees DataBase memory
    virtual void freeMemory( DataBase& dataBase ) = 0;

    //! allocates DataBase memory
    virtual void allocateMemory( SPtr<DataBase> dataBase) = 0;

    //! copies mesh information from GksMeshAdapter to DataBase
    virtual void copyMesh( SPtr<DataBase> dataBase, GksMeshAdapter& adapter ) = 0;

    //! copies the flow state data from DataBase::dataHost to DataBase::data
    virtual void copyDataHostToDevice( SPtr<DataBase> dataBase ) = 0;

    //! copies the flow state data from DataBase::data to dataHost
    virtual void copyDataDeviceToHost( SPtr<DataBase> dataBase, real* hostData ) = 0;

    //! \return the crash cell index from the device
    virtual int  getCrashCellIndex( SPtr<DataBase> dataBase ) = 0;                      

    //////////////////////////////////////////////////////////////////////////

    //! frees boundary condition memory
    virtual void freeMemory( BoundaryCondition& boundaryCondition ) = 0;

    //! allocates boundary condition memory and copies data to device
    virtual void allocateMemory( SPtr<BoundaryCondition> boundaryCondition, std::vector<uint> ghostCells, std::vector<uint> domainCells, std::vector<uint> secondCells ) = 0;

    //////////////////////////////////////////////////////////////////////////

    ~DataBaseAllocator();

    // \return the device type as string
    virtual std::string getDeviceType() = 0;

protected:

    DataBaseAllocator();
    DataBaseAllocator( const DataBaseAllocator& orig );

};


#endif