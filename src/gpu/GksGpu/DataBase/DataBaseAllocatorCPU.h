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
//! \file DataBaseAllocatorCPU.h
//! \ingroup DataBase
//! \author Stephan Lenz
//=======================================================================================
#ifndef DataBaseAllocatorCPU_H
#define DatabaseAllocatorCPU_H

#include "Core/DataTypes.h"
#include "Core/PointerDefinitions.h"

#include "DataBaseAllocator.h"

#include "VirtualFluidsDefinitions.h"

//! DataBaseAllocator for CPU memory
class GKSGPU_EXPORT DataBaseAllocatorCPU : public DataBaseAllocator {

public:

    //! frees DataBase memory
    virtual void freeMemory( DataBase& dataBase ) override;

    //! allocates DataBase memory
    virtual void allocateMemory( SPtr<DataBase> dataBase ) override;

    //! copies mesh information from GksMeshAdapter to DataBase
    virtual void copyMesh( SPtr<DataBase> dataBase, GksMeshAdapter& adapter ) override;

    //! copies the flow state data from DataBase::dataHost to DataBase::data
    virtual void copyDataHostToDevice( SPtr<DataBase> dataBase ) override;
    
    //! copies the flow state data from DataBase::data to dataHost
    virtual void copyDataDeviceToHost( SPtr<DataBase> dataBase, real* dataHost ) override;

    //! \return the crash cell index from the device
    virtual int  getCrashCellIndex( SPtr<DataBase> dataBase ) override;

    //////////////////////////////////////////////////////////////////////////

    //! frees boundary condition memory
    virtual void freeMemory( BoundaryCondition& boundaryCondition ) override;

    //! allocates boundary condition memory and copies data to device
    virtual void allocateMemory( SPtr<BoundaryCondition> boundaryCondition, std::vector<uint> ghostCells, std::vector<uint> domainCells, std::vector<uint> secondCells ) override;

    //////////////////////////////////////////////////////////////////////////

    //! \return "CPU"
    virtual std::string getDeviceType() override;
};


#endif