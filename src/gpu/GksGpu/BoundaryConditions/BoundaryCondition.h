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
//! \file BoundaryCondition.h
//! \ingroup BoundaryCondition
//! \author Stephan Lenz
//=======================================================================================
#ifndef BoundaryCondition_H
#define BoundaryCondition_H

#include <functional>

#include <memory>
#include <vector>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "Parameters/Parameters.h"

#include "GksGpu_export.h"

class  GksMeshAdapter;
class  DataBaseAllocator;
struct DataBase;
struct BoundaryConditionStruct;

//! Simple struct that carries variables to the boundary condition kernel
struct BoundaryConditionStruct
{
    uint  numberOfCells;    //!< number of boundary cells

    uint* ghostCells;       //!< pointer to device or host memory with ghost cell indices
    uint* domainCells;      //!< pointer to device or host memory with domain cell indices
    uint* secondCells;      //!< pointer to device or host memory with second cell indices
};

//! \brief purely virtual base class from which specific boundary conditions derive
//! 
//! It implements basic boundary condition functionality, such as finding boundary cells
struct GKSGPU_EXPORT BoundaryCondition : virtual public BoundaryConditionStruct, public std::enable_shared_from_this<BoundaryCondition>
{
    SPtr<DataBaseAllocator> myAllocator;        //!< shared pointer to an \ref DataBaseAllocator, used for memory management

    std::vector<uint> numberOfCellsPerLevel;    
    std::vector<uint> startOfCellsPerLevel;     

    //! constructor that initializes pointers
    //! \param dataBase   shared pointer to a \ref DataBase
    BoundaryCondition( SPtr<DataBase> dataBase );

    //! destructor that initiates memory deallocation
    ~BoundaryCondition();

    //! \brief searches for boundary cells and allocates memory
    //! \param adapter          reference to a <b>GksMeshAdapter</b>, in which the boundary cells are searched
    //! \param allowGhostCells  when true, ghost cells are allowed as domain cells is they are outside the boundary region
    //! \param boundaryFinder   a function that defines the boundary region based on the ghost cell center
    virtual void findBoundaryCells( GksMeshAdapter& adapter, bool allowGhostCells, std::function<bool(Vec3)> boundaryFinder);

    //! \return true if boundary condition is impermeable -> mass flux zero
    virtual bool isWall() = 0;

    //! \return true if boundary condition applies fluxes
    virtual bool isFluxBC();
    
    //! \return true if boundary condition is insulated -> heat flux zero
    virtual bool isInsulated();

    //! \return true is second cells need to be searched and stored
    virtual bool secondCellsNeeded();

    //! executes the boundary condition
    //! \param dataBase    shared pointer to a \ref DataBase
    //! \param parameters  \ref Parameters object
    //! \param level       grid level to which the boundary condition should be applied
    virtual void runBoundaryConditionKernel( const SPtr<DataBase> dataBase,
                                             const Parameters parameters, 
                                             const uint level ) = 0;

};

#endif
