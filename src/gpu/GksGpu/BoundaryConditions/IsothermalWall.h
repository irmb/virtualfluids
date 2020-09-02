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
//! \file IsothermalWall.h
//! \ingroup BoundaryCondition
//! \author Stephan Lenz
//=======================================================================================
#ifndef IsothermalWall_CUH
#define IsothermalWall_CUH

#include <memory>

#include "VirtualFluidsDefinitions.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "BoundaryConditions/BoundaryCondition.h"

#include "GksGpu_export.h"

//! Simple struct that carries variables to the boundary condition kernel
struct IsothermalWallStruct
{
    uint  numberOfCells;    //!< number of boundary cells
                            
    uint* ghostCells;       //!< pointer to device or host memory with ghost cell indices
    uint* domainCells;      //!< pointer to device or host memory with domain cell indices
    uint* secondCells;      //!< pointer to device or host memory with second cell indices

    Vec3 velocity;          //!< wall velocity
    real lambda;            //!< wall temperature in terms of lambda
    real S_1;               //!< wall value for passive scalar
    real S_2;               //!< wall value for passive scalar

    bool useSecondCells;    //!< true if linear pressure extrapolation should be used
};

//! Models a wall that is kept at a constant temperature
struct GKSGPU_EXPORT IsothermalWall : public BoundaryCondition //, public IsothermalWallStruct
{
    Vec3 velocity;          //!< wall velocity
    real lambda;            //!< wall temperature in terms of lambda
    real S_1;               //!< wall value for passive scalar
    real S_2;               //!< wall value for passive scalar

    bool useSecondCells;    //!< true if linear pressure extrapolation should be used

    //! constructor
    //! \param dataBase        shared pointer to a \ref DataBase
    //! \param velocity        wall velocity
    //! \param lambda          wall temperature in terms of lambda
    //! \param useSecondCells  true if linear pressure extrapolation should be used
    //! \param S_1             wall value for passive scalar
    //! \param S_2             wall value for passive scalar
    IsothermalWall( SPtr<DataBase> dataBase, Vec3 velocity, real lambda, bool useSecondCells, real S_1 = 0.0, real S_2 = 0.0 );

    //! \return true
    virtual bool isWall() override;
    
    //! \return true
    virtual bool secondCellsNeeded() override;

    //! executes the boundary condition
    //! \param dataBase    shared pointer to a \ref DataBase
    //! \param parameters  \ref Parameters object
    //! \param level       grid level to which the boundary condition should be applied
    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    // \return struct with variables for boundary condition kernel
    IsothermalWallStruct toStruct()
    {
        IsothermalWallStruct boundaryCondition;

        boundaryCondition.numberOfCells = this->numberOfCells;

        boundaryCondition.ghostCells      = this->ghostCells;
        boundaryCondition.domainCells     = this->domainCells;
        boundaryCondition.secondCells     = this->secondCells;

        boundaryCondition.velocity        = this->velocity;
        boundaryCondition.lambda          = this->lambda;
        boundaryCondition.S_1             = this->S_1;
        boundaryCondition.S_2             = this->S_2;

        boundaryCondition.useSecondCells  = this->useSecondCells;

        return boundaryCondition;
    }
};

#endif
