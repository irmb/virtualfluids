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
//! \file Initializer.h
//! \ingroup Initializer
//! \author Stephan Lenz
//=======================================================================================
#ifndef  Initializer_H
#define  Initializer_H

#include <string>
#include <memory>
#include <functional>

#include "VirtualFluidsDefinitions.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "DataBase/DataBase.h"
#include "FlowStateData/FlowStateData.cuh"

#include "GksGpu_export.h"

//! \brief Class for initial conditions and data initialization
class GKSGPU_EXPORT Initializer
{
public:

    //! \brief interprets an initial condition and set initial value in \ref DataBase::dataHost
    //!
    //! \param dataBase            shared pointer to \ref DataBase
    //! \param initialCondition    function that returns a \ref ConservedVariables object for every location
    static void interpret( SPtr<DataBase> dataBase, std::function<ConservedVariables(Vec3)> initialCondition );

    //! \brief initializes the \ref DataBase::dataUpdate array with zero
    static void initializeDataUpdate( SPtr<DataBase> dataBase );
};

#endif
