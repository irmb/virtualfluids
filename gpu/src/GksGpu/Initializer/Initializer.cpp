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
//! \file Initializer.cpp
//! \ingroup Initializer
//! \author Stephan Lenz
//=======================================================================================
#include "Initializer.h"

#include <sstream>
#define _USE_MATH_DEFINES
#include <math.h>

#include "Core/PointerDefinitions.h"
#include "Core/RealConstants.h"

#include "DataBase/DataBaseStruct.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

void Initializer::interpret(SPtr<DataBase> dataBase, std::function<ConservedVariables(Vec3)> initialCondition)
{
    for( uint cellIdx = 0; cellIdx < dataBase->numberOfCells; cellIdx++ ){

        Vec3 cellCenter = dataBase->getCellCenter( cellIdx );

        ConservedVariables cellCons = initialCondition(cellCenter);

        dataBase->dataHost[ RHO__(cellIdx, dataBase->numberOfCells) ] = cellCons.rho ;
        dataBase->dataHost[ RHO_U(cellIdx, dataBase->numberOfCells) ] = cellCons.rhoU;
        dataBase->dataHost[ RHO_V(cellIdx, dataBase->numberOfCells) ] = cellCons.rhoV;
        dataBase->dataHost[ RHO_W(cellIdx, dataBase->numberOfCells) ] = cellCons.rhoW;
        dataBase->dataHost[ RHO_E(cellIdx, dataBase->numberOfCells) ] = cellCons.rhoE;
    #ifdef USE_PASSIVE_SCALAR
	    dataBase->dataHost[ RHO_S_1(cellIdx, dataBase->numberOfCells) ] = cellCons.rhoS_1;
	    dataBase->dataHost[ RHO_S_2(cellIdx, dataBase->numberOfCells) ] = cellCons.rhoS_2;
    #endif // USE_PASSIVE_SCALAR
    }

    return;
}
