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
//! \file CellUpdate.h
//! \ingroup CellUpdate
//! \author Stephan Lenz
//=======================================================================================
#ifndef  CellUpdate_H
#define  CellUpdate_H

#include "VirtualFluidsDefinitions.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"
#include "Parameters/Parameters.h"

#include "GksGpu_export.h"

//! \brief Performs the cell update at the end of a time step
//! 
//! This includes adding the fluxes and the source that are caused by volume forces.
//! Further, it checks if the simulation crashed and stores the cell index of the crashed cell.
class GKSGPU_EXPORT CellUpdate
{
public:

    //! \brief executes the cell update
    //!
    //! This function is only the entry point of the cell update. The execution is 
    //! given to CPU or GPU based on the type of DataBaseAllocator by the \ref runKernel
    //! function which calls the \ref cellUpdateFunction. 
    //!
    //! \param dataBase    shared pointer to a \ref DataBase
    //! \param parameters  \ref Parameters object
    //! \param level       grid level to which the cell update should be applied
    static void run( SPtr<DataBase> dataBase, 
                     Parameters parameters, 
                     uint level );
};

#endif
