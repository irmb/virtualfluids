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
//! \file WriteBlocksCoProcessor.h
//! \ingroup CoProcessors
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef WriteBlocksCoProcessor_H_
#define WriteBlocksCoProcessor_H_

#include <PointerDefinitions.h>
#include <string>

#include "CoProcessor.h"

class Communicator;
class Grid3D;
class UbScheduler;
class WbWriter;

//! \class WriteBlocksCoProcessor
//! \brief A class writes a block grid to a VTK-file
class WriteBlocksCoProcessor: public CoProcessor 
{
public:
   //! \brief Construct WriteBlocksCoProcessor object.
   //! \pre The Grid3D and UbScheduler objects must exist.
   //! \param grid is observable Grid3D object
   //! \param s is UbScheduler object for scheduling of observer
   //! \param path is path of folder for output
   //! \param writer is WbWriter object
   //! \param comm is Communicator object
   WriteBlocksCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path, WbWriter* const writer, SPtr<Communicator> comm);
   virtual ~WriteBlocksCoProcessor();

   void process(double step) override;

protected:
   //! Collect data for VTK-file
   //! \param step is a time step
   void collectData(double step);

   std::string path;
   WbWriter* writer;
   SPtr<Communicator>  comm;
};


#endif 
