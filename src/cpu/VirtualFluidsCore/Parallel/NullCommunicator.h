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
//! \file NullCommunicator.h
//! \ingroup Parallel
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef NullCommunicator_H
#define NullCommunicator_H

#include "Communicator.h"

#include <PointerDefinitions.h>

//! \brief A class implements Communicator for shared memory.
//! \details NullCommunicator is only a place-holder. It is only one process in shared memory.
class NullCommunicator : public Communicator
{
public:
   NullCommunicator();
   ~NullCommunicator() override;
   int getBundleID() override;
   int getNumberOfBundles() override;
   int getProcessID() override;
   int getNumberOfProcesses() override;
   void* getNativeCommunicator() override;
   int getRoot() override;
   int getBundleRoot() override;
   int getProcessRoot() override;
   std::vector<std::string> gather(const std::string& str) override;
   std::vector<double> gatherDoubles(std::vector<double>& values); 
   void allGatherInts(std::vector<int>& svalues, std::vector<int>& rvalues);
   void sendSerializedObject(std::stringstream& ss, int target) override;
   void receiveSerializedObject(std::stringstream& ss, int source) override;
protected:
private:
};

#endif
