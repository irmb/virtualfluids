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
//! \file NullCommunicator.cpp
//! \ingroup Parallel
//! \author Konstantin Kutscher
//=======================================================================================

#include "NullCommunicator.h"


NullCommunicator::NullCommunicator()
{
}
//////////////////////////////////////////////////////////////////////////
NullCommunicator::~NullCommunicator()
{
}
//////////////////////////////////////////////////////////////////////////
int NullCommunicator::getBundleID() 
{
   return 0;
}
//////////////////////////////////////////////////////////////////////////
int NullCommunicator::getNumberOfBundles() 
{
   return 0;
}
//////////////////////////////////////////////////////////////////////////
int NullCommunicator::getProcessID() 
{
   return 0;
}
//////////////////////////////////////////////////////////////////////////
int NullCommunicator::getNumberOfProcesses()
{
   return 0;
}
//////////////////////////////////////////////////////////////////////////
void* NullCommunicator::getNativeCommunicator()
{
   return NULL;
}
//////////////////////////////////////////////////////////////////////////
int NullCommunicator::getRoot() 
{
   return 0;
}
//////////////////////////////////////////////////////////////////////////
int NullCommunicator::getBundleRoot() 
{
   return 0;
}
//////////////////////////////////////////////////////////////////////////
int NullCommunicator::getProcessRoot() 
{
   return 0;
}
//////////////////////////////////////////////////////////////////////////
std::vector<std::string> NullCommunicator::gather(const std::string& str)
{
   return std::vector<std::string>();
}
//////////////////////////////////////////////////////////////////////////
std::vector<double> NullCommunicator::gatherDoubles(std::vector<double>& values) 
{
   return std::vector<double>();
}
//////////////////////////////////////////////////////////////////////////
void NullCommunicator::allGatherInts(std::vector<int>& svalues, std::vector<int>& rvalues)
{

}
//////////////////////////////////////////////////////////////////////////
void NullCommunicator::sendSerializedObject(std::stringstream& ss, int target) 
{

}
//////////////////////////////////////////////////////////////////////////
void NullCommunicator::receiveSerializedObject(std::stringstream& ss, int source) 
{

}
