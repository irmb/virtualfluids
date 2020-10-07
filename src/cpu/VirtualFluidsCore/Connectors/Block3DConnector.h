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
//! \file Block3DConnector.h
//! \ingroup Connectors
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef BLOCK2DCONNECTOR_H
#define BLOCK2DCONNECTOR_H

#include <vector>
#include <string>

#include <basics/utilities/UbTuple.h>

#include <PointerDefinitions.h>

//! \brief   Abstract class of connectors  
//! \details Connector send and receive full distributions between two blocks in shared memory.
class Block3DConnector
{
public:
   Block3DConnector() = default;
   Block3DConnector(const int& sendDir) : sendDir(sendDir) {}
   virtual ~Block3DConnector() = default;
   //!Iniitializes connector
   virtual void init()=0;
   //!Synchronizes the send-buffer length
   virtual void sendTransmitterDataSize()=0;  
   //!Synchronizes the receive-buffer length
   virtual void receiveTransmitterDataSize()=0;

   //Send (should be called in given order!!!)
   virtual void prepareForSend()=0;
   virtual void fillSendVectors()=0;
   virtual void sendVectors()=0;
   
   //Receive (should be called in given order!!!)
   virtual void prepareForReceive()=0;
   virtual void receiveVectors()=0;
   virtual void distributeReceiveVectors()=0;

   //info section
   virtual bool isLocalConnector()  = 0;
   virtual bool isRemoteConnector() = 0;
   virtual bool isInterpolationConnectorCF() = 0;
   virtual bool isInterpolationConnectorFC() = 0;

   //grid refinement
   virtual int getSendDir() const { return sendDir; } 

protected:
   int  sendDir{-1};
};

#endif //BLOCK3DCONNECTOR_H
