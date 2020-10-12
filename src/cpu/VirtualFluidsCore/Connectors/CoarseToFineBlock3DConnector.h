//! \file CoarseToFineBlock3DConnector.h
//! \brief Base class for connectors that interpolates and sends data from coarse level to fine.  
//! \author Konstantin Kutscher
//! \date 18.05.2015

#ifndef CoarseToFineBlock3DConnector_H
#define CoarseToFineBlock3DConnector_H

#include "TransmitterType.h"
#include "Block3DConnector.h"
#include "D3Q27System.h"
#include "Block3D.h"
#include "InterpolationProcessor.h"
#include <PointerDefinitions.h>


class Block3D;

//! \class CoarseToFineBlock3DConnector
//! \brief Base class for connectors that interpolates and sends data from coarse level to fine. 
//! \details The data is copied in a vector (this is located in the transmitter). 
//! The vector is transmitted via transmitter. 
//! The transmitter can be a local, MPI, RCG, CTL or whatever 
//! which a transmitter that is derived from transmitter base class.
//!
//! four fine blocks inside a coarse block:
//!
//! |    |    |
//! |:--:|:---| 
//! | 01 | 11 | 
//! | 00 | 10 | 
//!
//! send direction:    
//!
//! |E<->W   |  N<->S  |  T<->B |
//! |--------|---------|--------|
//! |  x3    |   x3    |    x2  |
//! |  ^     |   ^     |    ^   |
//! |  +->x2 |  +->x1  |   +->x1|


class CoarseToFineBlock3DConnector : public Block3DConnector
{
public:
   CoarseToFineBlock3DConnector(SPtr<Block3D> block,
      VectorTransmitterPtr sender00, VectorTransmitterPtr receiver00,
      VectorTransmitterPtr sender01, VectorTransmitterPtr receiver01,
      VectorTransmitterPtr sender10, VectorTransmitterPtr receiver10,
      VectorTransmitterPtr sender11, VectorTransmitterPtr receiver11,
      int sendDir, InterpolationProcessorPtr iprocessor);

   ~CoarseToFineBlock3DConnector() override = default;

   bool isLocalConnector() override;
   bool isRemoteConnector() override;

   void init() override =0;

   void sendTransmitterDataSize() override;
   void receiveTransmitterDataSize() override;

   void prepareForSend() override;
   void sendVectors() override;

   void prepareForReceive() override;
   void receiveVectors() override;

   void fillSendVectors() override =0;
   void distributeReceiveVectors() override = 0;

   bool isInterpolationConnectorCF() override { return true; }
   bool isInterpolationConnectorFC() override { return false; }

   void prepareForSendX1() override {}
   void prepareForSendX2() override {}
   void prepareForSendX3() override {}

   void sendVectorsX1() override {}
   void sendVectorsX2() override {}
   void sendVectorsX3() override {}

   void prepareForReceiveX1() override {}
   void prepareForReceiveX2() override {}
   void prepareForReceiveX3() override {}

   void receiveVectorsX1() override {}
   void receiveVectorsX2() override {}
   void receiveVectorsX3() override {}

protected:
   WPtr<Block3D> block; //dieser nvd sendet daten und die empfangenen werden diesem nvd zugeordnet
   VectorTransmitterPtr sender00, receiver00,
                        sender01, receiver01,
                        sender10, receiver10,
                        sender11, receiver11;

   InterpolationProcessorPtr iprocessor;

};




#endif 

