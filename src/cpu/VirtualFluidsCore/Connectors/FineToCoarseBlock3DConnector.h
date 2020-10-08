//! \file FineToCoarseBlock3DConnector.h
//! \brief Base class for connectors that interpolates and sends data from fine level to coarse.  
//! \author Konstantin Kutscher
//! \date 21.05.2015

#ifndef FineToCoarseBlock3DConnector_H
#define FineToCoarseBlock3DConnector_H

#include "TransmitterType.h"
#include "Block3DConnector.h"
#include "D3Q27System.h"
#include "Block3D.h"
#include "LBMKernel.h"
#include "InterpolationProcessor.h"

#include <PointerDefinitions.h>


class Block3D;

//! \class FineToCoarseBlock3DConnector
//! \brief Base class for connectors that interpolates and sends data from fine level to coarse.  
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

class FineToCoarseBlock3DConnector : public Block3DConnector
{
public:
   enum CFconnectorType { Type00, Type10, Type01, Type11 };
public:
   FineToCoarseBlock3DConnector(SPtr<Block3D> block, VectorTransmitterPtr sender, VectorTransmitterPtr receiver, int sendDir, InterpolationProcessorPtr iprocessor, CFconnectorType connType);

   bool isLocalConnector() override;
   bool isRemoteConnector() override;
   
   void sendTransmitterDataSize() override;
   void receiveTransmitterDataSize() override;

   void prepareForSend() override;
   void sendVectors() override;

   void prepareForReceive() override;
   void receiveVectors() override;

   void init() override =0;
   void fillSendVectors() override =0;
   void distributeReceiveVectors() override =0;

   bool isInterpolationConnectorCF() override { return false; }
   bool isInterpolationConnectorFC() override { return true; }

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
   WPtr<Block3D> block; 
   VectorTransmitterPtr sender, receiver;
   InterpolationProcessorPtr iprocessor;
   CFconnectorType connType;
};



#endif 

