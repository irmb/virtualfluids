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

#include <memory>


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
   FineToCoarseBlock3DConnector(Block3DPtr block, VectorTransmitterPtr sender, VectorTransmitterPtr receiver, int sendDir, InterpolationProcessorPtr iprocessor, CFconnectorType connType);

   bool isLocalConnector();
   bool isRemoteConnector();
   
   void sendTransmitterDataSize();
   void receiveTransmitterDataSize();

   void prepareForSend();
   void sendVectors();

   void prepareForReceive();
   void receiveVectors();

   virtual void init()=0;
   virtual void fillSendVectors()=0;
   virtual void distributeReceiveVectors()=0;

   bool isInterpolationConnectorCF() { return false; }
   bool isInterpolationConnectorFC() { return true; }

   void prepareForSendX1() {}
   void prepareForSendX2() {}
   void prepareForSendX3() {}

   void sendVectorsX1() {}
   void sendVectorsX2() {}
   void sendVectorsX3() {}

   void prepareForReceiveX1() {}
   void prepareForReceiveX2() {}
   void prepareForReceiveX3() {}

   void receiveVectorsX1() {}
   void receiveVectorsX2() {}
   void receiveVectorsX3() {}

protected:
   std::weak_ptr<Block3D> block; 
   VectorTransmitterPtr sender, receiver;
   InterpolationProcessorPtr iprocessor;
   CFconnectorType connType;
};



#endif 

