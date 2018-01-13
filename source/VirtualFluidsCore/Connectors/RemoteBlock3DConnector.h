#ifndef RemoteBlock3DConnector_H
#define RemoteBlock3DConnector_H

#include <vector>

#include "TransmitterType.h"
#include "Block3DConnector.h"
#include "D3Q27System.h"
#include "Block3D.h"
#include "LBMKernel.h"
#include "EsoTwistD3Q27System.h"


//daten werden in einen vector (dieser befindet sich im transmitter) kopiert
//der vector wird via transmitter uebertragen
//transmitter kann ein lokal, MPI, RCG, CTL oder was auch immer fuer ein
//transmitter sein, der von Transmitter abgeleitet ist ;-)
class RemoteBlock3DConnector : public Block3DConnector
{
public:
   RemoteBlock3DConnector(SPtr<Block3D> block
      , VectorTransmitterPtr sender
      , VectorTransmitterPtr receiver
      , int sendDir);

   bool isLocalConnector();
   bool isRemoteConnector();

   virtual void init() = 0;

   void sendTransmitterDataSize();
   void receiveTransmitterDataSize();

   void prepareForSend();
   void sendVectors();

   void prepareForReceive();
   void receiveVectors();

   virtual void fillSendVectors() = 0;
   virtual void distributeReceiveVectors() = 0;

   bool isInterpolationConnectorCF() { return false; }
   bool isInterpolationConnectorFC() { return false; }

   double getSendRecieveTime() { return 0; }

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
   WPtr<Block3D> block; 
   VectorTransmitterPtr sender;
   VectorTransmitterPtr receiver;
};


#endif //RemoteBlock3DConnector_H

