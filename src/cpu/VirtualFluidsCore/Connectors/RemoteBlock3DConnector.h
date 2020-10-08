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

   bool isLocalConnector() override;
   bool isRemoteConnector() override;

   void init() override = 0;

   void sendTransmitterDataSize() override;
   void receiveTransmitterDataSize() override;

   void prepareForSend() override;
   void sendVectors() override;

   void prepareForReceive() override;
   void receiveVectors() override;

   void fillSendVectors() override = 0;
   void distributeReceiveVectors() override = 0;

   bool isInterpolationConnectorCF() override { return false; }
   bool isInterpolationConnectorFC() override { return false; }

   double getSendRecieveTime() { return 0; }

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
   VectorTransmitterPtr sender;
   VectorTransmitterPtr receiver;
};


#endif //RemoteBlock3DConnector_H

