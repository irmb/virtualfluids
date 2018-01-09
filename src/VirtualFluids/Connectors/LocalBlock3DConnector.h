#ifndef LocalBlock3DConnector_H
#define LocalBlock3DConnector_H



#include "Block3DConnector.h"
#include "Block3D.h"

class LocalBlock3DConnector : public Block3DConnector
{
public:
   LocalBlock3DConnector(Block3DPtr from, Block3DPtr to, int sendDir)
      : Block3DConnector(sendDir)
      , from(from)
      , to(to)
   {

   }
   virtual ~LocalBlock3DConnector() {}
   void sendTransmitterDataSize() {}
   void receiveTransmitterDataSize() {}
   virtual void init() = 0;
   void prepareForReceive() {}
   void prepareForSend() {}
   void fillSendVectors() {}
   virtual void sendVectors()=0;
   void receiveVectors() {}

   void distributeReceiveVectors() {}

   bool isLocalConnector() { return true; }
   bool isRemoteConnector() { return false; }
   bool isInterpolationConnectorCF() { return false; }
   bool isInterpolationConnectorFC() { return false; }

   double getSendRecieveTime();

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
   std::weak_ptr<Block3D> from;
   std::weak_ptr<Block3D> to;
};

#endif //LocalBlock3DConnector_H
