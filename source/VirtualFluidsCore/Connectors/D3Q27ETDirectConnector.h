#ifndef D3Q27DIRECTCONNECTOR_H
#define D3Q27DIRECTCONNECTOR_H

#include <boost/weak_ptr.hpp>

#include "Block3DConnector.h"
#include "Block3D.h"

class D3Q27ETDirectConnector : public Block3DConnector
{
public:
   D3Q27ETDirectConnector(Block3DPtr from, Block3DPtr to, const int& sendDir) 
      :  Block3DConnector(sendDir)
        , from(from)
        , to(to)
   {

   }
   
   void sendTransmitterDataSize()    { }  
   void receiveTransmitterDataSize() { }
   void init()                       { }
   void prepareForReceive()          { }
   void prepareForSend()             { }
   void fillSendVectors()            { }
   void sendVectors();//                { }
   void receiveVectors()             { }
   
   void distributeReceiveVectors()   { }

   bool isLocalConnector()  { return true;  }
   bool isRemoteConnector() { return false; }
   bool isInterpolationConnectorCF() { return false; }
   bool isInterpolationConnectorFC() { return false; }

   double getSendRecieveTime();

   void prepareForSendX1() {}
   void prepareForSendX2() {}
   void prepareForSendX3() {}

   void sendVectorsX1(){}
   void sendVectorsX2(){}
   void sendVectorsX3(){}

   void prepareForReceiveX1() {}
   void prepareForReceiveX2() {}
   void prepareForReceiveX3() {}

   void receiveVectorsX1() {}
   void receiveVectorsX2() {}
   void receiveVectorsX3() {}

protected:
   boost::weak_ptr<Block3D> from;
   boost::weak_ptr<Block3D> to;
};

#endif //D3Q27DIRECTCONNECTOR_H

