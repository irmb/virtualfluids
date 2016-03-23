#include "RemoteBlock3DConnector.h"

//////////////////////////////////////////////////////////////////////////
RemoteBlock3DConnector::RemoteBlock3DConnector(Block3DPtr block
   , VectorTransmitterPtr sender
   , VectorTransmitterPtr receiver
   , int sendDir)
   : Block3DConnector(sendDir)
   , block(block)
   , sender(sender)
   , receiver(receiver)
{
   if (!block || !sender || !receiver)
      UB_THROW(UbException(UB_EXARGS, "sender or receiver == NULL!!"));
}
//////////////////////////////////////////////////////////////////////////

bool RemoteBlock3DConnector::isLocalConnector()
{
   return !this->isRemoteConnector();
}
//////////////////////////////////////////////////////////////////////////

bool RemoteBlock3DConnector::isRemoteConnector()
{
   return ((sender && sender->isRemoteTransmitter())
      || (receiver && receiver->isRemoteTransmitter()));
}
//////////////////////////////////////////////////////////////////////////

void RemoteBlock3DConnector::sendTransmitterDataSize()
{
   assert(sender  !=NULL); sender->sendDataSize();
}
//////////////////////////////////////////////////////////////////////////

void RemoteBlock3DConnector::receiveTransmitterDataSize()
{
   assert(receiver!=NULL); receiver->receiveDataSize();
}
//////////////////////////////////////////////////////////////////////////

void RemoteBlock3DConnector::prepareForSend()
{
   assert(sender  !=NULL); sender->prepareForSend();
}
//////////////////////////////////////////////////////////////////////////

void RemoteBlock3DConnector::sendVectors()
{
   assert(sender  !=NULL); sender->sendData();
}
//////////////////////////////////////////////////////////////////////////

void RemoteBlock3DConnector::prepareForReceive()
{
   assert(receiver!=NULL); receiver->prepareForReceive();
}
//////////////////////////////////////////////////////////////////////////

void RemoteBlock3DConnector::receiveVectors()
{
   assert(receiver!=NULL); receiver->receiveData();
}
