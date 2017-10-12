#include "FineToCoarseBlock3DConnector.h"

////////////////////////////////////////////////////////////////////////////
FineToCoarseBlock3DConnector::FineToCoarseBlock3DConnector(Block3DPtr block, VectorTransmitterPtr sender,
   VectorTransmitterPtr receiver, int sendDir,
   InterpolationProcessorPtr iprocessor,
   CFconnectorType connType)
   : Block3DConnector(sendDir)
   , block(block)
   , sender(sender)
   , receiver(receiver)
   , iprocessor(iprocessor)
   , connType(connType)
{
   if (!(sendDir==D3Q27System::E  || sendDir==D3Q27System::W  || sendDir==D3Q27System::N  || sendDir==D3Q27System::S  || sendDir==D3Q27System::T || sendDir==D3Q27System::B
      ||  sendDir==D3Q27System::NE || sendDir==D3Q27System::SW || sendDir==D3Q27System::SE || sendDir==D3Q27System::NW
      ||  sendDir==D3Q27System::TE || sendDir==D3Q27System::BW ||  sendDir==D3Q27System::BE || sendDir==D3Q27System::TW
      ||  sendDir==D3Q27System::TN || sendDir==D3Q27System::BS ||  sendDir==D3Q27System::BN || sendDir==D3Q27System::TS

      ||  sendDir==D3Q27System::TNE || sendDir==D3Q27System::TNW ||  sendDir==D3Q27System::TSE || sendDir==D3Q27System::TSW
      ||  sendDir==D3Q27System::BNE || sendDir==D3Q27System::BNW ||  sendDir==D3Q27System::BSE || sendDir==D3Q27System::BSW

      ))
   {
      throw UbException(UB_EXARGS, "invalid constructor for this direction");
   }
}
//////////////////////////////////////////////////////////////////////////
bool FineToCoarseBlock3DConnector::isLocalConnector()
{
   return !this->isRemoteConnector();
}
//////////////////////////////////////////////////////////////////////////
bool FineToCoarseBlock3DConnector::isRemoteConnector()
{
   return sender->isRemoteTransmitter()  ||  receiver->isRemoteTransmitter();
}
//////////////////////////////////////////////////////////////////////////
void FineToCoarseBlock3DConnector::sendTransmitterDataSize()
{
   if (sender)
   {
      UBLOG(logDEBUG5, "FineToCoarseBlock3DConnector<VectorTransmitter>::sendTransmitterDataSize()"<<block.lock()->toString()+"sendDir="<<sendDir);
      sender->sendDataSize();
   }
}
//////////////////////////////////////////////////////////////////////////
void FineToCoarseBlock3DConnector::receiveTransmitterDataSize()
{
   if (receiver)
   {
      UBLOG(logDEBUG5, "FineToCoarseBlock3DConnector<VectorTransmitter>::receiveTransmitterDataSize()"<<block.lock()->toString()<<"sendDir="<<sendDir);
      receiver->receiveDataSize();
   }
}
//////////////////////////////////////////////////////////////////////////
void FineToCoarseBlock3DConnector::prepareForSend()
{
   if (sender) sender->prepareForSend();
}
//////////////////////////////////////////////////////////////////////////
void FineToCoarseBlock3DConnector::sendVectors()
{
   if (sender) sender->sendData();
}
//////////////////////////////////////////////////////////////////////////
void FineToCoarseBlock3DConnector::prepareForReceive()
{
   if (receiver) receiver->prepareForReceive();
}
//////////////////////////////////////////////////////////////////////////
void FineToCoarseBlock3DConnector::receiveVectors()
{
   if (receiver) receiver->receiveData();
}


