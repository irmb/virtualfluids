#include "CoarseToFineBlock3DConnector.h"



////////////////////////////////////////////////////////////////////////////

CoarseToFineBlock3DConnector::CoarseToFineBlock3DConnector(SPtr<Block3D> block,
   VectorTransmitterPtr sender00, VectorTransmitterPtr receiver00,
   VectorTransmitterPtr sender01, VectorTransmitterPtr receiver01,
   VectorTransmitterPtr sender10, VectorTransmitterPtr receiver10,
   VectorTransmitterPtr sender11, VectorTransmitterPtr receiver11,
   int sendDir, InterpolationProcessorPtr iprocessor) : Block3DConnector(sendDir)
   , block(block)
   , sender00(sender00)
   , sender01(sender01)
   , sender10(sender10)
   , sender11(sender11)
   , receiver00(receiver00)
   , receiver01(receiver01)
   , receiver10(receiver10)
   , receiver11(receiver11)
   , iprocessor(iprocessor)
{
   if (!(sendDir==D3Q27System::E  || sendDir==D3Q27System::W  || sendDir==D3Q27System::N  || sendDir==D3Q27System::S  || sendDir==D3Q27System::T || sendDir==D3Q27System::B
      ||  sendDir==D3Q27System::NE || sendDir==D3Q27System::SW || sendDir==D3Q27System::SE || sendDir==D3Q27System::NW
      ||  sendDir==D3Q27System::TE || sendDir==D3Q27System::BW || sendDir==D3Q27System::BE || sendDir==D3Q27System::TW
      ||  sendDir==D3Q27System::TN || sendDir==D3Q27System::BS || sendDir==D3Q27System::BN || sendDir==D3Q27System::TS
      ||  sendDir==D3Q27System::TNE || sendDir==D3Q27System::TNW || sendDir==D3Q27System::TSE || sendDir==D3Q27System::TSW
      ||  sendDir==D3Q27System::BNE || sendDir==D3Q27System::BNW || sendDir==D3Q27System::BSE || sendDir==D3Q27System::BSW
      ))
   {
      throw UbException(UB_EXARGS, "invalid constructor for this direction");
   }
}
//////////////////////////////////////////////////////////////////////////

bool CoarseToFineBlock3DConnector::isLocalConnector()
{
   return !this->isRemoteConnector();
}
//////////////////////////////////////////////////////////////////////////

bool CoarseToFineBlock3DConnector::isRemoteConnector()
{
   return ((sender11 && sender11->isRemoteTransmitter()) || (receiver11 && receiver11->isRemoteTransmitter())
      || (sender00 && sender00->isRemoteTransmitter()) || (receiver00 && receiver00->isRemoteTransmitter())
      || (sender01 && sender01->isRemoteTransmitter()) || (receiver01 && receiver01->isRemoteTransmitter())
      || (sender10 && sender10->isRemoteTransmitter()) || (receiver10 && receiver10->isRemoteTransmitter()));
}
//////////////////////////////////////////////////////////////////////////

void CoarseToFineBlock3DConnector::sendTransmitterDataSize()
{
   if (sender00)
   {
      UBLOG(logDEBUG5, "CoarseToFineBlock3DConnector<VectorTransmitter>::sendTransmitterDataSize()-sender00 "<<block.lock()->toString()<<" sendDir="<<sendDir);
      sender00->sendDataSize();
   }
   if (sender01)
   {
      UBLOG(logDEBUG5, "CoarseToFineBlock3DConnector<VectorTransmitter>::sendTransmitterDataSize()-sender01 "<<block.lock()->toString()<<"sendDir="<<sendDir);
      sender01->sendDataSize();
   }
   if (sender10)
   {
      UBLOG(logDEBUG5, "CoarseToFineBlock3DConnector<VectorTransmitter>::sendTransmitterDataSize()-sender10 "<<block.lock()->toString()+"sendDir="<<sendDir);
      sender10->sendDataSize();
   }
   if (sender11)
   {
      UBLOG(logDEBUG5, "CoarseToFineBlock3DConnector<VectorTransmitter>::sendTransmitterDataSize()-sender11 "<<block.lock()->toString()<<"sendDir="<<sendDir);
      sender11->sendDataSize();
   }
}
//////////////////////////////////////////////////////////////////////////

void CoarseToFineBlock3DConnector::receiveTransmitterDataSize()
{
   if (receiver00)
   {
      UBLOG(logDEBUG5, "CoarseToFineBlock3DConnector<VectorTransmitter>::receiveTransmitterDataSize()-receiver00 "<<block.lock()->toString()<<"sendDir="<<sendDir);
      receiver00->receiveDataSize();
   }
   if (receiver01)
   {
      UBLOG(logDEBUG5, "CoarseToFineBlock3DConnector<VectorTransmitter>::receiveTransmitterDataSize()-receiver01 "<<block.lock()->toString()<<"sendDir="<<sendDir);
      receiver01->receiveDataSize();
   }
   if (receiver10)
   {
      UBLOG(logDEBUG5, "CoarseToFineBlock3DConnector<VectorTransmitter>::receiveTransmitterDataSize()-receiver10 "<<block.lock()->toString()<<"sendDir="<<sendDir);
      receiver10->receiveDataSize();
   }
   if (receiver11)
   {
      UBLOG(logDEBUG5, "CoarseToFineBlock3DConnector<VectorTransmitter>::receiveTransmitterDataSize()-receiver11 "<<block.lock()->toString()<<"sendDir="<<sendDir);
      receiver11->receiveDataSize();
   }
}
//////////////////////////////////////////////////////////////////////////

void CoarseToFineBlock3DConnector::prepareForSend()
{
   if (sender00) sender00->prepareForSend();
   if (sender01) sender01->prepareForSend();
   if (sender10) sender10->prepareForSend();
   if (sender11) sender11->prepareForSend();
}
//////////////////////////////////////////////////////////////////////////

void CoarseToFineBlock3DConnector::sendVectors()
{
   if (sender00) sender00->sendData();
   if (sender01) sender01->sendData();
   if (sender10) sender10->sendData();
   if (sender11) sender11->sendData();
}
//////////////////////////////////////////////////////////////////////////

void CoarseToFineBlock3DConnector::prepareForReceive()
{
   if (receiver00) receiver00->prepareForReceive();
   if (receiver01) receiver01->prepareForReceive();
   if (receiver10) receiver10->prepareForReceive();
   if (receiver11) receiver11->prepareForReceive();
}
//////////////////////////////////////////////////////////////////////////

void CoarseToFineBlock3DConnector::receiveVectors()
{
   if (receiver00) receiver00->receiveData();
   if (receiver01) receiver01->receiveData();
   if (receiver10) receiver10->receiveData();
   if (receiver11) receiver11->receiveData();
}
