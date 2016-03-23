/**
* @file D3Q27ETCFVectorConnector.h
* @brief Interpolation from coarse level to fine.
* @author Kostyantyn Kucher
* @date 08.06.2011
*/
#ifndef D3Q27ETCFVECTORCONNECTOR_H
#define D3Q27ETCFVECTORCONNECTOR_H

#include <vector>

#include "basics/transmitter/TbTransmitter.h"
#include "Block3DConnector.h"
#include "D3Q27System.h"
#include "Block3D.h"
#include "LBMKernelETD3Q27.h"
#include "D3Q27InterpolationProcessor.h"
#include "MathUtil.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

class Block3D;

//daten werden in einen vector (dieser befindet sich im transmitter) kopiert
//der vector wird via transmitter uebertragen
//transmitter kann ein lokal, MPI, RCG, CTL oder was auch immer fuer ein
//transmitter sein, der von Transmitter abgeleitet ist ;-)

//sendrichtung:    E<->W     N<->S    T<->B
//  ---------       x3       x3        x2
// | NW | NE |      ^        ^         ^
// |----+----|      +-> x2   +->x1     +->x1
// | SW | SE |     
//  ---------
// NW==even-odd, SW==even-even, SE==odd-even, NE==odd-odd

template< typename VectorTransmitter >
class D3Q27ETCFVectorConnector : public Block3DConnector
{
protected:
   typedef typename VectorTransmitter::value_type  vector_type;
   typedef boost::shared_ptr< VectorTransmitter > VectorTransmitterPtr;
public:
   D3Q27ETCFVectorConnector(  Block3DPtr block,
      VectorTransmitterPtr senderEvenEvenSW, VectorTransmitterPtr receiverEvenEvenSW, 
      VectorTransmitterPtr senderEvenOddNW,  VectorTransmitterPtr receiverEvenOddNW, 
      VectorTransmitterPtr senderOddEvenSE,  VectorTransmitterPtr receiverOddEvenSE, 
      VectorTransmitterPtr senderOddOddNE,   VectorTransmitterPtr receiverOddOddNE,
      int sendDir, D3Q27InterpolationProcessorPtr iprocessor); 

   bool isLocalConnector();
   bool isRemoteConnector();
   void init();

   void sendTransmitterDataSize();
   void receiveTransmitterDataSize();

   void prepareForSend();
   void sendVectors();

   void prepareForReceive();
   void receiveVectors();

   void fillSendVectors();
   void distributeReceiveVectors();

   bool isInterpolationConnectorCF() { return true; }
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
   boost::weak_ptr<Block3D> block; //dieser nvd sendet daten und die empfangenen werden diesem nvd zugeordnet

   VectorTransmitterPtr senderEvenEvenSW, receiverEvenEvenSW, 
                        senderEvenOddNW,  receiverEvenOddNW, 
                        senderOddEvenSE,  receiverOddEvenSE, 
                        senderOddOddNE,   receiverOddOddNE;
   
   D3Q27InterpolationProcessorPtr iprocessor;

   void readICellC(DistributionArray3DPtr f, D3Q27ICell& icellC, const int& x1, const int& x2, const int& x3) ;
   void writeICellFtoData(vector_type& data, int& index, D3Q27ICell& icellF);
   void writeNodeToVector(vector_type& data, int& index, LBMReal* inode);
   void getLocalMinMax(const int& gMin, const int& gMax, const bool& even, int& lMin, int& lMax, const bool& dataDistribution);
   void getLocalMinMax(int& minX1, int& minX2, int& minX3, int& maxX1, int& maxX2, int& maxX3);
   void fillSendVectorExt(DistributionArray3DPtr fFrom, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index);

   void distributeReceiveVector(DistributionArray3DPtr fTo, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index);
   void writeICellC(DistributionArray3DPtr f, LBMReal* icellC, const int& x1, const int& x2, const int& x3);
   void readICellCfromData(vector_type& data, int& index, LBMReal* icellC);
};

//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
D3Q27ETCFVectorConnector<VectorTransmitter>::D3Q27ETCFVectorConnector(  Block3DPtr block,
                         VectorTransmitterPtr senderEvenEvenSW, VectorTransmitterPtr receiverEvenEvenSW, 
                         VectorTransmitterPtr senderEvenOddNW,  VectorTransmitterPtr receiverEvenOddNW, 
                         VectorTransmitterPtr senderOddEvenSE,  VectorTransmitterPtr receiverOddEvenSE, 
                         VectorTransmitterPtr senderOddOddNE,   VectorTransmitterPtr receiverOddOddNE,
                         int sendDir, D3Q27InterpolationProcessorPtr iprocessor) :  Block3DConnector(sendDir)
                         , block(block)
                         , senderEvenEvenSW(senderEvenEvenSW)
                         , senderEvenOddNW(senderEvenOddNW)
                         , senderOddEvenSE(senderOddEvenSE)
                         , senderOddOddNE(senderOddOddNE)
                         , receiverEvenEvenSW(receiverEvenEvenSW)
                         , receiverEvenOddNW(receiverEvenOddNW)
                         , receiverOddEvenSE(receiverOddEvenSE)
                         , receiverOddOddNE(receiverOddOddNE)
                         , iprocessor(iprocessor)
{
   if( !(   sendDir==D3Q27System::E || sendDir==D3Q27System::W || sendDir==D3Q27System::N 
      || sendDir==D3Q27System::S || sendDir==D3Q27System::T || sendDir==D3Q27System::B ) )
   {
      throw UbException(UB_EXARGS,"invalid constructor for this direction");
   }
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
bool D3Q27ETCFVectorConnector<VectorTransmitter>::isLocalConnector()
{ 
   return !this->isRemoteConnector(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
bool D3Q27ETCFVectorConnector<VectorTransmitter>::isRemoteConnector() 
{ 
   return (   ( senderOddOddNE && senderOddOddNE->isRemoteTransmitter() ) ||  ( receiverOddOddNE && receiverOddOddNE->isRemoteTransmitter() )
      || ( senderEvenEvenSW && senderEvenEvenSW->isRemoteTransmitter() ) ||  ( receiverEvenEvenSW && receiverEvenEvenSW->isRemoteTransmitter() )
      || ( senderEvenOddNW && senderEvenOddNW->isRemoteTransmitter() ) ||  ( receiverEvenOddNW && receiverEvenOddNW->isRemoteTransmitter() )
      || ( senderOddEvenSE && senderOddEvenSE->isRemoteTransmitter() ) ||  ( receiverOddEvenSE && receiverOddEvenSE->isRemoteTransmitter() ) );
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFVectorConnector<VectorTransmitter>::sendTransmitterDataSize()  
{ 
   if(senderEvenEvenSW) senderEvenEvenSW->sendDataSize(); 
   if(senderEvenOddNW) senderEvenOddNW->sendDataSize(); 
   if(senderOddEvenSE) senderOddEvenSE->sendDataSize(); 
   if(senderOddOddNE) senderOddOddNE->sendDataSize(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFVectorConnector<VectorTransmitter>::receiveTransmitterDataSize()
{ 
   if(receiverEvenEvenSW) receiverEvenEvenSW->receiveDataSize(); 
   if(receiverEvenOddNW) receiverEvenOddNW->receiveDataSize(); 
   if(receiverOddEvenSE) receiverOddEvenSE->receiveDataSize(); 
   if(receiverOddOddNE) receiverOddOddNE->receiveDataSize(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFVectorConnector<VectorTransmitter>::prepareForSend()
{ 
   if(senderEvenEvenSW) senderEvenEvenSW->prepareForSend(); 
   if(senderEvenOddNW) senderEvenOddNW->prepareForSend(); 
   if(senderOddEvenSE) senderOddEvenSE->prepareForSend(); 
   if(senderOddOddNE) senderOddOddNE->prepareForSend(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFVectorConnector<VectorTransmitter>::sendVectors()     
{ 
   if(senderEvenEvenSW) senderEvenEvenSW->sendData();
   if(senderEvenOddNW) senderEvenOddNW->sendData();
   if(senderOddEvenSE) senderOddEvenSE->sendData();
   if(senderOddOddNE) senderOddOddNE->sendData();
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFVectorConnector<VectorTransmitter>::prepareForReceive()     
{ 
   if(receiverEvenEvenSW) receiverEvenEvenSW->prepareForReceive(); 
   if(receiverEvenOddNW) receiverEvenOddNW->prepareForReceive(); 
   if(receiverOddEvenSE) receiverOddEvenSE->prepareForReceive(); 
   if(receiverOddOddNE) receiverOddOddNE->prepareForReceive(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFVectorConnector<VectorTransmitter>::receiveVectors() 
{ 
   if(receiverEvenEvenSW) receiverEvenEvenSW->receiveData(); 
   if(receiverEvenOddNW) receiverEvenOddNW->receiveData();  
   if(receiverOddEvenSE) receiverOddEvenSE->receiveData(); 
   if(receiverOddOddNE) receiverOddOddNE->receiveData(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFVectorConnector<VectorTransmitter>::init()
{
   using namespace D3Q27System;

   int maxX1 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX1();
   int maxX2 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX2();
   int maxX3 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX3();

   int       sendSize  = 0;
   LBMReal initValue = -999.0;

   int sendDataPerNode = 27/*f*/;
   int iCellSize = 8; //size of interpolation cell

   switch(this->sendDir)
   {		                  
   case E : case W : sendSize = maxX2*maxX3*sendDataPerNode*iCellSize; break; 
   case N : case S : sendSize = maxX1*maxX3*sendDataPerNode*iCellSize; break; 
   case T : case B : sendSize = maxX1*maxX2*sendDataPerNode*iCellSize; break; 
   default: throw UbException(UB_EXARGS,"direction not allowed in this constructor");
   }
   senderEvenEvenSW->getData().resize(sendSize, initValue);
   senderEvenOddNW->getData().resize(sendSize, initValue);
   senderOddEvenSE->getData().resize(sendSize, initValue);
   senderOddOddNE->getData().resize(sendSize, initValue);
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFVectorConnector< VectorTransmitter>::fillSendVectors() 
{ 
   using namespace D3Q27System;

   DistributionArray3DPtr  fFrom = block.lock()->getKernel()->getDataSet()->getFdistributions();
   int maxX1 = (int)fFrom->getNX1();
   int maxX2 = (int)fFrom->getNX2();
   int maxX3 = (int)fFrom->getNX3();
   int minX1 = 0;
   int minX2 = 0;
   int minX3 = 0;

   int indexEvEv = 0;
   int indexEvOd = 0;
   int indexOdEv = 0;
   int indexOdOd = 0;
   
   vector_type& dataEvEv = this->senderEvenEvenSW->getData();
   vector_type& dataEvOd = this->senderEvenOddNW->getData();
   vector_type& dataOdEv = this->senderOddEvenSE->getData();
   vector_type& dataOdOd = this->senderOddOddNE->getData();

   int lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3;
   
   switch(sendDir)
   {
   case E: 
      lMinX1 = maxX1-3;
      lMaxX1 = lMinX1 + 1;
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
      break;
   case W:
      lMinX1 = 1;
      lMaxX1 = lMinX1 + 1;
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
      break;
   case N:
      lMinX2 = maxX2-3;
      lMaxX2 = lMinX2 + 1;
      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
      break;
   case S:
      lMinX2 = 1;
      lMaxX2 = lMinX2 + 1;
      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
      break;
   case T:
      lMinX3 = maxX3-3;
      lMaxX3 = lMinX3 + 1;
      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
      break;
   case B:
      lMinX3 = 1;
      lMaxX3 = lMinX3 + 1;
      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
      break;
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETCFVectorConnector< VectorTransmitter>::getLocalMinMax(const int& gMin, const int& gMax, const bool& even, int& lMin, int& lMax, const bool& dataDistribution)
{
   int halfEven = 0;
   int halfOdd = 0;
   int dCoef = 0;

   if (dataDistribution)
      dCoef = 1;
   
   if (Utilities::isOdd(gMax))
   {
      halfEven = gMax/2;
      halfOdd =  gMax/2;
   }
   if (Utilities::isEven(gMax))
   {
      halfEven = gMax/2;
      halfOdd =  gMax/2 - 1 + dCoef;
   }

   switch (even)
   {
   case true :
      lMin = gMin + dCoef;
      lMax = lMin + halfEven - dCoef;
      break;
   case false :
      lMin = gMin + halfOdd;
      lMax = gMax-1;
      break;
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETCFVectorConnector< VectorTransmitter>::fillSendVectorExt(DistributionArray3DPtr fFrom, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index)
{
   int ix1, ix2, ix3;
   for (ix3=lMinX3; ix3<lMaxX3; ix3++)
   {
      for (ix2=lMinX2; ix2<lMaxX2; ix2++)
      {
         for (ix1=lMinX1; ix1<lMaxX1; ix1++)
         {
            D3Q27ICell icellC;
            D3Q27ICell icellF;
            this->readICellC(fFrom, icellC, ix1, ix2, ix3);
            iprocessor->interpolateCoarseToFine(icellC, icellF);
            this->writeICellFtoData(data, index, icellF);
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETCFVectorConnector< VectorTransmitter>::readICellC(DistributionArray3DPtr f, D3Q27ICell& icellC, const int& x1, const int& x2, const int& x3) 
{
   f->getDistribution(icellC.BSW, x1, x2, x3);
   f->getDistribution(icellC.BSE, x1+1, x2, x3);
   f->getDistribution(icellC.BNW, x1, x2+1, x3);
   f->getDistribution(icellC.BNE, x1+1, x2+1, x3);
   f->getDistribution(icellC.TSW, x1, x2, x3+1);
   f->getDistribution(icellC.TSE, x1+1, x2, x3+1);
   f->getDistribution(icellC.TNW, x1, x2+1, x3+1);
   f->getDistribution(icellC.TNE, x1+1, x2+1, x3+1);
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETCFVectorConnector< VectorTransmitter>::writeICellFtoData(vector_type& data, int& index, D3Q27ICell& icellF) 
{
   writeNodeToVector(data, index, icellF.BSW);
   writeNodeToVector(data, index, icellF.BSE);
   writeNodeToVector(data, index, icellF.BNW);
   writeNodeToVector(data, index, icellF.BNE);
   writeNodeToVector(data, index, icellF.TSW);
   writeNodeToVector(data, index, icellF.TSE);
   writeNodeToVector(data, index, icellF.TNW);
   writeNodeToVector(data, index, icellF.TNE);
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETCFVectorConnector< VectorTransmitter>::writeNodeToVector(vector_type& data, int& index, LBMReal* inode)
{
   for (int i = D3Q27System::STARTF; i < D3Q27System::ENDF+1; i++)
   {
      data[index++] = inode[i];
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETCFVectorConnector< VectorTransmitter>::distributeReceiveVectors() 
{
   using namespace D3Q27System;

   DistributionArray3DPtr  fTo = block.lock()->getKernel()->getDataSet()->getFdistributions();
   int maxX1 = (int)fTo->getNX1();
   int maxX2 = (int)fTo->getNX2();
   int maxX3 = (int)fTo->getNX3();
   int minX1 = 0;
   int minX2 = 0;
   int minX3 = 0;


   int indexEvEv = 0;
   int indexEvOd = 0;
   int indexOdEv = 0;
   int indexOdOd = 0;

   vector_type& dataEvEv = this->receiverEvenEvenSW->getData();
   vector_type& dataEvOd = this->receiverEvenOddNW->getData();
   vector_type& dataOdEv = this->receiverOddEvenSE->getData();
   vector_type& dataOdOd = this->receiverOddOddNE->getData();

   int lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3;
   int dummy;

   switch(sendDir)
   {
   case E: 
      lMinX1 = maxX1-4;
      lMaxX1 = lMinX1 + 1;
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, lMinX2, lMinX3, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, lMinX3, dummy, lMaxX2, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
      break;
   case W:
      lMinX1 = 3;
      lMaxX1 = lMinX1 + 1;
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, lMinX2, lMinX3, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, lMinX3, dummy, lMaxX2, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
      break;
   case N:
      lMinX2 = maxX2-4;
      lMaxX2 = lMinX2 + 1;
      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
      getLocalMinMax(lMinX1, dummy, lMinX3, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
      getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, lMinX3, lMaxX1, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
      break;
   case S:
      lMinX2 = 3;
      lMaxX2 = lMinX2 + 1;
      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
      getLocalMinMax(lMinX1, dummy, lMinX3, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
      getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, lMinX3, lMaxX1, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
      break;
   case T:
      lMinX3 = maxX3-4;
      lMaxX3 = lMinX3 + 1;
      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
      getLocalMinMax(lMinX1, lMinX2, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
      getLocalMinMax(lMinX1, dummy, dummy, dummy, lMaxX2, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, lMinX2, dummy, lMaxX1, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, dummy, dummy, lMaxX1, lMaxX2, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
      break;
   case B:
      lMinX3 = 3;
      lMaxX3 = lMinX3 + 1;
      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
      getLocalMinMax(lMinX1, lMinX2, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
      getLocalMinMax(lMinX1, dummy, dummy, dummy, lMaxX2, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, lMinX2, dummy, lMaxX1, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, dummy, dummy, lMaxX1, lMaxX2, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
      break;
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETCFVectorConnector< VectorTransmitter>::distributeReceiveVector(DistributionArray3DPtr fTo, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index)
{
   int ix1, ix2, ix3;
   for (ix3=lMinX3; ix3<lMaxX3; ix3++)
   {
      for (ix2=lMinX2; ix2<lMaxX2; ix2++)
      {
         for (ix1=lMinX1; ix1<lMaxX1; ix1++)
         {
            LBMReal icellC[27];
            this->readICellCfromData(data, index, icellC);
            this->writeICellC(fTo, icellC, ix1, ix2, ix3);
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETCFVectorConnector< VectorTransmitter>::readICellCfromData(vector_type& data, int& index, LBMReal* icellC) 
{
   for (int i = D3Q27System::STARTF; i < D3Q27System::ENDF+1; i++)
   {
      icellC[i] = data[index++];
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETCFVectorConnector< VectorTransmitter>::writeICellC(DistributionArray3DPtr f, LBMReal* icellC, const int& x1, const int& x2, const int& x3) 
{
   f->setDistributionInv(icellC, x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETCFVectorConnector< VectorTransmitter>::getLocalMinMax(int& minX1, int& minX2, int& minX3, int& maxX1, int& maxX2, int& maxX3)
{
   using namespace D3Q27System;

   if(block.lock()->hasInterpolationFlagCF(E))
   {
      maxX1 -= 2;
   }
   if(block.lock()->hasInterpolationFlagCF(W))
   {
      minX1 += 2;
   }
   if(block.lock()->hasInterpolationFlagCF(N))
   {
      maxX2 -= 2;  
   }
   if(block.lock()->hasInterpolationFlagCF(S))
   {
      minX2 += 2;
   }
   if(block.lock()->hasInterpolationFlagCF(T))
   {
      maxX3 -= 2;
   }
   if(block.lock()->hasInterpolationFlagCF(B))
   {
      minX3 += 2;
   }
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
double D3Q27ETCFVectorConnector<VectorTransmitter>::getSendRecieveTime()
{
   return 0;
}

#endif 
