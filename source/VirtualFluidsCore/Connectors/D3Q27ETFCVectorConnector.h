/**
* @file D3Q27ETFCVectorConnector.h
* @brief Interpolation from fine level to coarse.
* @author Kostyantyn Kucher
* @date 08.06.2011
*/
#ifndef D3Q27ETFCVECTORCONNECTOR_H
#define D3Q27ETFCVECTORCONNECTOR_H

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
enum CFconnectorType {EvenOddNW, EvenEvenSW, OddEvenSE, OddOddNE, none};

//daten werden in einen vector (dieser befindet sich im transmitter) kopiert
//der vector wird via transmitter uebertragen
//transmitter kann ein lokal, MPI, RCG, CTL oder was auch immer fuer ein
//transmitter sein, der von Transmitter abgeleitet ist ;-)

template< typename VectorTransmitter >
class D3Q27ETFCVectorConnector : public Block3DConnector
{
public:

protected:
   typedef typename VectorTransmitter::value_type  vector_type;
   typedef boost::shared_ptr< VectorTransmitter > VectorTransmitterPtr;
public:
   D3Q27ETFCVectorConnector(Block3DPtr block, VectorTransmitterPtr sender, VectorTransmitterPtr receiver, int sendDir, D3Q27InterpolationProcessorPtr iprocessor, CFconnectorType connType); 

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

   bool isInterpolationConnectorCF() { return false; }
   bool isInterpolationConnectorFC() { return true; }

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
   //gegenstelle muss "inversen" connector besitzen

   VectorTransmitterPtr sender, receiver;

   D3Q27InterpolationProcessorPtr iprocessor;

   CFconnectorType connType;

   void readICellF(DistributionArray3DPtr f, D3Q27ICell& icellF, const int& x1, const int& x2, const int& x3) ;
   void writeICellCtoData(vector_type& data, int& index, LBMReal* icellC);
   void writeNodeToVector(vector_type& data, int& index, LBMReal* inode);
   void getLocalMinMax(int& minX1, int& minX2, int& minX3, int& maxX1, int& maxX2, int& maxX3);
   void getLocalMinMaxCF(int gMax, int& lMin, int& lMax);
   void fillSendVector(DistributionArray3DPtr fFrom, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index);

   void distributeReceiveVector(DistributionArray3DPtr fTo, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index);
   void readICellFfromData(vector_type& data, int& index, D3Q27ICell& icellF);
   void writeICellF(DistributionArray3DPtr f, D3Q27ICell& icellF, const int& x1, const int& x2, const int& x3);
   void readNodeFromVector(vector_type& data, int& index, LBMReal* inode);
   void getLocalOffsets(const int& gMax, int& oMin);
   void getLocalMins(int& minX1, int& minX2, int& minX3, const int& oMinX1, const int& oMinX2, const int& oMinX3);
};
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
D3Q27ETFCVectorConnector<VectorTransmitter>::D3Q27ETFCVectorConnector(Block3DPtr block, VectorTransmitterPtr sender, 
                                                                      VectorTransmitterPtr receiver, int sendDir, 
                                                                      D3Q27InterpolationProcessorPtr iprocessor,
                                                                      CFconnectorType connType) 
                                                                      : Block3DConnector(sendDir)
                                                                      , block(block)
                                                                      , sender(sender)
                                                                      , receiver(receiver)
                                                                      , iprocessor(iprocessor)
                                                                      , connType(connType)
{
   if( !(   sendDir==D3Q27System::E || sendDir==D3Q27System::W || sendDir==D3Q27System::N 
      || sendDir==D3Q27System::S || sendDir==D3Q27System::T || sendDir==D3Q27System::B ) )
   {
      throw UbException(UB_EXARGS,"invalid constructor for this direction");
   }
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
bool D3Q27ETFCVectorConnector<VectorTransmitter>::isLocalConnector()
{ 
   return !this->isRemoteConnector(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
bool D3Q27ETFCVectorConnector<VectorTransmitter>::isRemoteConnector() 
{ 
   return sender->isRemoteTransmitter()  ||  receiver->isRemoteTransmitter();
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETFCVectorConnector<VectorTransmitter>::sendTransmitterDataSize()  
{ 
   if(sender) sender->sendDataSize(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETFCVectorConnector<VectorTransmitter>::receiveTransmitterDataSize()
{ 
   if(receiver) receiver->receiveDataSize(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETFCVectorConnector<VectorTransmitter>::prepareForSend()
{ 
   if(sender) sender->prepareForSend(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETFCVectorConnector<VectorTransmitter>::sendVectors()     
{ 
   if(sender) sender->sendData();
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETFCVectorConnector<VectorTransmitter>::prepareForReceive()     
{ 
   if(receiver) receiver->prepareForReceive(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETFCVectorConnector<VectorTransmitter>::receiveVectors() 
{ 
   if(receiver) receiver->receiveData(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETFCVectorConnector<VectorTransmitter>::init()
{
   using namespace D3Q27System;

   int maxX1 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX1();
   int maxX2 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX2();
   int maxX3 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX3();

   int       sendSize  = 0;
   LBMReal initValue = -999.0;

   int sendDataPerNode = 27/*f*/;
   int iCellSize = 1; //size of interpolation cell

   switch(this->sendDir)
   {		                  
   case E : case W : sendSize = maxX2*maxX3*sendDataPerNode*iCellSize; break; 
   case N : case S : sendSize = maxX1*maxX3*sendDataPerNode*iCellSize; break; 
   case T : case B : sendSize = maxX1*maxX2*sendDataPerNode*iCellSize; break; 
   default: throw UbException(UB_EXARGS,"direction not allowed in this constructor");
   }
    sender->getData().resize(sendSize, initValue);
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETFCVectorConnector< VectorTransmitter>::fillSendVectors() 
{ 
   using namespace D3Q27System;

   DistributionArray3DPtr  fFrom = block.lock()->getKernel()->getDataSet()->getFdistributions();
   int maxX1 = (int)fFrom->getNX1();
   int maxX2 = (int)fFrom->getNX2();
   int maxX3 = (int)fFrom->getNX3();
   int minX1 = 0;
   int minX2 = 0;
   int minX3 = 0;

   int oMinX1, oMinX2, oMinX3; 
   getLocalOffsets(maxX1, oMinX1);
   getLocalOffsets(maxX2, oMinX2);
   getLocalOffsets(maxX3, oMinX3);

   int lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3;
   int index = 0;
   vector_type& data = sender->getData();

   lMinX1 = minX1+1; lMinX2 = minX2+1; lMinX3 = minX3+1;
   lMaxX1 = maxX1-2; lMaxX2 = maxX2-2; lMaxX3 = maxX3-2;

   switch(sendDir)
   {
   case E: 
      getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
      getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
      lMinX1 = maxX1-7;
      lMaxX1 = lMinX1 + 1;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
   case W: 
      getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
      getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
      lMinX1 = 5;
      lMaxX1 = lMinX1 + 1;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;  
   case N:
      getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
      getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
      lMinX2 = maxX2-7;
      lMaxX2 = lMinX2 + 1;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
   case S:
      getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
      getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
      lMinX2 = 5;
      lMaxX2 = lMinX2 + 1;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
   case T:
      getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
      getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
      lMinX3 = maxX3-7;
      lMaxX3 = lMinX3 + 1;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
   case B:
      getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
      getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
      lMinX3 = 5;
      lMaxX3 = lMinX3 + 1;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCVectorConnector< VectorTransmitter>::fillSendVector(DistributionArray3DPtr fFrom, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index)
{
   int ix1, ix2, ix3;
   for (ix3=lMinX3; ix3<lMaxX3; ix3+=2)
   {
      for (ix2=lMinX2; ix2<lMaxX2; ix2+=2)
      {
         for (ix1=lMinX1; ix1<lMaxX1; ix1+=2)
         {
            LBMReal icellC[27];
            D3Q27ICell icellF;
            this->readICellF(fFrom, icellF, ix1, ix2, ix3);
            iprocessor->interpolateFineToCoarse(icellF, icellC);
            this->writeICellCtoData(data, index, icellC);
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCVectorConnector< VectorTransmitter>::readICellF(DistributionArray3DPtr f, D3Q27ICell& icellF, const int& x1, const int& x2, const int& x3) 
{
   f->getDistribution(icellF.BSW, x1, x2, x3);
   f->getDistribution(icellF.BSE, x1+1, x2, x3);
   f->getDistribution(icellF.BNW, x1, x2+1, x3);
   f->getDistribution(icellF.BNE, x1+1, x2+1, x3);
   f->getDistribution(icellF.TSW, x1, x2, x3+1);
   f->getDistribution(icellF.TSE, x1+1, x2, x3+1);
   f->getDistribution(icellF.TNW, x1, x2+1, x3+1);
   f->getDistribution(icellF.TNE, x1+1, x2+1, x3+1);
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCVectorConnector< VectorTransmitter>::writeICellCtoData(vector_type& data, int& index, LBMReal* icellC) 
{
   for (int i = D3Q27System::STARTF; i < D3Q27System::ENDF+1; i++)
   {
      data[index++] = icellC[i];
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCVectorConnector< VectorTransmitter>::getLocalMinMaxCF(int gMax, int& lMin, int& lMax)
{
   if (Utilities::isOdd(gMax))
   {
      if(connType == OddEvenSE || connType == OddOddNE)
      {
         lMin = 1;
         lMax = gMax;
      }
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCVectorConnector< VectorTransmitter>::distributeReceiveVectors() 
{
   using namespace D3Q27System;

   DistributionArray3DPtr  fTo = block.lock()->getKernel()->getDataSet()->getFdistributions();
   int maxX1 = (int)fTo->getNX1();
   int maxX2 = (int)fTo->getNX2();
   int maxX3 = (int)fTo->getNX3();
   int minX1 = 0;
   int minX2 = 0;
   int minX3 = 0;

   int lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3;
   int index = 0;
   vector_type& data = receiver->getData();

   lMinX1 = minX1; lMinX2 = minX2; lMinX3 = minX3;
   lMaxX1 = maxX1-1; lMaxX2 = maxX2-1; lMaxX3 = maxX3-1;

   switch(sendDir)
   {
   case E: 
      lMinX1 = maxX1-4;
      lMaxX1 = lMinX1 + 1;
      getLocalMinMaxCF(maxX2, lMinX2, lMaxX2);
      getLocalMinMaxCF(maxX3, lMinX3, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
   case W: 
      lMinX1 = 2;
      lMaxX1 = lMinX1 + 1;
      getLocalMinMaxCF(maxX2, lMinX2, lMaxX2);
      getLocalMinMaxCF(maxX3, lMinX3, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;  
   case N:
      lMinX2 = maxX2-4;
      lMaxX2 = lMinX2 + 1;
      getLocalMinMaxCF(maxX1, lMinX1, lMaxX1);
      getLocalMinMaxCF(maxX3, lMinX3, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
   case S:
      lMinX2 = 2;
      lMaxX2 = lMinX2 + 1;
      getLocalMinMaxCF(maxX1, lMinX1, lMaxX1);
      getLocalMinMaxCF(maxX3, lMinX3, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
   case T:
      lMinX3 = maxX3-4;
      lMaxX3 = lMinX3 + 1;
      getLocalMinMaxCF(maxX1, lMinX1, lMaxX1);
      getLocalMinMaxCF(maxX2, lMinX2, lMaxX2);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
   case B:
      lMinX3 = 2;
      lMaxX3 = lMinX3 + 1;
      getLocalMinMaxCF(maxX1, lMinX1, lMaxX1);
      getLocalMinMaxCF(maxX2, lMinX2, lMaxX2);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCVectorConnector< VectorTransmitter>::distributeReceiveVector(DistributionArray3DPtr fTo, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index)
{
   int ix1, ix2, ix3;
   for (ix3=lMinX3; ix3<lMaxX3; ix3+=2)
   {
      for (ix2=lMinX2; ix2<lMaxX2; ix2+=2)
      {
         for (ix1=lMinX1; ix1<lMaxX1; ix1+=2)
         {
            D3Q27ICell icellF;
            this->readICellFfromData(data, index, icellF);
            this->writeICellF(fTo, icellF, ix1, ix2, ix3);
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCVectorConnector< VectorTransmitter>::writeICellF(DistributionArray3DPtr f, D3Q27ICell& icellF, const int& x1, const int& x2, const int& x3) 
{
   f->setDistributionInv(icellF.BSW, x1, x2, x3);
   f->setDistributionInv(icellF.BSE, x1+1, x2, x3);
   f->setDistributionInv(icellF.BNW, x1, x2+1, x3);
   f->setDistributionInv(icellF.BNE, x1+1, x2+1, x3);
   f->setDistributionInv(icellF.TSW, x1, x2, x3+1);
   f->setDistributionInv(icellF.TSE, x1+1, x2, x3+1);
   f->setDistributionInv(icellF.TNW, x1, x2+1, x3+1);
   f->setDistributionInv(icellF.TNE, x1+1, x2+1, x3+1);
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCVectorConnector< VectorTransmitter>::readICellFfromData(vector_type& data, int& index, D3Q27ICell& icellF) 
{
   readNodeFromVector(data, index, icellF.BSW);
   readNodeFromVector(data, index, icellF.BSE);
   readNodeFromVector(data, index, icellF.BNW);
   readNodeFromVector(data, index, icellF.BNE);
   readNodeFromVector(data, index, icellF.TSW);
   readNodeFromVector(data, index, icellF.TSE);
   readNodeFromVector(data, index, icellF.TNW);
   readNodeFromVector(data, index, icellF.TNE);
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCVectorConnector< VectorTransmitter>::readNodeFromVector(vector_type& data, int& index, LBMReal* inode)
{
   for (int i = D3Q27System::STARTF; i < D3Q27System::ENDF+1; i++)
   {
      inode[i] = data[index++];
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCVectorConnector< VectorTransmitter>::getLocalMinMax(int& minX1, int& minX2, int& minX3, int& maxX1, int& maxX2, int& maxX3)
{
   using namespace D3Q27System;

   if(block.lock()->hasInterpolationFlagFC(E))
   {
      maxX1 -= 3;
   }
   if(block.lock()->hasInterpolationFlagFC(W))
   {
      minX1 += 4;
   }
   if(block.lock()->hasInterpolationFlagFC(N))
   {
      maxX2 -= 3;
   }
   if(block.lock()->hasInterpolationFlagFC(S))
   {
      minX2 += 4;
   }
   if(block.lock()->hasInterpolationFlagFC(T))
   {
      maxX3 -= 3;
   }
   if(block.lock()->hasInterpolationFlagFC(B))
   {
      minX3 += 4;
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCVectorConnector< VectorTransmitter>::getLocalOffsets(const int& gMax, int& oMin)
{
   if (Utilities::isEven(gMax))
   {
      oMin = 0;
   }
   if (Utilities::isOdd(gMax))
   {
      oMin = -1;
   }

}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCVectorConnector< VectorTransmitter>::getLocalMins(int& minX1, int& minX2, int& minX3, const int& oMinX1, const int& oMinX2, const int& oMinX3)
{
   using namespace D3Q27System;

   switch(sendDir)
   {
   case E: case W:
      if(connType == OddEvenSE)
         minX2 += oMinX2;
      if(connType == OddOddNE)
      {
         minX2 += oMinX2;
         minX3 += oMinX3;
      }
      if(connType == EvenOddNW)
         minX3 += oMinX3;
      break;
   case N: case S:
      if(connType == OddEvenSE)
         minX1 += oMinX1;
      if(connType == OddOddNE)
      {
         minX1 += oMinX1;
         minX3 += oMinX3;
      }
      if(connType == EvenOddNW)
         minX3 += oMinX3;
      break;
   case T: case B:
      if(connType == OddEvenSE)
         minX1 += oMinX1;
      if(connType == OddOddNE)
      {
         minX1 += oMinX1;
         minX2 += oMinX2;
      }
      if(connType == EvenOddNW)
         minX2 += oMinX2;
      break;
   }
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
double D3Q27ETFCVectorConnector<VectorTransmitter>::getSendRecieveTime()
{
   return 0;
}

#endif
