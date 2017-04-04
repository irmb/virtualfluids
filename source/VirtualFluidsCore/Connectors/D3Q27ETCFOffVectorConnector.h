/**
* @file D3Q27ETCFOffVectorConnector.h
* @class D3Q27ETCFOffVectorConnector
* @brief Interpolation from coarse level to fine.
* @author Kostyantyn Kucher and Ehsan Fard
* @date 08.06.2011
*/
#ifndef D3Q27ETCFOffVectorConnector_H
#define D3Q27ETCFOffVectorConnector_H

#include <vector>

#include "basics/transmitter/TbTransmitter.h"
#include "basics/transmitter/TbTransmitterLocal.h"
#include "basics/container/CbVector.h"
#include "Block3DConnector.h"
#include "D3Q27System.h"
#include "Block3D.h"
#include "LBMKernel.h"
#include "InterpolationProcessor.h"
#include "MathUtil.hpp"
#include "Grid3D.h"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include "D3Q27ETFCOffVectorConnector.h"
#include "BCProcessor.h"

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
class D3Q27ETCFOffVectorConnector : public Block3DConnector
{
public:
   typedef typename VectorTransmitter::value_type  vector_type;
   typedef boost::shared_ptr< VectorTransmitter > VectorTransmitterPtr;
public:
   D3Q27ETCFOffVectorConnector(Block3DPtr block,
      VectorTransmitterPtr senderEvenEvenSW, VectorTransmitterPtr receiverEvenEvenSW,
      VectorTransmitterPtr senderEvenOddNW, VectorTransmitterPtr receiverEvenOddNW,
      VectorTransmitterPtr senderOddEvenSE, VectorTransmitterPtr receiverOddEvenSE,
      VectorTransmitterPtr senderOddOddNE, VectorTransmitterPtr receiverOddOddNE,
      int sendDir, InterpolationProcessorPtr iprocessor);

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
   boost::weak_ptr<Block3D> block; //dieser nvd sendet daten und die empfangenen werden diesem nvd zugeordnet
   VectorTransmitterPtr senderEvenEvenSW, receiverEvenEvenSW,
      senderEvenOddNW, receiverEvenOddNW,
      senderOddEvenSE, receiverOddEvenSE,
      senderOddOddNE, receiverOddOddNE;

   InterpolationProcessorPtr iprocessor;

   void writeICellFtoData(vector_type& data, int& index, D3Q27ICell& icellF);
   void writeNodeToVector(vector_type& data, int& index, LBMReal* inode);
   void getLocalMinMax(const int& gMin, const int& gMax, const bool& even, int& lMin, int& lMax, const bool& dataDistribution);
   void getLocalMinMax(int& minX1, int& minX2, int& minX3, int& maxX1, int& maxX2, int& maxX3);
   void getLocalMinMax(int& minX1, int& minX2, int& minX3, int& maxX1, int& maxX2, int& maxX3, CFconnectorType connType);
   void fillSendVectorExt(DistributionArray3DPtr fFrom, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index);

   void distributeReceiveVector(DistributionArray3DPtr fTo, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index);
   void readICellCfromData(vector_type& data, int& index, LBMReal* icellC);

   void findCFnodes();
   void findCFnodes(DistributionArray3DPtr fFrom, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index);

   int bMaxX1, bMaxX2, bMaxX3;
};

////////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
D3Q27ETCFOffVectorConnector<VectorTransmitter>::D3Q27ETCFOffVectorConnector(Block3DPtr block,
   VectorTransmitterPtr senderEvenEvenSW, VectorTransmitterPtr receiverEvenEvenSW,
   VectorTransmitterPtr senderEvenOddNW, VectorTransmitterPtr receiverEvenOddNW,
   VectorTransmitterPtr senderOddEvenSE, VectorTransmitterPtr receiverOddEvenSE,
   VectorTransmitterPtr senderOddOddNE, VectorTransmitterPtr receiverOddOddNE,
   int sendDir, InterpolationProcessorPtr iprocessor) : Block3DConnector(sendDir)
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
   if (!(sendDir == D3Q27System::E || sendDir == D3Q27System::W || sendDir == D3Q27System::N || sendDir == D3Q27System::S || sendDir == D3Q27System::T || sendDir == D3Q27System::B
      || sendDir == D3Q27System::NE || sendDir == D3Q27System::SW || sendDir == D3Q27System::SE || sendDir == D3Q27System::NW
      || sendDir == D3Q27System::TE || sendDir == D3Q27System::BW || sendDir == D3Q27System::BE || sendDir == D3Q27System::TW
      || sendDir == D3Q27System::TN || sendDir == D3Q27System::BS || sendDir == D3Q27System::BN || sendDir == D3Q27System::TS
      || sendDir == D3Q27System::TNE || sendDir == D3Q27System::TNW || sendDir == D3Q27System::TSE || sendDir == D3Q27System::TSW
      || sendDir == D3Q27System::BNE || sendDir == D3Q27System::BNW || sendDir == D3Q27System::BSE || sendDir == D3Q27System::BSW
      ))
   {
      throw UbException(UB_EXARGS, "invalid constructor for this direction");
   }
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
bool D3Q27ETCFOffVectorConnector<VectorTransmitter>::isLocalConnector()
{
   return !this->isRemoteConnector();
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
bool D3Q27ETCFOffVectorConnector<VectorTransmitter>::isRemoteConnector()
{
   return ((senderOddOddNE && senderOddOddNE->isRemoteTransmitter()) || (receiverOddOddNE && receiverOddOddNE->isRemoteTransmitter())
      || (senderEvenEvenSW && senderEvenEvenSW->isRemoteTransmitter()) || (receiverEvenEvenSW && receiverEvenEvenSW->isRemoteTransmitter())
      || (senderEvenOddNW && senderEvenOddNW->isRemoteTransmitter()) || (receiverEvenOddNW && receiverEvenOddNW->isRemoteTransmitter())
      || (senderOddEvenSE && senderOddEvenSE->isRemoteTransmitter()) || (receiverOddEvenSE && receiverOddEvenSE->isRemoteTransmitter()));
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFOffVectorConnector<VectorTransmitter>::sendTransmitterDataSize()
{
   if (senderEvenEvenSW)
   {
      UBLOG(logDEBUG5, "D3Q27ETCFOffVectorConnector<VectorTransmitter>::sendTransmitterDataSize()-senderEvenEvenSW " << block.lock()->toString() << " sendDir=" << sendDir);
      senderEvenEvenSW->sendDataSize();
   }
   if (senderEvenOddNW)
   {
      UBLOG(logDEBUG5, "D3Q27ETCFOffVectorConnector<VectorTransmitter>::sendTransmitterDataSize()-senderEvenOddNW " << block.lock()->toString() << "sendDir=" << sendDir);
      senderEvenOddNW->sendDataSize();
   }
   if (senderOddEvenSE)
   {
      UBLOG(logDEBUG5, "D3Q27ETCFOffVectorConnector<VectorTransmitter>::sendTransmitterDataSize()-senderOddEvenSE " << block.lock()->toString() + "sendDir=" << sendDir);
      senderOddEvenSE->sendDataSize();
   }
   if (senderOddOddNE)
   {
      UBLOG(logDEBUG5, "D3Q27ETCFOffVectorConnector<VectorTransmitter>::sendTransmitterDataSize()-senderOddOddNE " << block.lock()->toString() << "sendDir=" << sendDir);
      senderOddOddNE->sendDataSize();
   }
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFOffVectorConnector<VectorTransmitter>::receiveTransmitterDataSize()
{
   if (receiverEvenEvenSW)
   {
      UBLOG(logDEBUG5, "D3Q27ETCFOffVectorConnector<VectorTransmitter>::receiveTransmitterDataSize()-receiverEvenEvenSW " << block.lock()->toString() << "sendDir=" << sendDir);
      receiverEvenEvenSW->receiveDataSize();
   }
   if (receiverEvenOddNW)
   {
      UBLOG(logDEBUG5, "D3Q27ETCFOffVectorConnector<VectorTransmitter>::receiveTransmitterDataSize()-receiverEvenOddNW " << block.lock()->toString() << "sendDir=" << sendDir);
      receiverEvenOddNW->receiveDataSize();
   }
   if (receiverOddEvenSE)
   {
      UBLOG(logDEBUG5, "D3Q27ETCFOffVectorConnector<VectorTransmitter>::receiveTransmitterDataSize()-receiverOddEvenSE " << block.lock()->toString() << "sendDir=" << sendDir);
      receiverOddEvenSE->receiveDataSize();
   }
   if (receiverOddOddNE)
   {
      UBLOG(logDEBUG5, "D3Q27ETCFOffVectorConnector<VectorTransmitter>::receiveTransmitterDataSize()-receiverOddOddNE " << block.lock()->toString() << "sendDir=" << sendDir);
      receiverOddOddNE->receiveDataSize();
   }
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFOffVectorConnector<VectorTransmitter>::prepareForSend()
{
   if (senderEvenEvenSW) senderEvenEvenSW->prepareForSend();
   if (senderEvenOddNW) senderEvenOddNW->prepareForSend();
   if (senderOddEvenSE) senderOddEvenSE->prepareForSend();
   if (senderOddOddNE) senderOddOddNE->prepareForSend();
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFOffVectorConnector<VectorTransmitter>::sendVectors()
{
   if (senderEvenEvenSW) senderEvenEvenSW->sendData();
   if (senderEvenOddNW) senderEvenOddNW->sendData();
   if (senderOddEvenSE) senderOddEvenSE->sendData();
   if (senderOddOddNE) senderOddOddNE->sendData();
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFOffVectorConnector<VectorTransmitter>::prepareForReceive()
{
   if (receiverEvenEvenSW) receiverEvenEvenSW->prepareForReceive();
   if (receiverEvenOddNW) receiverEvenOddNW->prepareForReceive();
   if (receiverOddEvenSE) receiverOddEvenSE->prepareForReceive();
   if (receiverOddOddNE) receiverOddOddNE->prepareForReceive();
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFOffVectorConnector<VectorTransmitter>::receiveVectors()
{
   if (receiverEvenEvenSW) receiverEvenEvenSW->receiveData();
   if (receiverEvenOddNW) receiverEvenOddNW->receiveData();
   if (receiverOddEvenSE) receiverOddEvenSE->receiveData();
   if (receiverOddOddNE) receiverOddOddNE->receiveData();
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFOffVectorConnector<VectorTransmitter>::init()
{
   using namespace D3Q27System;

   bMaxX1 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX1();
   bMaxX2 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX2();
   bMaxX3 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX3();

   int       sendSize = 0;
   LBMReal initValue = -999.0;

   int sendDataPerNode = 27/*f*/;
   int iCellSize = 8; //size of interpolation cell

   switch (this->sendDir)
   {
   case E: case W: sendSize = bMaxX2*bMaxX3*sendDataPerNode*iCellSize; break;
   case N: case S: sendSize = bMaxX1*bMaxX3*sendDataPerNode*iCellSize; break;
   case T: case B: sendSize = bMaxX1*bMaxX2*sendDataPerNode*iCellSize; break;
   case NE: case SW:case SE: case NW: sendSize = 2 * bMaxX3*sendDataPerNode*iCellSize; break;
   case TE: case BW:case BE: case TW: sendSize = 2 * bMaxX2*sendDataPerNode*iCellSize; break;
   case TN: case BS:case BN: case TS: sendSize = 2 * bMaxX1*sendDataPerNode*iCellSize; break;
   case TNE: case TNW:case TSE: case TSW:case BNE: case BNW:case BSE: case BSW: sendSize = 6 * bMaxX1*sendDataPerNode*iCellSize; break;
   default: throw UbException(UB_EXARGS, "direction not allowed in this constructor");
   }
   if (senderEvenEvenSW) senderEvenEvenSW->getData().resize(sendSize, initValue);
   else senderEvenEvenSW = VectorTransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());
   if (senderEvenOddNW)  senderEvenOddNW->getData().resize(sendSize, initValue);
   else senderEvenOddNW = VectorTransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());
   if (senderOddEvenSE)  senderOddEvenSE->getData().resize(sendSize, initValue);
   else senderOddEvenSE = VectorTransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());
   if (senderOddOddNE)   senderOddOddNE->getData().resize(sendSize, initValue);
   else senderOddOddNE = VectorTransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());

   if (!receiverEvenEvenSW) receiverEvenEvenSW = VectorTransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());
   if (!receiverEvenOddNW)  receiverEvenOddNW = VectorTransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());
   if (!receiverOddEvenSE)  receiverOddEvenSE = VectorTransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());
   if (!receiverOddOddNE)   receiverOddOddNE = VectorTransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());

   //findCFnodes();
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFOffVectorConnector< VectorTransmitter>::fillSendVectors()
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
   //int lMinX1_2, lMinX2_2, lMinX3_2, lMaxX1_2, lMaxX2_2, lMaxX3_2;

   //for coners
   int lMinX1W = 1;
   int lMaxX1W = 2;

   int lMinX1E = maxX1 - 3;
   int lMaxX1E = maxX1 - 2;

   int lMinX2S = 1;
   int lMaxX2S = 2;

   int lMinX2N = maxX2 - 3;
   int lMaxX2N = maxX2 - 2;

   int lMinX3B = 1;
   int lMaxX3B = 2;

   int lMinX3T = maxX3 - 3;
   int lMaxX3T = maxX3 - 2;


   switch (sendDir)
   {
   case E:
      lMinX1 = maxX1 - 3;
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
      ///////////////////////////////////////
      ///DEBUG
      //if (block.lock()->getGlobalID() == 5780)
      //{
      //   int test = 0;
      //}
      //////////////
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
      lMinX2 = maxX2 - 3;
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
      lMinX3 = maxX3 - 3;
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
      ///N-S-E-W
   case NE:
      lMinX1 = maxX1 - 3;
      lMaxX1 = lMinX1 + 2;
      lMinX2 = maxX2 - 3;
      lMaxX2 = lMinX2 + 2;

      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case SW:
      lMinX1 = 0;
      lMaxX1 = lMinX1 + 2;
      lMinX2 = 0;
      lMaxX2 = lMinX2 + 2;

      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case SE:
      lMinX1 = maxX1 - 3;
      lMaxX1 = lMinX1 + 2;
      lMinX2 = 0;
      lMaxX2 = lMinX2 + 2;

      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case NW:
      lMinX1 = 0;
      lMaxX1 = lMinX1 + 2;
      lMinX2 = maxX2 - 3;
      lMaxX2 = lMinX2 + 2;

      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;
      /////T-B-E-W
   case TE:
      lMinX1 = maxX1 - 3;
      lMaxX1 = lMinX1 + 2;
      lMinX3 = maxX3 - 3;
      lMaxX3 = lMinX3 + 2;

      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case BW:
      lMinX1 = 0;
      lMaxX1 = lMinX1 + 2;
      lMinX3 = 0;
      lMaxX3 = lMinX3 + 2;

      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case BE:
      lMinX1 = maxX1 - 3;
      lMaxX1 = lMinX1 + 2;
      lMinX3 = 0;
      lMaxX3 = lMinX3 + 2;

      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case TW:
      lMinX1 = 0;
      lMaxX1 = lMinX1 + 2;
      lMinX3 = maxX3 - 3;
      lMaxX3 = lMinX3 + 2;

      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;
      ////
      /////T-B-N-S
   case TN:
      lMinX2 = maxX2 - 3;
      lMaxX2 = lMinX2 + 2;
      lMinX3 = maxX3 - 3;
      lMaxX3 = lMinX3 + 2;

      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case BS:
      lMinX2 = 0;
      lMaxX2 = lMinX2 + 2;
      lMinX3 = 0;
      lMaxX3 = lMinX3 + 2;

      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case BN:
      lMinX2 = maxX2 - 3;
      lMaxX2 = lMinX2 + 2;
      lMinX3 = 0;
      lMaxX3 = lMinX3 + 2;

      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case TS:
      lMinX2 = 0;
      lMaxX2 = lMinX2 + 2;
      lMinX3 = maxX3 - 3;
      lMaxX3 = lMinX3 + 2;

      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, false);
      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;


      //TNE
   case TNE:
      lMinX1 = maxX1 - 3;
      lMaxX1 = maxX1 - 1;
      lMinX2 = maxX2 - 3;
      lMaxX2 = maxX2 - 1;
      lMinX3 = maxX3 - 3;
      lMaxX3 = maxX3 - 1;

      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
      break;
      //   TNW
   case TNW:
      lMinX1 = 0;
      lMaxX1 = 2;
      lMinX2 = maxX2 - 3;
      lMaxX2 = maxX2 - 1;
      lMinX3 = maxX3 - 3;
      lMaxX3 = maxX3 - 1;

      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
      break;
      //   TSE
   case TSE:
      lMinX1 = maxX1 - 3;
      lMaxX1 = maxX1 - 1;
      lMinX2 = 0;
      lMaxX2 = 2;
      lMinX3 = maxX3 - 3;
      lMaxX3 = maxX3 - 1;

      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
      break;
      //   TSW
   case TSW:
      lMinX1 = 0;
      lMaxX1 = 2;
      lMinX2 = 0;
      lMaxX2 = 2;
      lMinX3 = maxX3 - 3;
      lMaxX3 = maxX3 - 1;

      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
      break;
      //   BNE
   case BNE:
      lMinX1 = maxX1 - 3;
      lMaxX1 = maxX1 - 1;
      lMinX2 = maxX2 - 3;
      lMaxX2 = maxX2 - 1;
      lMinX3 = 0;
      lMaxX3 = 2;

      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
      break;
      //   BNW
   case BNW:
      lMinX1 = 0;
      lMaxX1 = 2;
      lMinX2 = maxX2 - 3;
      lMaxX2 = maxX2 - 1;
      lMinX3 = 0;
      lMaxX3 = 2;

      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
      break;
      //   BSE
   case BSE:
      lMinX1 = maxX1 - 3;
      lMaxX1 = maxX1 - 1;
      lMinX2 = 0;
      lMaxX2 = 2;
      lMinX3 = 0;
      lMaxX3 = 2;

      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
      break;
      //   BSW
   case BSW:
      lMinX1 = 0;
      lMaxX1 = 2;
      lMinX2 = 0;
      lMaxX2 = 2;
      lMinX3 = 0;
      lMaxX3 = 2;

      fillSendVectorExt(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
      break;
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  >
void D3Q27ETCFOffVectorConnector< VectorTransmitter>::getLocalMinMax(const int& gMin, const int& gMax, const bool& even, int& lMin, int& lMax, const bool& dataDistribution)
{
   int halfEven = 0;
   int halfOdd = 0;
   int dCoef = 0;

   if (dataDistribution)
      dCoef = 1;

   if (Utilities::isOdd(gMax))
   {
      halfEven = gMax / 2;
      halfOdd = gMax / 2;
   }
   if (Utilities::isEven(gMax))
   {
      halfEven = gMax / 2;
      halfOdd = gMax / 2 - 1 + dCoef;
   }

   switch (even)
   {
   case true:
      lMin = gMin + dCoef;
      lMax = lMin + halfEven - dCoef;
      break;
   case false:
      lMin = gMin + halfOdd;
      lMax = gMax - 1;
      break;
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  >
void D3Q27ETCFOffVectorConnector< VectorTransmitter>::fillSendVectorExt(DistributionArray3DPtr fFrom, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index)
{
   if (data.size() == 0) return;
   int ix1, ix2, ix3;
   LBMReal xoff, yoff, zoff;
   BCArray3D& bcArray = block.lock()->getKernel()->getBCProcessor()->getBCArray();

   for (ix3 = lMinX3; ix3 < lMaxX3; ix3++)
   {
      for (ix2 = lMinX2; ix2 < lMaxX2; ix2++)
      {
         for (ix1 = lMinX1; ix1 < lMaxX1; ix1++)
         {
            D3Q27ICell icellC;
            D3Q27ICell icellF;

            int howManySolids = iprocessor->iCellHowManySolids(bcArray, ix1, ix2, ix3);

            if (howManySolids == 0 || howManySolids == 8)
            {
               iprocessor->readICell(fFrom, icellC, ix1, ix2, ix3);
               xoff = 0.0;
               yoff = 0.0;
               zoff = 0.0;
            }
            else
            {
               if (!iprocessor->findNeighborICell(bcArray, fFrom, icellC, bMaxX1, bMaxX2, bMaxX3, ix1, ix2, ix3, xoff, yoff, zoff))
               {
                  std::string err = "For " + block.lock()->toString() + " x1=" + UbSystem::toString(ix1) + ", x2=" + UbSystem::toString(ix2) + ", x3=" + UbSystem::toString(ix3) +
                     " interpolation is not implemented for other direction" +
                     " by using in: " + (std::string)typeid(*this).name() +
                     " or maybe you have a solid on the block boundary";
                  UB_THROW(UbException(UB_EXARGS, err));
               }
            }

            iprocessor->interpolateCoarseToFine(icellC, icellF, xoff, yoff, zoff);
            this->writeICellFtoData(data, index, icellF);
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  >
void D3Q27ETCFOffVectorConnector< VectorTransmitter>::writeICellFtoData(vector_type& data, int& index, D3Q27ICell& icellF)
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
void D3Q27ETCFOffVectorConnector< VectorTransmitter>::writeNodeToVector(vector_type& data, int& index, LBMReal* inode)
{
   for (int i = D3Q27System::STARTF; i < D3Q27System::ENDF + 1; i++)
   {
      data[index++] = inode[i];
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  >
void D3Q27ETCFOffVectorConnector< VectorTransmitter>::distributeReceiveVectors()
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

   //for coners
   int lMinX1W = 3;
   int lMaxX1W = 3;

   int lMinX1E = maxX1 - 3;
   int lMaxX1E = maxX1 - 2;

   int lMinX2S = 1;
   int lMaxX2S = 3;

   int lMinX2N = maxX2 - 3;
   int lMaxX2N = maxX2 - 2;

   int lMinX3B = 1;
   int lMaxX3B = 3;

   int lMinX3T = maxX3 - 3;
   int lMaxX3T = maxX3 - 2;

   ///////////////////////////////////////
   ///DEBUG
   //if (block.lock()->getGlobalID() == 5780)
   //{
   //   int test = 0;
   //}
   //////////////

   switch (sendDir)
   {
   case E:
      lMinX1 = maxX1 - 4;
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
      ///////////////////////////////////////
      ///DEBUG
      //if (block.lock()->getGlobalID() == 5780)
      //{
      //   int test = 0;
      //}
      //////////////
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
      lMinX2 = maxX2 - 4;
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
      lMinX3 = maxX3 - 4;
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

      //	/////E-W-N-S
   case NE:
      lMinX1 = maxX1 - 4;
      lMaxX1 = lMinX1 + 3;
      lMinX2 = maxX2 - 4;
      lMaxX2 = lMinX2 + 1;
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, lMinX3, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, dummy, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      lMinX1 = maxX1 - 4;
      lMaxX1 = lMinX1 + 1;
      lMinX2 = maxX2 - 4;
      lMaxX2 = lMinX2 + 3;
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, lMinX3, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, dummy, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case SW:
      lMinX1 = 1;
      lMaxX1 = lMinX1 + 3;
      lMinX2 = 3;
      lMaxX2 = lMinX2 + 1;
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, lMinX3, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, dummy, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      lMinX1 = 3;
      lMaxX1 = lMinX1 + 1;
      lMinX2 = 1;
      lMaxX2 = lMinX2 + 3;
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, lMinX3, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, dummy, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case SE:
      lMinX1 = maxX1 - 4;
      lMaxX1 = lMinX1 + 3;
      lMinX2 = 3;
      lMaxX2 = lMinX2 + 1;
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, lMinX3, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, dummy, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      lMinX1 = maxX1 - 4;
      lMaxX1 = lMinX1 + 1;
      lMinX2 = 1;
      lMaxX2 = lMinX2 + 3;
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, lMinX3, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, dummy, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case NW:
      lMinX1 = 1;
      lMaxX1 = lMinX1 + 3;
      lMinX2 = maxX2 - 4;
      lMaxX2 = lMinX2 + 1;
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, lMinX3, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, dummy, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      lMinX1 = 3;
      lMaxX1 = lMinX1 + 1;
      lMinX2 = maxX2 - 4;
      lMaxX2 = lMinX2 + 3;
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, lMinX3, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, dummy, lMaxX3);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;
      //		/////T-B-E-W
   case TE:
      lMinX1 = maxX1 - 4;
      lMaxX1 = lMinX1 + 3;
      lMinX3 = maxX3 - 4;
      lMaxX3 = lMinX3 + 1;
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      lMinX1 = maxX1 - 4;
      lMaxX1 = lMinX1 + 1;
      lMinX3 = maxX3 - 4;
      lMaxX3 = lMinX3 + 3;
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case BW:
      lMinX1 = 1;
      lMaxX1 = lMinX1 + 3;
      lMinX3 = 3;
      lMaxX3 = lMinX3 + 1;
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      lMinX1 = 3;
      lMaxX1 = lMinX1 + 1;
      lMinX3 = 1;
      lMaxX3 = lMinX3 + 3;
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case BE:
      lMinX1 = maxX1 - 4;
      lMaxX1 = lMinX1 + 3;
      lMinX3 = 3;
      lMaxX3 = lMinX3 + 1;
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      lMinX1 = maxX1 - 4;
      lMaxX1 = lMinX1 + 1;
      lMinX3 = 1;
      lMaxX3 = lMinX3 + 3;
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case TW:
      lMinX1 = 1;
      lMaxX1 = lMinX1 + 3;
      lMinX3 = maxX3 - 4;
      lMaxX3 = lMinX3 + 1;
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      lMinX1 = 3;
      lMaxX1 = lMinX1 + 1;
      lMinX3 = maxX3 - 4;
      lMaxX3 = lMinX3 + 3;
      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, lMinX2, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, true);
      getLocalMinMax(dummy, dummy, dummy, dummy, lMaxX2, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

      /////////////////////////T-N-B-S
   case TN:
      lMinX2 = maxX2 - 4;
      lMaxX2 = lMinX2 + 3;
      lMinX3 = maxX3 - 4;
      lMaxX3 = lMinX3 + 1;
      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
      getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
      getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      lMinX2 = maxX2 - 4;
      lMaxX2 = lMinX2 + 1;
      lMinX3 = maxX3 - 4;
      lMaxX3 = lMinX3 + 3;
      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
      getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
      getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case BS:
      lMinX2 = 1;
      lMaxX2 = lMinX2 + 3;
      lMinX3 = 3;
      lMaxX3 = lMinX3 + 1;
      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
      getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
      getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      lMinX2 = 3;
      lMaxX2 = lMinX2 + 1;
      lMinX3 = 1;
      lMaxX3 = lMinX3 + 3;
      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
      getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
      getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;


   case BN:
      lMinX2 = maxX2 - 4;
      lMaxX2 = lMinX2 + 3;
      lMinX3 = 3;
      lMaxX3 = lMinX3 + 1;
      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
      getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
      getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      lMinX2 = maxX2 - 4;
      lMaxX2 = lMinX2 + 1;
      lMinX3 = 1;
      lMaxX3 = lMinX3 + 3;
      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
      getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
      getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

   case TS:
      lMinX2 = 1;
      lMaxX2 = lMinX2 + 3;
      lMinX3 = maxX3 - 4;
      lMaxX3 = lMinX3 + 1;
      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
      getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
      getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      lMinX2 = 3;
      lMaxX2 = lMinX2 + 1;
      lMinX3 = maxX3 - 4;
      lMaxX3 = lMinX3 + 3;
      getLocalMinMax(minX1, maxX1, true, lMinX1, lMaxX1, true);
      getLocalMinMax(lMinX1, dummy, dummy, dummy, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1, maxX1, false, lMinX1, lMaxX1, true);
      getLocalMinMax(dummy, dummy, dummy, lMaxX1, dummy, dummy);
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      break;

      //TNE
   case TNE:
      lMinX1 = maxX1 - 4;
      lMaxX1 = maxX1 - 3;
      lMinX2 = maxX2 - 4;
      lMaxX2 = maxX2 - 1;
      lMinX3 = maxX3 - 4;
      lMaxX3 = maxX3 - 1;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      lMinX1 = maxX1 - 4;
      lMaxX1 = maxX1 - 1;
      lMinX2 = maxX2 - 4;
      lMaxX2 = maxX2 - 3;
      lMinX3 = maxX3 - 4;
      lMaxX3 = maxX3 - 1;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      lMinX1 = maxX1 - 4;
      lMaxX1 = maxX1 - 1;
      lMinX2 = maxX2 - 4;
      lMaxX2 = maxX2 - 1;
      lMinX3 = maxX3 - 4;
      lMaxX3 = maxX3 - 3;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      break;
      //   TNW
   case TNW:
      lMinX1 = 3;
      lMaxX1 = 4;
      lMinX2 = maxX2 - 4;
      lMaxX2 = maxX2 - 1;
      lMinX3 = maxX3 - 4;
      lMaxX3 = maxX3 - 1;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      lMinX1 = 1;
      lMaxX1 = 4;
      lMinX2 = maxX2 - 4;
      lMaxX2 = maxX2 - 3;
      lMinX3 = maxX3 - 4;
      lMaxX3 = maxX3 - 1;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      lMinX1 = 1;
      lMaxX1 = 4;
      lMinX2 = maxX2 - 4;
      lMaxX2 = maxX2 - 1;
      lMinX3 = maxX3 - 4;
      lMaxX3 = maxX3 - 3;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      break;
      //   TSE
   case TSE:
      lMinX1 = maxX1 - 4;
      lMaxX1 = maxX1 - 3;
      lMinX2 = 1;
      lMaxX2 = 4;
      lMinX3 = maxX3 - 4;
      lMaxX3 = maxX3 - 1;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      lMinX1 = maxX1 - 4;
      lMaxX1 = maxX1 - 1;
      lMinX2 = 3;
      lMaxX2 = 4;
      lMinX3 = maxX3 - 4;
      lMaxX3 = maxX3 - 1;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      lMinX1 = maxX1 - 4;
      lMaxX1 = maxX1 - 1;
      lMinX2 = 1;
      lMaxX2 = 4;
      lMinX3 = maxX3 - 4;
      lMaxX3 = maxX3 - 3;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
      break;
      //   TSW
   case TSW:
      lMinX1 = 3;
      lMaxX1 = 4;
      lMinX2 = 1;
      lMaxX2 = 4;
      lMinX3 = maxX3 - 4;
      lMaxX3 = maxX3 - 1;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      lMinX1 = 1;
      lMaxX1 = 4;
      lMinX2 = 3;
      lMaxX2 = 4;
      lMinX3 = maxX3 - 4;
      lMaxX3 = maxX3 - 1;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      lMinX1 = 1;
      lMaxX1 = 4;
      lMinX2 = 1;
      lMaxX2 = 4;
      lMinX3 = maxX3 - 4;
      lMaxX3 = maxX3 - 3;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
      break;
      //   BNE
   case BNE:
      lMinX1 = maxX1 - 4;
      lMaxX1 = maxX1 - 3;
      lMinX2 = maxX2 - 4;
      lMaxX2 = maxX2 - 1;
      lMinX3 = 1;
      lMaxX3 = 4;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      lMinX1 = maxX1 - 4;
      lMaxX1 = maxX1 - 1;
      lMinX2 = maxX2 - 4;
      lMaxX2 = maxX2 - 3;
      lMinX3 = 1;
      lMaxX3 = 4;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      lMinX1 = maxX1 - 4;
      lMaxX1 = maxX1 - 1;
      lMinX2 = maxX2 - 4;
      lMaxX2 = maxX2 - 1;
      lMinX3 = 3;
      lMaxX3 = 4;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      break;
      //   BNW
   case BNW:
      lMinX1 = 3;
      lMaxX1 = 4;
      lMinX2 = maxX2 - 4;
      lMaxX2 = maxX2 - 1;
      lMinX3 = 1;
      lMaxX3 = 4;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      lMinX1 = 1;
      lMaxX1 = 4;
      lMinX2 = maxX2 - 4;
      lMaxX2 = maxX2 - 3;
      lMinX3 = 1;
      lMaxX3 = 4;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      lMinX1 = 1;
      lMaxX1 = 4;
      lMinX2 = maxX2 - 4;
      lMaxX2 = maxX2 - 1;
      lMinX3 = 3;
      lMaxX3 = 4;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
      break;
      //   BSE
   case BSE:
      lMinX1 = maxX1 - 4;
      lMaxX1 = maxX1 - 3;
      lMinX2 = 1;
      lMaxX2 = 4;
      lMinX3 = 1;
      lMaxX3 = 4;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      lMinX1 = maxX1 - 4;
      lMaxX1 = maxX1 - 1;
      lMinX2 = 3;
      lMaxX2 = 4;
      lMinX3 = 1;
      lMaxX3 = 4;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      lMinX1 = maxX1 - 4;
      lMaxX1 = maxX1 - 1;
      lMinX2 = 1;
      lMaxX2 = 4;
      lMinX3 = 3;
      lMaxX3 = 4;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);
      break;
      //   BSW
   case BSW:
      lMinX1 = 3;
      lMaxX1 = 4;
      lMinX2 = 1;
      lMaxX2 = 4;
      lMinX3 = 1;
      lMaxX3 = 4;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      lMinX1 = 1;
      lMaxX1 = 4;
      lMinX2 = 3;
      lMaxX2 = 4;
      lMinX3 = 1;
      lMaxX3 = 4;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      lMinX1 = 1;
      lMaxX1 = 4;
      lMinX2 = 1;
      lMaxX2 = 4;
      lMinX3 = 3;
      lMaxX3 = 4;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      break;
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  >
void D3Q27ETCFOffVectorConnector< VectorTransmitter>::distributeReceiveVector(DistributionArray3DPtr fTo, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index)
{
   if (data.size() == 0) return;

   int ix1, ix2, ix3;
   for (ix3 = lMinX3; ix3 < lMaxX3; ix3++)
   {
      for (ix2 = lMinX2; ix2 < lMaxX2; ix2++)
      {
         for (ix1 = lMinX1; ix1 < lMaxX1; ix1++)
         {
            LBMReal icellC[27];
            this->readICellCfromData(data, index, icellC);
            iprocessor->writeINode(fTo, icellC, ix1, ix2, ix3);
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  >
void D3Q27ETCFOffVectorConnector< VectorTransmitter>::readICellCfromData(vector_type& data, int& index, LBMReal* icellC)
{
   for (int i = D3Q27System::STARTF; i < D3Q27System::ENDF + 1; i++)
   {
      icellC[i] = data[index++];
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  >
void D3Q27ETCFOffVectorConnector< VectorTransmitter>::getLocalMinMax(int& minX1, int& minX2, int& minX3, int& maxX1, int& maxX2, int& maxX3)
{
   using namespace D3Q27System;
   int TminX1 = minX1; int TminX2 = minX2; int TminX3 = minX3; int TmaxX1 = maxX1; int TmaxX2 = maxX2; int TmaxX3 = maxX3;

   if (block.lock()->hasInterpolationFlagCF(E))
   {
      if (maxX1 == TmaxX1) maxX1 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(W))
   {
      if (minX1 == TminX1) minX1 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(N))
   {
      if (maxX2 == TmaxX2)  maxX2 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(S))
   {
      if (minX2 == TminX2)  minX2 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(T))
   {
      if (maxX3 == TmaxX3)  maxX3 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(B))
   {
      if (minX3 == TminX3)  minX3 += 2;
   }

   //E-W-N-S
   if (block.lock()->hasInterpolationFlagCF(NE) && !block.lock()->hasInterpolationFlagCF(N) && !block.lock()->hasInterpolationFlagCF(E))
   {
      if (maxX1 == TmaxX1) maxX1 -= 2;
      if (maxX2 == TmaxX2) maxX2 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(SW) && !block.lock()->hasInterpolationFlagCF(W) && !block.lock()->hasInterpolationFlagCF(S))
   {
      if (minX1 == TminX1) minX1 += 2;
      if (minX2 == TminX2) minX2 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(SE) && !block.lock()->hasInterpolationFlagCF(E) && !block.lock()->hasInterpolationFlagCF(S))
   {
      if (maxX1 == TmaxX1) maxX1 -= 2;
      if (minX2 == TminX2) minX2 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(NW) && !block.lock()->hasInterpolationFlagCF(N) && !block.lock()->hasInterpolationFlagCF(W))
   {
      if (minX1 == TminX1) minX1 += 2;
      if (maxX2 == TmaxX2) maxX2 -= 2;
   }

   //	////T-B-E-W
   if (block.lock()->hasInterpolationFlagCF(TE) && !block.lock()->hasInterpolationFlagCF(E) && !block.lock()->hasInterpolationFlagCF(T))
   {
      if (maxX1 == TmaxX1) maxX1 -= 2;
      if (maxX3 == TmaxX3) maxX3 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(BW) && !block.lock()->hasInterpolationFlagCF(W) && !block.lock()->hasInterpolationFlagCF(B))
   {
      if (minX1 == TminX1) minX1 += 2;
      if (minX3 == TminX3) minX3 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(BE) && !block.lock()->hasInterpolationFlagCF(E) && !block.lock()->hasInterpolationFlagCF(B))
   {
      if (maxX1 == TmaxX1) maxX1 -= 2;
      if (minX3 == TminX3) minX3 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(TW) && !block.lock()->hasInterpolationFlagCF(W) && !block.lock()->hasInterpolationFlagCF(T))
   {
      if (minX1 == TminX1) minX1 += 2;
      if (maxX3 == TmaxX3) maxX3 -= 2;
   }


   ////T-B-N-S
   if (block.lock()->hasInterpolationFlagCF(TN) && !block.lock()->hasInterpolationFlagCF(N) && !block.lock()->hasInterpolationFlagCF(T))
   {
      if (maxX2 == TmaxX2) maxX2 -= 2;
      if (maxX3 == TmaxX3) maxX3 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(BS) && !block.lock()->hasInterpolationFlagCF(S) && !block.lock()->hasInterpolationFlagCF(B))
   {
      if (minX2 == TminX2) minX2 += 2;
      if (minX3 == TminX3) minX3 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(BN) && !block.lock()->hasInterpolationFlagCF(N) && !block.lock()->hasInterpolationFlagCF(B))
   {
      if (maxX2 == TmaxX2) maxX2 -= 2;
      if (minX3 == TminX3) minX3 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(TS) && !block.lock()->hasInterpolationFlagCF(S) && !block.lock()->hasInterpolationFlagCF(T))
   {
      if (minX2 == TminX2) minX2 += 2;
      if (maxX3 == TmaxX3) maxX3 -= 2;
   }

   //if (block.lock()->hasInterpolationFlagCF(D3Q27System::TNE)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::TE)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::TN)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::NE)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::T)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::N) && !block.lock()->hasInterpolationFlagCF(D3Q27System::E))
   //if (!block.lock()->hasInterpolationFlagCF(D3Q27System::TE)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::T)&& !block.lock()->hasInterpolationFlagCF(D3Q27System::E))
   //{
   //   if (maxX1==TmaxX1) maxX1 -= 2;
   //   if (maxX2==TmaxX2) maxX2 -= 2;
   //   if (maxX3==TmaxX3) maxX3 -= 2;
   //}
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFOffVectorConnector< VectorTransmitter>::getLocalMinMax(int& minX1, int& minX2, int& minX3, int& maxX1, int& maxX2, int& maxX3, CFconnectorType connType)
{
   using namespace D3Q27System;
   int TminX1 = minX1; int TminX2 = minX2; int TminX3 = minX3; int TmaxX1 = maxX1; int TmaxX2 = maxX2; int TmaxX3 = maxX3;

   if (block.lock()->hasInterpolationFlagCF(E))
   {
      if (maxX1 == TmaxX1) maxX1 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(W))
   {
      if (minX1 == TminX1) minX1 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(N))
   {
      if (maxX2 == TmaxX2)  maxX2 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(S))
   {
      if (minX2 == TminX2)  minX2 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(T))
   {
      if (maxX3 == TmaxX3)  maxX3 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(B))
   {
      if (minX3 == TminX3)  minX3 += 2;
   }

   //E-W-N-S
   if (block.lock()->hasInterpolationFlagCF(NE) && !block.lock()->hasInterpolationFlagCF(N) && !block.lock()->hasInterpolationFlagCF(E))
   {
      if (maxX1 == TmaxX1) maxX1 -= 2;
      if (maxX2 == TmaxX2) maxX2 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(SW) && !block.lock()->hasInterpolationFlagCF(W) && !block.lock()->hasInterpolationFlagCF(S))
   {
      if (minX1 == TminX1) minX1 += 2;
      if (minX2 == TminX2) minX2 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(SE) && !block.lock()->hasInterpolationFlagCF(E) && !block.lock()->hasInterpolationFlagCF(S))
   {
      if (maxX1 == TmaxX1) maxX1 -= 2;
      if (minX2 == TminX2) minX2 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(NW) && !block.lock()->hasInterpolationFlagCF(N) && !block.lock()->hasInterpolationFlagCF(W))
   {
      if (minX1 == TminX1) minX1 += 2;
      if (maxX2 == TmaxX2) maxX2 -= 2;
   }

   //	////T-B-E-W
   if (block.lock()->hasInterpolationFlagCF(TE) && !block.lock()->hasInterpolationFlagCF(E) && !block.lock()->hasInterpolationFlagCF(T))
   {
      if (maxX1 == TmaxX1) maxX1 -= 2;
      if (maxX3 == TmaxX3) maxX3 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(BW) && !block.lock()->hasInterpolationFlagCF(W) && !block.lock()->hasInterpolationFlagCF(B))
   {
      if (minX1 == TminX1) minX1 += 2;
      if (minX3 == TminX3) minX3 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(BE) && !block.lock()->hasInterpolationFlagCF(E) && !block.lock()->hasInterpolationFlagCF(B))
   {
      if (maxX1 == TmaxX1) maxX1 -= 2;
      if (minX3 == TminX3) minX3 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(TW) && !block.lock()->hasInterpolationFlagCF(W) && !block.lock()->hasInterpolationFlagCF(T))
   {
      if (minX1 == TminX1) minX1 += 2;
      if (maxX3 == TmaxX3) maxX3 -= 2;
   }


   ////T-B-N-S
   if (block.lock()->hasInterpolationFlagCF(TN) && !block.lock()->hasInterpolationFlagCF(N) && !block.lock()->hasInterpolationFlagCF(T))
   {
      if (maxX2 == TmaxX2) maxX2 -= 2;
      if (maxX3 == TmaxX3) maxX3 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(BS) && !block.lock()->hasInterpolationFlagCF(S) && !block.lock()->hasInterpolationFlagCF(B))
   {
      if (minX2 == TminX2) minX2 += 2;
      if (minX3 == TminX3) minX3 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(BN) && !block.lock()->hasInterpolationFlagCF(N) && !block.lock()->hasInterpolationFlagCF(B))
   {
      if (maxX2 == TmaxX2) maxX2 -= 2;
      if (minX3 == TminX3) minX3 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(TS) && !block.lock()->hasInterpolationFlagCF(S) && !block.lock()->hasInterpolationFlagCF(T))
   {
      if (minX2 == TminX2) minX2 += 2;
      if (maxX3 == TmaxX3) maxX3 -= 2;
   }

   //if (block.lock()->hasInterpolationFlagCF(D3Q27System::TNE)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::TE)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::TN)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::NE)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::T)&&!block.lock()->hasInterpolationFlagCF(D3Q27System::N) && !block.lock()->hasInterpolationFlagCF(D3Q27System::E))
   //{
   //   if (maxX1==TmaxX1) maxX1 -= 2;
   //   if (maxX2==TmaxX2) maxX2 -= 2;
   //   if (maxX3==TmaxX3) maxX3 -= 2;
   //}
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETCFOffVectorConnector< VectorTransmitter>::findCFnodes()
{
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

   using namespace D3Q27System;
   if (block.lock()->hasInterpolationFlagCF(W))
   {
      lMinX1 = 1;
      lMaxX1 = lMinX1 + 1;

      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
      findCFnodes(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX2, maxX2, true, lMinX2, lMaxX2, false);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
      findCFnodes(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvOd, indexEvOd);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
      getLocalMinMax(minX3, maxX3, true, lMinX3, lMaxX3, false);
      findCFnodes(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);

      getLocalMinMax(minX2, maxX2, false, lMinX2, lMaxX2, false);
      getLocalMinMax(minX3, maxX3, false, lMinX3, lMaxX3, false);
      findCFnodes(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdOd, indexOdOd);
   }
   if (block.lock()->hasInterpolationFlagCF(TN) && !block.lock()->hasInterpolationFlagCF(N) && !block.lock()->hasInterpolationFlagCF(T))
   {
      lMinX2 = maxX2 - 3;
      lMaxX2 = lMinX2 + 1;
      lMinX3 = maxX3 - 3;
      lMaxX3 = lMinX3 + 1;

      getLocalMinMax(minX1 + 1, maxX1, true, lMinX1, lMaxX1, false);
      findCFnodes(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataEvEv, indexEvEv);

      getLocalMinMax(minX1 + 1, maxX1, false, lMinX1, lMaxX1, false);
      findCFnodes(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, dataOdEv, indexOdEv);
   }
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  >
void D3Q27ETCFOffVectorConnector< VectorTransmitter>::findCFnodes(DistributionArray3DPtr fFrom, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index)
{
   if (data.size() == 0) return;
   int ix1, ix2, ix3;
   LBMReal xoff, yoff, zoff;
   BCArray3D& bcArray = block.lock()->getKernel()->getBCProcessor()->getBCArray();

   for (ix3 = lMinX3; ix3 < lMaxX3; ix3++)
   {
      for (ix2 = lMinX2; ix2 < lMaxX2; ix2++)
      {
         for (ix1 = lMinX1; ix1 < lMaxX1; ix1++)
         {
            D3Q27ICell icellC;
            D3Q27ICell icellF;

            int howManySolids = iprocessor->iCellHowManySolids(bcArray, ix1, ix2, ix3);

            if (howManySolids == 0 || howManySolids == 8)
            {
               iprocessor->readICell(fFrom, icellC, ix1, ix2, ix3);
               xoff = 0.0;
               yoff = 0.0;
               zoff = 0.0;
            }
            else
            {
               if (!iprocessor->findNeighborICell(bcArray, fFrom, icellC, bMaxX1, bMaxX2, bMaxX3, ix1, ix2, ix3, xoff, yoff, zoff))
               {
                  std::string err = "For " + block.lock()->toString() + " x1=" + UbSystem::toString(ix1) + ", x2=" + UbSystem::toString(ix2) + ", x3=" + UbSystem::toString(ix3) +
                     " interpolation is not implemented for other direction" +
                     " by using in: " + (std::string)typeid(*this).name() +
                     " or maybe you have a solid on the block boundary";
                  UBLOG(logINFO, err);
                  //UB_THROW(UbException(UB_EXARGS, err));
               }
            }

            iprocessor->interpolateCoarseToFine(icellC, icellF, xoff, yoff, zoff);
            this->writeICellFtoData(data, index, icellF);
            //for (int iix3 = ix3; iix3<=ix3+1; iix3++)
            //{
            //   for (int iix2 = ix2; iix2<=ix2+1; iix2++)
            //   {
            //      for (int iix1 = ix1; iix1<=ix1+1; iix1++)
            //      {
            //         bcArray.setInterfaceCF(iix1, iix2, iix3);
            //      }
            //   }
            //}

         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
double D3Q27ETCFOffVectorConnector<VectorTransmitter>::getSendRecieveTime()
{
   return 0;
}

#endif 
