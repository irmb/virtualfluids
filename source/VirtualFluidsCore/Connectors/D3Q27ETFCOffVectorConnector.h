/**
* @file D3Q27ETFCOffVectorConnector.h
* @class D3Q27ETFCVectorConnector
* @brief Interpolation from fine level to coarse.
* @author Kostyantyn Kucher and Ehsan Fard
* @date 08.06.2011
*/
#ifndef D3Q27ETFCOffVectorConnector_H
#define D3Q27ETFCOffVectorConnector_H

#include <vector>

#include "basics/transmitter/TbTransmitter.h"
#include "Block3DConnector.h"
#include "D3Q27System.h"
#include "Block3D.h"
#include "Grid3D.h"
#include "LBMKernelETD3Q27.h"
#include "D3Q27InterpolationProcessor.h"
#include "MathUtil.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>


class Block3D;

enum CFconnectorType {EvenOddNW, EvenEvenSW, OddEvenSE, OddOddNE};

//daten werden in einen vector (dieser befindet sich im transmitter) kopiert
//der vector wird via transmitter uebertragen
//transmitter kann ein lokal, MPI, RCG, CTL oder was auch immer fuer ein
//transmitter sein, der von Transmitter abgeleitet ist ;-)

template< typename VectorTransmitter >
class D3Q27ETFCOffVectorConnector : public Block3DConnector
{
public:

protected:
	typedef typename VectorTransmitter::value_type  vector_type;
	typedef boost::shared_ptr< VectorTransmitter > VectorTransmitterPtr;
public:
   D3Q27ETFCOffVectorConnector(Block3DPtr block, VectorTransmitterPtr sender, VectorTransmitterPtr receiver, int sendDir, 
      D3Q27InterpolationProcessorPtr iprocessor, CFconnectorType connType);

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

	void writeICellCtoData(vector_type& data, int& index, LBMReal* icellC);
	void writeNodeToVector(vector_type& data, int& index, LBMReal* inode);
	void getLocalMinMax(int& minX1, int& minX2, int& minX3, int& maxX1, int& maxX2, int& maxX3);
   void getLocalMinMax(int& minX1, int& minX2, int& minX3, int& maxX1, int& maxX2, int& maxX3, CFconnectorType connType);
	void getLocalMinMaxCF(int gMax, int& lMin, int& lMax);
	void fillSendVector(DistributionArray3DPtr fFrom, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index);

	void distributeReceiveVector(DistributionArray3DPtr fTo, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index);
	void readICellFfromData(vector_type& data, int& index, D3Q27ICell& icellF);
	void readNodeFromVector(vector_type& data, int& index, LBMReal* inode);
	void getLocalOffsets(const int& gMax, int& oMin);
	void getLocalMins(int& minX1, int& minX2, int& minX3, const int& oMinX1, const int& oMinX2, const int& oMinX3);

	int bMaxX1, bMaxX2, bMaxX3;
};
////////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
D3Q27ETFCOffVectorConnector<VectorTransmitter>::D3Q27ETFCOffVectorConnector(Block3DPtr block, VectorTransmitterPtr sender, 
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
	if( !(   sendDir==D3Q27System::E  || sendDir==D3Q27System::W  || sendDir==D3Q27System::N  || sendDir==D3Q27System::S  || sendDir==D3Q27System::T || sendDir==D3Q27System::B 
		||  sendDir==D3Q27System::NE || sendDir==D3Q27System::SW || sendDir==D3Q27System::SE || sendDir==D3Q27System::NW
		||  sendDir==D3Q27System::TE || sendDir==D3Q27System::BW ||  sendDir==D3Q27System::BE || sendDir==D3Q27System::TW
		||  sendDir==D3Q27System::TN || sendDir==D3Q27System::BS ||  sendDir==D3Q27System::BN || sendDir==D3Q27System::TS 

		||  sendDir==D3Q27System::TNE || sendDir==D3Q27System::TNW ||  sendDir==D3Q27System::TSE || sendDir==D3Q27System::TSW
		||  sendDir==D3Q27System::BNE || sendDir==D3Q27System::BNW ||  sendDir==D3Q27System::BSE || sendDir==D3Q27System::BSW 
		
		) )
	{
		throw UbException(UB_EXARGS,"invalid constructor for this direction");
	}
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
bool D3Q27ETFCOffVectorConnector<VectorTransmitter>::isLocalConnector()
{ 
	return !this->isRemoteConnector(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
bool D3Q27ETFCOffVectorConnector<VectorTransmitter>::isRemoteConnector() 
{ 
	return sender->isRemoteTransmitter()  ||  receiver->isRemoteTransmitter();
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETFCOffVectorConnector<VectorTransmitter>::sendTransmitterDataSize()  
{ 
	if(sender)
	{
		UBLOG(logDEBUG5, "D3Q27ETFCOffVectorConnector<VectorTransmitter>::sendTransmitterDataSize()"<<block.lock()->toString()+"sendDir="<<sendDir);
		sender->sendDataSize(); 
	}
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETFCOffVectorConnector<VectorTransmitter>::receiveTransmitterDataSize()
{ 
	if(receiver)
	{
		UBLOG(logDEBUG5, "D3Q27ETFCOffVectorConnector<VectorTransmitter>::receiveTransmitterDataSize()"<<block.lock()->toString()<<"sendDir="<<sendDir);
		receiver->receiveDataSize();
	}
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETFCOffVectorConnector<VectorTransmitter>::prepareForSend()
{ 
	if(sender) sender->prepareForSend(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETFCOffVectorConnector<VectorTransmitter>::sendVectors()     
{ 
	if(sender) sender->sendData();
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETFCOffVectorConnector<VectorTransmitter>::prepareForReceive()     
{ 
	if(receiver) receiver->prepareForReceive(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETFCOffVectorConnector<VectorTransmitter>::receiveVectors() 
{ 
	if(receiver) receiver->receiveData(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETFCOffVectorConnector<VectorTransmitter>::init()
{
	using namespace D3Q27System;

	bMaxX1 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX1();
	bMaxX2 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX2();
	bMaxX3 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX3();

	int       sendSize  = 0;
	LBMReal initValue = -999.0;

	int sendDataPerNode = 27/*f*/;
	int iCellSize = 1; //size of interpolation cell

	switch(this->sendDir)
	{		                  
	case E : case W : sendSize = (bMaxX2-1)/2*(bMaxX3-1)/2*sendDataPerNode*iCellSize; break; 
	case N : case S : sendSize = (bMaxX1-1)/2*(bMaxX3-1)/2*sendDataPerNode*iCellSize; break; 
	case T : case B : sendSize = (bMaxX1-1)/2*(bMaxX2-1)/2*sendDataPerNode*iCellSize; break; 		  
	case NE : case SW :case SE : case NW : sendSize = (3*bMaxX3-3)*sendDataPerNode*iCellSize; break; // buffer overhead, should be (3*bMaxX3-6) for even bMax3		
	case TE : case BW :case BE : case TW : sendSize = (3*bMaxX2-3)*sendDataPerNode*iCellSize; break; 
	case TN : case BS :case BN : case TS : sendSize = (3*bMaxX1-3)*sendDataPerNode*iCellSize; break;	
   case TNE: case TNW:case TSE: case TSW:case BNE: case BNW:case BSE: case BSW: sendSize = 3*(3*bMaxX1-3)*sendDataPerNode*iCellSize; break;
	default: throw UbException(UB_EXARGS,"direction not allowed in this constructor");
	}
	sender->getData().resize(sendSize, initValue);
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETFCOffVectorConnector< VectorTransmitter>::fillSendVectors() 
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

   ///////////////////////////////////////
   ///DEBUG
#ifdef _DEBUG
   if (block.lock()->getGlobalID() == 2558)
   {
      int test = 0;
   }
#endif
   //////////////

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
      ///////////////////////////////////////
      ///DEBUG
#ifdef _DEBUG
      if (block.lock()->getGlobalID() == 2516)
      {
         int test = 0;
      }
#endif
      //////////////
		//getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, none);
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

	//	////N-S-E-W
	case NE: 
		getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
		getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
		lMinX1 = maxX1-7;
		lMaxX1 = lMinX1 +5;
		lMinX2 = maxX2-7;
		lMaxX2 = lMinX2 + 1;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		
		lMinX1 = maxX1-7;
		lMaxX1 = lMinX1 + 1;
		lMinX2 = maxX2-7;
		lMaxX2 = lMinX2 + 5;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;
	case SW: 
		
		getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
		getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
		lMinX1 = 1;
		lMaxX1 = lMinX1 + 5;
		lMinX2 = 5;
		lMaxX2 = lMinX2 + 1;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

		lMinX1 = 5;
		lMaxX1 = lMinX1 + 1;
		lMinX2 = 1;
		lMaxX2 = lMinX2 + 5;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;

	case SE: 
		getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
		getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
		lMinX1 = maxX1-7;
		lMaxX1 = lMinX1 +5;
		lMinX2 = 5;
		lMaxX2 = lMinX2 + 1;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		
		lMinX1 = maxX1-7;
		lMaxX1 = lMinX1 + 1;
		lMinX2 = 1;
		lMaxX2 = lMinX2 + 5;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

		break;

	case NW: 
		getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
		getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
		lMinX1 = 1;
		lMaxX1 = lMinX1 + 5;
		lMinX2 = maxX2-7;
		lMaxX2 = lMinX2 + 1;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

		lMinX1 = 5;
		lMaxX1 = lMinX1 + 1;
		lMinX2 = maxX2-7;
		lMaxX2 = lMinX2 + 5;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;
//////T-B-E-W
	case TE: 
		getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
		getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
		lMinX1 = maxX1-7;
		lMaxX1 = lMinX1 +5;
		lMinX3 = maxX3-7;
		lMaxX3 = lMinX3 + 1;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

		lMinX1 = maxX1-7;
		lMaxX1 = lMinX1 + 1;
		lMinX3 = maxX3-7;
		lMaxX3 = lMinX3 + 5;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;

	case BW: 
		getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
		getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
		lMinX1 = 1;
		lMaxX1 = lMinX1 + 5;
		lMinX3 = 5;
		lMaxX3 = lMinX3 + 1;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

		lMinX1 = 5;
		lMaxX1 = lMinX1 + 1;
		lMinX3 = 1;
		lMaxX3 = lMinX3 + 5;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;

	case BE: 
		getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
		getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
		lMinX1 = maxX1-7;
		lMaxX1 = lMinX1 +5;
		lMinX3 = 5;
		lMaxX3 = lMinX3 + 1;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

		lMinX1 = maxX1-7;
		lMaxX1 = lMinX1 + 1;
		lMinX3 = 1;
		lMaxX3 = lMinX3 + 5;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;

	case TW: 
		getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
		getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
		lMinX1 = 1;
		lMaxX1 = lMinX1 + 5;
		lMinX3 = maxX3-7;
		lMaxX3 = lMinX3 + 1;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

		lMinX1 = 5;
		lMaxX1 = lMinX1 + 1;
		lMinX3 = maxX3-7;
		lMaxX3 = lMinX3 + 5;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;
///////////////T-B-N-S
//
	case TN: 
		getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
		getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
		lMinX2 = maxX2-7;
		lMaxX2 = lMinX2 + 5;
		lMinX3 = maxX3-7;
		lMaxX3 = lMinX3 + 1;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

		lMinX2 = maxX2-7;
		lMaxX2 = lMinX2 + 1;
		lMinX3 = maxX3-7;
		lMaxX3 = lMinX3 + 5;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;

	case BS: 
		getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
		getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
		lMinX2 = 1;
		lMaxX2 = lMinX2 + 5;
		lMinX3 = 5;
		lMaxX3 = lMinX3 + 1;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

		lMinX2 = 5;
		lMaxX2 = lMinX2 + 1;
		lMinX3 = 1;
		lMaxX3 = lMinX3 + 5;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;

	case BN: 
		getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
		getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
		lMinX2 = maxX2-7;
		lMaxX2 = lMinX2 + 5;
		lMinX3 = 5;
		lMaxX3 = lMinX3 + 1;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

		lMinX2 = maxX2-7;
		lMaxX2 = lMinX2 + 1;
		lMinX3 = 1;
		lMaxX3 = lMinX3 + 5;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;

	case TS: 
		getLocalMinMax(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3);
		getLocalMins(lMinX1, lMinX2, lMinX3, oMinX1, oMinX2, oMinX3);
		lMinX2 = 1;
		lMaxX2 = lMinX2 + 5;
		lMinX3 = maxX3-7;
		lMaxX3 = lMinX3 + 1;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

		lMinX2 = 5;
		lMaxX2 = lMinX2 + 1;
		lMinX3 = maxX3-7;
		lMaxX3 = lMinX3 + 5;
		fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;
   
   //TNE
   case TNE:
      lMinX1 = maxX1-7;
      lMaxX1 = maxX1-6;
      lMinX2 = maxX2-7;
      lMaxX2 = maxX2-2;
      lMinX3 = maxX3-7;
      lMaxX3 = maxX3-2;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      lMinX1 = maxX1-7;
      lMaxX1 = maxX1-2;
      lMinX2 = maxX2-7;
      lMaxX2 = maxX2-6;
      lMinX3 = maxX3-7;
      lMaxX3 = maxX3-2;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      lMinX1 = maxX1-7;
      lMaxX1 = maxX1-2;
      lMinX2 = maxX2-7;
      lMaxX2 = maxX2-2;
      lMinX3 = maxX3-7;
      lMaxX3 = maxX3-6;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;


   //TNW
   case TNW:
      lMinX1 = 5;
      lMaxX1 = 6;
      lMinX2 = maxX2-7;
      lMaxX2 = maxX2-2;
      lMinX3 = maxX3-7;
      lMaxX3 = maxX3-2;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      lMinX1 = 1;
      lMaxX1 = 6;
      lMinX2 = maxX2-7;
      lMaxX2 = maxX2-6;
      lMinX3 = maxX3-7;
      lMaxX3 = maxX3-2;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      lMinX1 = 1;
      lMaxX1 = 6;
      lMinX2 = maxX2-7;
      lMaxX2 = maxX2-2;
      lMinX3 = maxX3-7;
      lMaxX3 = maxX3-6;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;

      break;
   
   //      TSE
   case TSE:
      lMinX1 = maxX1-7;
      lMaxX1 = maxX1-6;
      lMinX2 = 1;
      lMaxX2 = 6;
      lMinX3 = maxX3-7;
      lMaxX3 = maxX3-2;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      lMinX1 = maxX1-7;
      lMaxX1 = maxX1-2;
      lMinX2 = 5;
      lMaxX2 = 6;
      lMinX3 = maxX3-7;
      lMaxX3 = maxX3-2;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      lMinX1 = maxX1-7;
      lMaxX1 = maxX1-2;
      lMinX2 = 1;
      lMaxX2 = 6;
      lMinX3 = maxX3-7;
      lMaxX3 = maxX3-6;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      break;
   //      TSW
   case TSW:
      lMinX1 = 5;
      lMaxX1 = 6;
      lMinX2 = 1;
      lMaxX2 = 6;
      lMinX3 = maxX3-7;
      lMaxX3 = maxX3-2;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      lMinX1 = 1;
      lMaxX1 = 6;
      lMinX2 = 5;
      lMaxX2 = 6;
      lMinX3 = maxX3-7;
      lMaxX3 = maxX3-2;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      lMinX1 = 1;
      lMaxX1 = 6;
      lMinX2 = 1;
      lMaxX2 = 6;
      lMinX3 = maxX3-7;
      lMaxX3 = maxX3-6;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      break;
   //      BNE
   case BNE:
      lMinX1 = maxX1-7;
      lMaxX1 = maxX1-6;
      lMinX2 = maxX2-7;
      lMaxX2 = maxX2-2;
      lMinX3 = 1;
      lMaxX3 = 6;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      lMinX1 = maxX1-7;
      lMaxX1 = maxX1-2;
      lMinX2 = maxX2-7;
      lMaxX2 = maxX2-6;
      lMinX3 = 1;
      lMaxX3 = 6;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      lMinX1 = maxX1-7;
      lMaxX1 = maxX1-2;
      lMinX2 = maxX2-7;
      lMaxX2 = maxX2-2;
      lMinX3 = 5;
      lMaxX3 = 6;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      break;
   //      BNW
   case BNW:
      lMinX1 = 5;
      lMaxX1 = 6;
      lMinX2 = maxX2-7;
      lMaxX2 = maxX2-2;
      lMinX3 = 1;
      lMaxX3 = 6;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      lMinX1 = 1;
      lMaxX1 = 6;
      lMinX2 = maxX2-7;
      lMaxX2 = maxX2-6;
      lMinX3 = 1;
      lMaxX3 = 6;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      lMinX1 = 1;
      lMaxX1 = 6;
      lMinX2 = maxX2-7;
      lMaxX2 = maxX2-2;
      lMinX3 = 5;
      lMaxX3 = 6;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      break;


   //      BSE
   case BSE:
      lMinX1 = maxX1-7;
      lMaxX1 = maxX1-6;
      lMinX2 = 1;
      lMaxX2 = 6;
      lMinX3 = 1;
      lMaxX3 = 6;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      lMinX1 = maxX1-7;
      lMaxX1 = maxX1-2;
      lMinX2 = 5;
      lMaxX2 = 6;
      lMinX3 = 1;
      lMaxX3 = 6;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      lMinX1 = maxX1-7;
      lMaxX1 = maxX1-2;
      lMinX2 = 1;
      lMaxX2 = 6;
      lMinX3 = 5;
      lMaxX3 = 6;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      break;

   //BSW
   case BSW:
      lMinX1 = 5;
      lMaxX1 = 6;
      lMinX2 = 1;
      lMaxX2 = 6;
      lMinX3 = 1;
      lMaxX3 = 6;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      
      lMinX1 = 1;
      lMaxX1 = 6;
      lMinX2 = 5;
      lMaxX2 = 6;
      lMinX3 = 1;
      lMaxX3 = 6;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      lMinX1 = 1;
      lMaxX1 = 6;
      lMinX2 = 1;
      lMaxX2 = 6;
      lMinX3 = 5;
      lMaxX3 = 6;
      fillSendVector(fFrom, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);

      break;
	}
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCOffVectorConnector< VectorTransmitter>::fillSendVector(DistributionArray3DPtr fFrom, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index)
{
	int ix1, ix2, ix3;
	LBMReal xoff, yoff, zoff;
	BCArray3D<D3Q27BoundaryCondition>& bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(block.lock()->getKernel()->getBCProcessor())->getBCArray();

	for (ix3=lMinX3; ix3<lMaxX3; ix3+=2)
	{
		for (ix2=lMinX2; ix2<lMaxX2; ix2+=2)
		{
			for (ix1=lMinX1; ix1<lMaxX1; ix1+=2)
			{
				LBMReal icellC[27];
				D3Q27ICell icellF;

				int howManySolids= iprocessor->iCellHowManySolids(bcArray, ix1, ix2, ix3);

				if(howManySolids == 0 || howManySolids == 8)
				{
					iprocessor->readICell(fFrom, icellF, ix1, ix2, ix3);
					xoff=0.0; 
					yoff=0.0;
					zoff=0.0;
				}
				else
				{
					if(!iprocessor->findNeighborICell(bcArray, fFrom, icellF, bMaxX1, bMaxX2, bMaxX3, ix1, ix2, ix3, xoff, yoff, zoff))
					{
						std::string err = "For "+block.lock()->toString()+" x1="+UbSystem::toString(ix1)+", x2=" + UbSystem::toString(ix2)+", x3=" + UbSystem::toString(ix3)+
							" interpolation is not implemented for other direction"+
							" by using in: "+(std::string)typeid(*this).name()+ 
							" or maybe you have a solid on the block boundary";
						UB_THROW(UbException(UB_EXARGS, err));
					}
				}

				iprocessor->interpolateFineToCoarse(icellF, icellC, xoff, yoff, zoff);
				this->writeICellCtoData(data, index, icellC);
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCOffVectorConnector< VectorTransmitter>::writeICellCtoData(vector_type& data, int& index, LBMReal* icellC) 
{
	for (int i = D3Q27System::STARTF; i < D3Q27System::ENDF+1; i++)
	{
		data[index++] = icellC[i];
	}
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCOffVectorConnector< VectorTransmitter>::getLocalMinMaxCF(int gMax, int& lMin, int& lMax)
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
void D3Q27ETFCOffVectorConnector< VectorTransmitter>::distributeReceiveVectors() 
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
      ///////////////////////////////////////
      ///DEBUG
      //if (block.lock()->getGlobalID() == 2554)
      //{
      //   int test = 0;
      //}
      //////////////
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

		/////E-W-N-S
	case NE: 
		lMinX1 = maxX1-4;
		lMaxX1 = lMinX1 + 3;
		lMinX2 = maxX2-4;
		lMaxX2 = lMinX2 + 3;
		getLocalMinMaxCF(maxX3, lMinX3, lMaxX3);
		distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;

	case SW: 
		lMinX1 = 0;
		lMaxX1 = lMinX1 + 3;
		lMinX2 = 0;
		lMaxX2 = lMinX2 + 3;
		getLocalMinMaxCF(maxX3, lMinX3, lMaxX3);
		distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;

	case SE: 
		lMinX1 = maxX1-4;
		lMaxX1 = lMinX1 + 3;
		lMinX2 = 0;
		lMaxX2 = lMinX2 + 3;
		getLocalMinMaxCF(maxX3, lMinX3, lMaxX3);
		distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;

	case NW: 
		lMinX1 = 0;
		lMaxX1 = lMinX1 + 3;
		lMinX2 = maxX2-4;
		lMaxX2 = lMinX2 + 3;
		getLocalMinMaxCF(maxX3, lMinX3, lMaxX3);
		distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;
	//	
	//	/////T-B-E-W
	case TE:
		lMinX1 = maxX1-4;
		lMaxX1 = lMinX1 + 3;
		lMinX3 = maxX3-4;
		lMaxX3 = lMinX3 + 3;
		getLocalMinMaxCF(maxX2, lMinX2, lMaxX2);
		distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;
	
	case BW:
		lMinX1 = 0;
		lMaxX1 = lMinX1 + 3;
		lMinX3 = 0;
		lMaxX3 = lMinX3 + 3;
		getLocalMinMaxCF(maxX2, lMinX2, lMaxX2);
		distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;
	
	case BE:
		lMinX1 = maxX1-4;
		lMaxX1 = lMinX1 + 3;
		lMinX3 = 0;
		lMaxX3 = lMinX3 + 3;
		getLocalMinMaxCF(maxX2, lMinX2, lMaxX2);
		distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;

	case TW:
		lMinX1 = 0;
		lMaxX1 = lMinX1 + 3;
		lMinX3 = maxX3-4;
		lMaxX3 = lMinX3 + 3;
		getLocalMinMaxCF(maxX2, lMinX2, lMaxX2);
		distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;

	//	////////////////T-B-N-S
	//	
	case TN:
		lMinX2 = maxX2-4;
		lMaxX2 = lMinX2 + 3;
		lMinX3 = maxX3-4;
		lMaxX3 = lMinX3 + 3;
		getLocalMinMaxCF(maxX1, lMinX1, lMaxX1);
		distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;

	case BS:
		lMinX2 = 0;
		lMaxX2 = lMinX2 + 3;
		lMinX3 = 0;
		lMaxX3 = lMinX3 + 3;
		getLocalMinMaxCF(maxX1, lMinX1, lMaxX1);
		distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;

	case BN:
		lMinX2 = maxX2-4;
		lMaxX2 = lMinX2 + 3;
		lMinX3 = 0;
		lMaxX3 = lMinX3 + 3;
		getLocalMinMaxCF(maxX1, lMinX1, lMaxX1);
		distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;

	case TS:
		lMinX2 = 0;
		lMaxX2 = lMinX2 + 3;
		lMinX3 = maxX3-4;
		lMaxX3 = lMinX3 + 3;
		getLocalMinMaxCF(maxX1, lMinX1, lMaxX1);
		distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
		break;

   //   //TNE
   case TNE:
      lMinX1 = maxX1-4;
      lMaxX1 = maxX1-1;
      lMinX2 = maxX2-4;
      lMaxX2 = maxX2-1;
      lMinX3 = maxX3-4;
      lMaxX3 = maxX3-1;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
      //   TNW
   case TNW:
      lMinX1 = 0;
      lMaxX1 = 3;
      lMinX2 = maxX2-4;
      lMaxX2 = maxX2-1;
      lMinX3 = maxX3-4;
      lMaxX3 = maxX3-1;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
      //   TSE
   case TSE:
      lMinX1 = maxX1-4;
      lMaxX1 = maxX1-1;
      lMinX2 = 0;
      lMaxX2 = 3;
      lMinX3 = maxX3-4;
      lMaxX3 = maxX3-1;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
      //   TSW
   case TSW:
      lMinX1 = 0;
      lMaxX1 = 3;
      lMinX2 = 0;
      lMaxX2 = 3;
      lMinX3 = maxX3-4;
      lMaxX3 = maxX3-1;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
      //   BNE
   case BNE:
      lMinX1 = maxX1-4;
      lMaxX1 = maxX1-1;
      lMinX2 = maxX2-4;
      lMaxX2 = maxX2-1;
      lMinX3 = 0;
      lMaxX3 = 3;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
      //   BNW
   case BNW:
      lMinX1 = 0;
      lMaxX1 = 3;
      lMinX2 = maxX2-4;
      lMaxX2 = maxX2-1;
      lMinX3 = 0;
      lMaxX3 = 3;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
      //   BSE
   case BSE:
      lMinX1 = maxX1-4;
      lMaxX1 = maxX1-1;
      lMinX2 = 0;
      lMaxX2 = 3;
      lMinX3 = 0;
      lMaxX3 = 3;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
         //BSW
   case BSW:
      lMinX1 = 0;
      lMaxX1 = 3;
      lMinX2 = 0;
      lMaxX2 = 3;
      lMinX3 = 0;
      lMaxX3 = 3;
      distributeReceiveVector(fTo, lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, data, index);
      break;
	}
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCOffVectorConnector< VectorTransmitter>::distributeReceiveVector(DistributionArray3DPtr fTo, const int& lMinX1, const int& lMinX2, const int& lMinX3, const int& lMaxX1, const int& lMaxX2, const int& lMaxX3, vector_type& data, int& index)
{
	if(data.size() == 0) return; 

	int ix1, ix2, ix3;
	for (ix3=lMinX3; ix3<lMaxX3; ix3+=2)
	{
		for (ix2=lMinX2; ix2<lMaxX2; ix2+=2)
		{
			for (ix1=lMinX1; ix1<lMaxX1; ix1+=2)
			{
				D3Q27ICell icellF;
				this->readICellFfromData(data, index, icellF);
				iprocessor->writeICell(fTo, icellF, ix1, ix2, ix3);
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCOffVectorConnector< VectorTransmitter>::readICellFfromData(vector_type& data, int& index, D3Q27ICell& icellF) 
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
void D3Q27ETFCOffVectorConnector< VectorTransmitter>::readNodeFromVector(vector_type& data, int& index, LBMReal* inode)
{
	for (int i = D3Q27System::STARTF; i < D3Q27System::ENDF+1; i++)
	{
		inode[i] = data[index++];
	}
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCOffVectorConnector< VectorTransmitter>::getLocalMinMax(int& minX1, int& minX2, int& minX3, int& maxX1, int& maxX2, int& maxX3)
{
	using namespace D3Q27System;
    int TminX1=minX1; int TminX2=minX2; int TminX3=minX3; int TmaxX1=maxX1; int TmaxX2=maxX2; int TmaxX3=maxX3;

	if(block.lock()->hasInterpolationFlagFC(E))
	{
		if (maxX1==TmaxX1) maxX1 -= 3;
	}
	if(block.lock()->hasInterpolationFlagFC(W))
	{
		if (minX1==TminX1) minX1 += 4;
	}
	if(block.lock()->hasInterpolationFlagFC(N))
	{
		if (maxX2==TmaxX2) maxX2 -= 3;
	}
	if(block.lock()->hasInterpolationFlagFC(S))
	{
		if (minX2==TminX2) minX2 += 4;
	}
	if(block.lock()->hasInterpolationFlagFC(T))
	{
		if (maxX3==TmaxX3) maxX3 -= 3;
	}
	if(block.lock()->hasInterpolationFlagFC(B))
	{
		if (minX3==TminX3) minX3 += 4;
	}

	////////////
	/////E-W-N-S
	if(block.lock()->hasInterpolationFlagFC(NE)&& !block.lock()->hasInterpolationFlagFC(N) && !block.lock()->hasInterpolationFlagFC(E))
	{
		if (maxX1==TmaxX1) maxX1 -= 3;
		if (maxX2==TmaxX2) maxX2 -= 3;
	}
	if(block.lock()->hasInterpolationFlagFC(SW)&& !block.lock()->hasInterpolationFlagFC(W) && !block.lock()->hasInterpolationFlagFC(S))
	{
		if (minX1==TminX1) minX1 += 4;
		if (minX2==TminX2) minX2 += 4;
	}
	if(block.lock()->hasInterpolationFlagFC(SE)&& !block.lock()->hasInterpolationFlagFC(E) && !block.lock()->hasInterpolationFlagFC(S))
	{
		if (maxX1==TmaxX1) maxX1 -= 3;
		if (minX2==TminX2) minX2 += 4;
	}
	if(block.lock()->hasInterpolationFlagFC(NW)&& !block.lock()->hasInterpolationFlagFC(N) && !block.lock()->hasInterpolationFlagFC(W))
	{
		if (minX1==TminX1) minX1 += 4;
		if (maxX2==TmaxX2) maxX2 -= 3;
	}

	//////T-B-E-W
	if(block.lock()->hasInterpolationFlagFC(TE) && !block.lock()->hasInterpolationFlagFC(E) && !block.lock()->hasInterpolationFlagFC(T))
	{
		if (maxX1==TmaxX1) maxX1 -= 3;
		if (maxX3==TmaxX3) maxX3 -= 3;
	}
	if(block.lock()->hasInterpolationFlagFC(BW)&& !block.lock()->hasInterpolationFlagFC(W) && !block.lock()->hasInterpolationFlagFC(B))
	{
		if (minX1==TminX1) minX1 += 4;
		if (minX3==TminX3) minX3 += 4;
	}
	if(block.lock()->hasInterpolationFlagFC(BE)&& !block.lock()->hasInterpolationFlagFC(E) && !block.lock()->hasInterpolationFlagFC(B))
	{
		if (maxX1==TmaxX1) maxX1 -= 3;
		if (minX3==TminX3) minX3 += 4;
	}
	if(block.lock()->hasInterpolationFlagFC(TW)&& !block.lock()->hasInterpolationFlagFC(W) && !block.lock()->hasInterpolationFlagFC(T))
	{
		if (minX1==TminX1) minX1 += 4;
		if (maxX3==TmaxX3) maxX3 -= 3;
	}


	////T-B-N-S
	if(block.lock()->hasInterpolationFlagFC(TN)&& !block.lock()->hasInterpolationFlagFC(N) && !block.lock()->hasInterpolationFlagFC(T))
	{
		if (maxX2==TmaxX2) maxX2 -= 3; 
		if (maxX3==TmaxX3) maxX3 -= 3;
	}
	if(block.lock()->hasInterpolationFlagFC(BS)&& !block.lock()->hasInterpolationFlagFC(S) && !block.lock()->hasInterpolationFlagFC(B))
	{
		if (minX2==TminX2) minX2 += 4;
		if (minX3==TminX3) minX3 += 4;
	}
	if(block.lock()->hasInterpolationFlagFC(BN)&& !block.lock()->hasInterpolationFlagFC(N) && !block.lock()->hasInterpolationFlagFC(B))
	{
		if (maxX2==TmaxX2) maxX2 -= 3; 
		if (minX3==TminX3) minX3 += 4;
	}
	if(block.lock()->hasInterpolationFlagFC(TS) && !block.lock()->hasInterpolationFlagFC(S) && !block.lock()->hasInterpolationFlagFC(T))
	{
		if (minX2==TminX2) minX2 += 4;
		if (maxX3==TmaxX3) maxX3 -= 3;
	}

   //if (block.lock()->hasInterpolationFlagFC(D3Q27System::TNE)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::TE)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::TN)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::NE)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::T)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::N) && !block.lock()->hasInterpolationFlagFC(D3Q27System::E))
   //if (!block.lock()->hasInterpolationFlagFC(D3Q27System::TE)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::T) && !block.lock()->hasInterpolationFlagFC(D3Q27System::E))
   //{
   //   if (maxX1==TmaxX1) maxX1 -= 3;
   //   if (maxX2==TmaxX2) maxX2 -= 3;
   //   if (maxX3==TmaxX3) maxX3 -= 3;
   //}
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  >
void D3Q27ETFCOffVectorConnector< VectorTransmitter>::getLocalMinMax(int& minX1, int& minX2, int& minX3, int& maxX1, int& maxX2, int& maxX3, CFconnectorType connType)
{
   using namespace D3Q27System;
   int TminX1 = minX1; int TminX2 = minX2; int TminX3 = minX3; int TmaxX1 = maxX1; int TmaxX2 = maxX2; int TmaxX3 = maxX3;

   if (block.lock()->hasInterpolationFlagFC(E))
   {
      if (maxX1==TmaxX1) maxX1 -= 3;
   }
   if (block.lock()->hasInterpolationFlagFC(W))
   {
      if (minX1==TminX1) minX1 += 4;
   }
   if (block.lock()->hasInterpolationFlagFC(N))
   {
      if (maxX2==TmaxX2) maxX2 -= 3;
   }
   if (block.lock()->hasInterpolationFlagFC(S))
   {
      if (minX2==TminX2) minX2 += 4;
   }
   if (block.lock()->hasInterpolationFlagFC(T))
   {
      if (maxX3==TmaxX3) maxX3 -= 3;
   }
   if (block.lock()->hasInterpolationFlagFC(B))
   {
      if (minX3==TminX3) minX3 += 4;
   }

   ////////////
   /////E-W-N-S
   if (block.lock()->hasInterpolationFlagFC(NE)&& !block.lock()->hasInterpolationFlagFC(N) && !block.lock()->hasInterpolationFlagFC(E))
   {
      if (maxX1==TmaxX1) maxX1 -= 3;
      if (maxX2==TmaxX2) maxX2 -= 3;
   }
   if (block.lock()->hasInterpolationFlagFC(SW)&& !block.lock()->hasInterpolationFlagFC(W) && !block.lock()->hasInterpolationFlagFC(S))
   {
      if (minX1==TminX1) minX1 += 4;
      if (minX2==TminX2) minX2 += 4;
   }
   if (block.lock()->hasInterpolationFlagFC(SE)&& !block.lock()->hasInterpolationFlagFC(E) && !block.lock()->hasInterpolationFlagFC(S))
   {
      if (maxX1==TmaxX1) maxX1 -= 3;
      if (minX2==TminX2) minX2 += 4;
   }
   if (block.lock()->hasInterpolationFlagFC(NW)&& !block.lock()->hasInterpolationFlagFC(N) && !block.lock()->hasInterpolationFlagFC(W))
   {
      if (minX1==TminX1) minX1 += 4;
      if (maxX2==TmaxX2) maxX2 -= 3;
   }

   //////T-B-E-W
   if (block.lock()->hasInterpolationFlagFC(TE) && !block.lock()->hasInterpolationFlagFC(E) && !block.lock()->hasInterpolationFlagFC(T))
   {
      if (maxX1==TmaxX1) maxX1 -= 3;
      if (maxX3==TmaxX3) maxX3 -= 3;
   }
   if (block.lock()->hasInterpolationFlagFC(BW)&& !block.lock()->hasInterpolationFlagFC(W) && !block.lock()->hasInterpolationFlagFC(B))
   {
      if (minX1==TminX1) minX1 += 4;
      if (minX3==TminX3) minX3 += 4;
   }
   if (block.lock()->hasInterpolationFlagFC(BE)&& !block.lock()->hasInterpolationFlagFC(E) && !block.lock()->hasInterpolationFlagFC(B))
   {
      if (maxX1==TmaxX1) maxX1 -= 3;
      if (minX3==TminX3) minX3 += 4;
   }
   if (block.lock()->hasInterpolationFlagFC(TW)&& !block.lock()->hasInterpolationFlagFC(W) && !block.lock()->hasInterpolationFlagFC(T))
   {
      if (minX1==TminX1) minX1 += 4;
      if (maxX3==TmaxX3) maxX3 -= 3;
   }


   ////T-B-N-S
   if (block.lock()->hasInterpolationFlagFC(TN)&& !block.lock()->hasInterpolationFlagFC(N) && !block.lock()->hasInterpolationFlagFC(T))
   {
      if (maxX2==TmaxX2) maxX2 -= 3;
      if (maxX3==TmaxX3) maxX3 -= 3;
   }
   if (block.lock()->hasInterpolationFlagFC(BS)&& !block.lock()->hasInterpolationFlagFC(S) && !block.lock()->hasInterpolationFlagFC(B))
   {
      if (minX2==TminX2) minX2 += 4;
      if (minX3==TminX3) minX3 += 4;
   }
   if (block.lock()->hasInterpolationFlagFC(BN)&& !block.lock()->hasInterpolationFlagFC(N) && !block.lock()->hasInterpolationFlagFC(B))
   {
      if (maxX2==TmaxX2) maxX2 -= 3;
      if (minX3==TminX3) minX3 += 4;
   }
   if (block.lock()->hasInterpolationFlagFC(TS) && !block.lock()->hasInterpolationFlagFC(S) && !block.lock()->hasInterpolationFlagFC(T))
   {
      if (minX2==TminX2) minX2 += 4;
      if (maxX3==TmaxX3) maxX3 -= 3;
   }

   //if (block.lock()->hasInterpolationFlagFC(D3Q27System::TNE)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::TE)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::TN)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::NE)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::T)&&!block.lock()->hasInterpolationFlagFC(D3Q27System::N) && !block.lock()->hasInterpolationFlagFC(D3Q27System::E))
   //{
   //   if (maxX1==TmaxX1) maxX1 -= 3;
   //   if (maxX2==TmaxX2) maxX2 -= 3;
   //   if (maxX3==TmaxX3) maxX3 -= 3;
   //}
}
//////////////////////////////////////////////////////////////////////////
template<  typename VectorTransmitter  > 
void D3Q27ETFCOffVectorConnector< VectorTransmitter>::getLocalOffsets(const int& gMax, int& oMin)
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
void D3Q27ETFCOffVectorConnector< VectorTransmitter>::getLocalMins(int& minX1, int& minX2, int& minX3, const int& oMinX1, const int& oMinX2, const int& oMinX3)
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

		/////
	case NE: case SW: case SE: case NW:
		//case SW:
		if(connType == OddEvenSE)
			//minX2 += oMinX2;
		if(connType == OddOddNE)
		{
			//minX2 += oMinX2;
			minX3 += oMinX3;
		}
		if(connType == EvenOddNW)
			minX3 += oMinX3;
		break;

		//////
	case TE: case BW: case BE: case TW:
		if(connType == OddEvenSE)
	//		minX1 += oMinX1;
		if(connType == OddOddNE)
		{
	//		minX1 += oMinX1;
			minX2 += oMinX2;
		}
		if(connType == EvenOddNW)
			minX2 += oMinX2;
		break;

	//	//////
	case TN: case BS: case BN: case TS:
		if(connType == OddEvenSE)
			minX1 += oMinX1;
		if(connType == OddOddNE)
		{
			minX1 += oMinX1;
			//minX3 += oMinX3;
		}
		if(connType == EvenOddNW)
			//minX3 += oMinX3;
		break;

	//	/////
	//	case TNE: case TNW: case TSE: case TSW: case BNE: case BNW: case BSE: case BSW:
	//	if(connType == OddEvenSE)
	//	//	minX1 += oMinX1;
	//	if(connType == OddOddNE)
	//	{
	//		//minX1 += oMinX1;
	//		//minX3 += oMinX3;
	//	}
	//	if(connType == EvenOddNW)
	//		//minX3 += oMinX3;
	//	break;
	}
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
double D3Q27ETFCOffVectorConnector<VectorTransmitter>::getSendRecieveTime()
{
	return 0;
}

#endif
