#ifndef D3Q27ETFULLVECTORCONNECTOR_H
#define D3Q27ETFULLVECTORCONNECTOR_H

#include <vector>

#include "RemoteBlock3DConnector.h"
#include "D3Q27System.h"
#include "Block3D.h"
#include "LBMKernel.h"
#include "EsoTwistD3Q27System.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"
#include "EsoTwist3D.h"


//daten werden in einen vector (dieser befindet sich im transmitter) kopiert
//der vector wird via transmitter uebertragen
//transmitter kann ein lokal, MPI, RCG, CTL oder was auch immer fuer ein
//transmitter sein, der von Transmitter abgeleitet ist ;-)
class D3Q27ETFullVectorConnector : public RemoteBlock3DConnector
{
public:
   D3Q27ETFullVectorConnector(Block3DPtr block
      , VectorTransmitterPtr sender
      , VectorTransmitterPtr receiver
      , int sendDir);

   void init();

   void fillSendVectors();
   void distributeReceiveVectors();

protected:
   inline void fillData(vector_type& sdata, int& index, int x1, int x2, int x3);
   inline void distributeData(vector_type& rdata, int& index, int x1, int x2, int x3);
private:
   int maxX1;
   int maxX2;
   int maxX3;

   CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
   CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
   CbArray3D <LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributions;

   EsoTwist3DPtr  fDis;

};

//////////////////////////////////////////////////////////////////////////
inline void D3Q27ETFullVectorConnector::fillData(vector_type& sdata, int& index, int x1, int x2, int x3)
{
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_NW, x1 + 1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TW, x1 + 1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TS, x1, x2 + 1, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TNW, x1 + 1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2 + 1, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TSW, x1 + 1, x2 + 1, x3);

   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_W, x1 + 1, x2, x3);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1 + 1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1 + 1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1 + 1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1 + 1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3 + 1);

   sdata[index++] = (*this->zeroDistributions)(x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
inline void D3Q27ETFullVectorConnector::distributeData(vector_type& rdata, int& index, int x1, int x2, int x3)
{
   (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_NW, x1 + 1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_TW, x1 + 1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_TS, x1, x2 + 1, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_TNW, x1 + 1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2 + 1, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_TSW, x1 + 1, x2 + 1, x3) = rdata[index++];

   (*this->nonLocalDistributions)(D3Q27System::ET_W, x1 + 1, x2, x3) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2 + 1, x3) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1 + 1, x2 + 1, x3) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2 + 1, x3) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1 + 1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1 + 1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1 + 1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3 + 1) = rdata[index++];

   (*this->zeroDistributions)(x1, x2, x3) = rdata[index++];
}


#endif //D3Q27VECTORCONNECTOR_H

