/**
* @file D3Q27ETFullDirectConnector.h
* @brief Connector send and receive full distribution in shared memory
*
* @author Konstantin Kutscher
* @date 28.04.2016
*/
#ifndef D3Q27ETFULLDIRECTCONNECTOR_H
#define D3Q27ETFULLDIRECTCONNECTOR_H

#include "LocalBlock3DConnector.h"
#include "Block3D.h"
#include "D3Q27System.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"

//! \brief   Exchange data between blocks. 
//! \details Connector send and receive full distributions between two blocks in shared memory.
//! \author  Konstantin Kutscher

class D3Q27ETFullDirectConnector : public LocalBlock3DConnector
{
public:
   D3Q27ETFullDirectConnector(Block3DPtr from, Block3DPtr to, int sendDir);
   void init();
   void sendVectors();

protected:
   inline void exchangeData(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To);
private:
   int maxX1;
   int maxX2;
   int maxX3;

   CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsFrom;
   CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsFrom;
   CbArray3D <LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsFrom;

   CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsTo;
   CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsTo;
   CbArray3D <LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsTo;

   EsoTwist3DPtr  fFrom;
   EsoTwist3DPtr  fTo;
};


//////////////////////////////////////////////////////////////////////////
inline void D3Q27ETFullDirectConnector::exchangeData(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To)
{
   (*this->localDistributionsTo)(D3Q27System::ET_E, x1To, x2To, x3To) = (*this->localDistributionsFrom)(D3Q27System::ET_E, x1From, x2From, x3From);
   (*this->localDistributionsTo)(D3Q27System::ET_N, x1To, x2To, x3To) = (*this->localDistributionsFrom)(D3Q27System::ET_N, x1From, x2From, x3From);
   (*this->localDistributionsTo)(D3Q27System::ET_T, x1To, x2To, x3To) = (*this->localDistributionsFrom)(D3Q27System::ET_T, x1From, x2From, x3From);
   (*this->localDistributionsTo)(D3Q27System::ET_NE, x1To, x2To, x3To) = (*this->localDistributionsFrom)(D3Q27System::ET_NE, x1From, x2From, x3From);
   (*this->localDistributionsTo)(D3Q27System::ET_NW, x1To + 1, x2To, x3To) = (*this->localDistributionsFrom)(D3Q27System::ET_NW, x1From + 1, x2From, x3From);
   (*this->localDistributionsTo)(D3Q27System::ET_TE, x1To, x2To, x3To) = (*this->localDistributionsFrom)(D3Q27System::ET_TE, x1From, x2From, x3From);
   (*this->localDistributionsTo)(D3Q27System::ET_TW, x1To + 1, x2To, x3To) = (*this->localDistributionsFrom)(D3Q27System::ET_TW, x1From + 1, x2From, x3From);
   (*this->localDistributionsTo)(D3Q27System::ET_TN, x1To, x2To, x3To) = (*this->localDistributionsFrom)(D3Q27System::ET_TN, x1From, x2From, x3From);
   (*this->localDistributionsTo)(D3Q27System::ET_TS, x1To, x2To + 1, x3To) = (*this->localDistributionsFrom)(D3Q27System::ET_TS, x1From, x2From + 1, x3From);
   (*this->localDistributionsTo)(D3Q27System::ET_TNE, x1To, x2To, x3To) = (*this->localDistributionsFrom)(D3Q27System::ET_TNE, x1From, x2From, x3From);
   (*this->localDistributionsTo)(D3Q27System::ET_TNW, x1To + 1, x2To, x3To) = (*this->localDistributionsFrom)(D3Q27System::ET_TNW, x1From + 1, x2From, x3From);
   (*this->localDistributionsTo)(D3Q27System::ET_TSE, x1To, x2To + 1, x3To) = (*this->localDistributionsFrom)(D3Q27System::ET_TSE, x1From, x2From + 1, x3From);
   (*this->localDistributionsTo)(D3Q27System::ET_TSW, x1To + 1, x2To + 1, x3To) = (*this->localDistributionsFrom)(D3Q27System::ET_TSW, x1From + 1, x2From + 1, x3From);

   (*this->nonLocalDistributionsTo)(D3Q27System::ET_W, x1To + 1, x2To, x3To) = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_W, x1From + 1, x2From, x3From);
   (*this->nonLocalDistributionsTo)(D3Q27System::ET_S, x1To, x2To + 1, x3To) = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_S, x1From, x2From + 1, x3From);
   (*this->nonLocalDistributionsTo)(D3Q27System::ET_B, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_B, x1From, x2From, x3From + 1);
   (*this->nonLocalDistributionsTo)(D3Q27System::ET_SW, x1To + 1, x2To + 1, x3To) = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_SW, x1From + 1, x2From + 1, x3From);
   (*this->nonLocalDistributionsTo)(D3Q27System::ET_SE, x1To, x2To + 1, x3To) = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_SE, x1From, x2From + 1, x3From);
   (*this->nonLocalDistributionsTo)(D3Q27System::ET_BW, x1To + 1, x2To, x3To + 1) = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BW, x1From + 1, x2From, x3From + 1);
   (*this->nonLocalDistributionsTo)(D3Q27System::ET_BE, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BE, x1From, x2From, x3From + 1);
   (*this->nonLocalDistributionsTo)(D3Q27System::ET_BS, x1To, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BS, x1From, x2From + 1, x3From + 1);
   (*this->nonLocalDistributionsTo)(D3Q27System::ET_BN, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BN, x1From, x2From, x3From + 1);
   (*this->nonLocalDistributionsTo)(D3Q27System::ET_BSW, x1To + 1, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BSW, x1From + 1, x2From + 1, x3From + 1);
   (*this->nonLocalDistributionsTo)(D3Q27System::ET_BSE, x1To, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BSE, x1From, x2From + 1, x3From + 1);
   (*this->nonLocalDistributionsTo)(D3Q27System::ET_BNW, x1To + 1, x2To, x3To + 1) = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BNW, x1From + 1, x2From, x3From + 1);
   (*this->nonLocalDistributionsTo)(D3Q27System::ET_BNE, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BNE, x1From, x2From, x3From + 1);

   (*this->zeroDistributionsTo)(x1To, x2To, x3To) = (*this->zeroDistributionsFrom)(x1From, x2From, x3From);
}
#endif 

