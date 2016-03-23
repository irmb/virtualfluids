/**
* @file D3Q27ETFullDirectConnector.h
* @brief Connector send and receive full distribution in shared memory
*         
* @author Kostyantyn Kucher
* @date 08.06.2011
*/
#ifndef D3Q27ETFULLDIRECTCONNECTOR2_H
#define D3Q27ETFULLDIRECTCONNECTOR2_H

#include <boost/weak_ptr.hpp>

#include "LocalBlock3DConnector.h"
#include "Block3D.h"
#include "D3Q27System.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"

class D3Q27ETFullDirectConnector2 : public LocalBlock3DConnector
{
public:
   //D3Q27ETFullDirectConnector2() {}
   D3Q27ETFullDirectConnector2(Block3DPtr from, Block3DPtr to, int sendDir);
   void init();
   void sendVectors();

 protected:
   //void fillData(EsoTwist3DPtr  fFrom, int x1, int x2, int x3);
   //void distributeData(EsoTwist3DPtr  fTo, int x1, int x2, int x3);

   //void fillData(int x1, int x2, int x3);
   //void distributeData(int x1, int x2, int x3);

   inline void fillData(int x1, int x2, int x3);
   inline void distributeData(int x1, int x2, int x3);
private:

   LBMReal f[D3Q27System::ENDF+1];

   int maxX1;
   int maxX2;
   int maxX3;

   CbArray4D <LBMReal,IndexerX4X3X2X1>::CbArray4DPtr localDistributionsFrom; 
   CbArray4D <LBMReal,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsFrom; 
   CbArray3D <LBMReal,IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsFrom;

   CbArray4D <LBMReal,IndexerX4X3X2X1>::CbArray4DPtr localDistributionsTo; 
   CbArray4D <LBMReal,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsTo; 
   CbArray3D <LBMReal,IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsTo;

   EsoTwist3DPtr  fFrom;
   EsoTwist3DPtr  fTo;  
};

//////////////////////////////////////////////////////////////////////////
inline void D3Q27ETFullDirectConnector2::fillData(int x1, int x2, int x3)
{
   if(invStep)
   {
      f[D3Q27System::E] = (*this->localDistributionsFrom)(D3Q27System::ET_E, x1,x2,x3);
      f[D3Q27System::N] = (*this->localDistributionsFrom)(D3Q27System::ET_N,x1,x2,x3);  
      f[D3Q27System::T] = (*this->localDistributionsFrom)(D3Q27System::ET_T,x1,x2,x3);
      f[D3Q27System::NE] = (*this->localDistributionsFrom)(D3Q27System::ET_NE,x1,x2,x3);
      f[D3Q27System::NW] = (*this->localDistributionsFrom)(D3Q27System::ET_NW,x1+1,x2,x3);
      f[D3Q27System::TE] = (*this->localDistributionsFrom)(D3Q27System::ET_TE,x1,x2,x3);
      f[D3Q27System::TW] = (*this->localDistributionsFrom)(D3Q27System::ET_TW, x1+1,x2,x3);
      f[D3Q27System::TN] = (*this->localDistributionsFrom)(D3Q27System::ET_TN,x1,x2,x3);
      f[D3Q27System::TS] = (*this->localDistributionsFrom)(D3Q27System::ET_TS,x1,x2+1,x3);
      f[D3Q27System::TNE] = (*this->localDistributionsFrom)(D3Q27System::ET_TNE,x1,x2,x3);
      f[D3Q27System::TNW] = (*this->localDistributionsFrom)(D3Q27System::ET_TNW,x1+1,x2,x3);
      f[D3Q27System::TSE] = (*this->localDistributionsFrom)(D3Q27System::ET_TSE,x1,x2+1,x3);
      f[D3Q27System::TSW] = (*this->localDistributionsFrom)(D3Q27System::ET_TSW,x1+1,x2+1,x3);

      f[D3Q27System::W ] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_W,x1+1,x2,x3  );
      f[D3Q27System::S ] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_S,x1,x2+1,x3  );
      f[D3Q27System::B ] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_B,x1,x2,x3+1  );
      f[D3Q27System::SW] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_SW,x1+1,x2+1,x3 );
      f[D3Q27System::SE] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_SE,x1,x2+1,x3 );
      f[D3Q27System::BW] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BW,x1+1,x2,x3+1 );
      f[D3Q27System::BE] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BE,x1,x2,x3+1 );
      f[D3Q27System::BS] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BS,x1,x2+1,x3+1 );
      f[D3Q27System::BN] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BN,x1,x2,x3+1 );
      f[D3Q27System::BSW] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1);
      f[D3Q27System::BSE] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BSE,x1,x2+1,x3+1);
      f[D3Q27System::BNW] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BNW,x1+1,x2,x3+1);
      f[D3Q27System::BNE] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BNE,x1,x2,x3+1);

      f[D3Q27System::ZERO] = (*this->zeroDistributionsFrom)(x1,x2,x3);
   }
   else
   {
      f[D3Q27System::INV_E] = (*this->localDistributionsFrom)(D3Q27System::ET_E, x1,x2,x3);
      f[D3Q27System::INV_N] = (*this->localDistributionsFrom)(D3Q27System::ET_N,x1,x2,x3);  
      f[D3Q27System::INV_T] = (*this->localDistributionsFrom)(D3Q27System::ET_T,x1,x2,x3);
      f[D3Q27System::INV_NE] = (*this->localDistributionsFrom)(D3Q27System::ET_NE,x1,x2,x3);
      f[D3Q27System::INV_NW] = (*this->localDistributionsFrom)(D3Q27System::ET_NW,x1+1,x2,x3);
      f[D3Q27System::INV_TE] = (*this->localDistributionsFrom)(D3Q27System::ET_TE,x1,x2,x3);
      f[D3Q27System::INV_TW] = (*this->localDistributionsFrom)(D3Q27System::ET_TW, x1+1,x2,x3);
      f[D3Q27System::INV_TN] = (*this->localDistributionsFrom)(D3Q27System::ET_TN,x1,x2,x3);
      f[D3Q27System::INV_TS] = (*this->localDistributionsFrom)(D3Q27System::ET_TS,x1,x2+1,x3);
      f[D3Q27System::INV_TNE] = (*this->localDistributionsFrom)(D3Q27System::ET_TNE,x1,x2,x3);
      f[D3Q27System::INV_TNW] = (*this->localDistributionsFrom)(D3Q27System::ET_TNW,x1+1,x2,x3);
      f[D3Q27System::INV_TSE] = (*this->localDistributionsFrom)(D3Q27System::ET_TSE,x1,x2+1,x3);
      f[D3Q27System::INV_TSW] = (*this->localDistributionsFrom)(D3Q27System::ET_TSW,x1+1,x2+1,x3);

      f[D3Q27System::INV_W ] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_W,x1+1,x2,x3  );
      f[D3Q27System::INV_S ] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_S,x1,x2+1,x3  );
      f[D3Q27System::INV_B ] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_B,x1,x2,x3+1  );
      f[D3Q27System::INV_SW] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_SW,x1+1,x2+1,x3 );
      f[D3Q27System::INV_SE] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_SE,x1,x2+1,x3 );
      f[D3Q27System::INV_BW] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BW,x1+1,x2,x3+1 );
      f[D3Q27System::INV_BE] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BE,x1,x2,x3+1 );
      f[D3Q27System::INV_BS] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BS,x1,x2+1,x3+1 );
      f[D3Q27System::INV_BN] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BN,x1,x2,x3+1 );
      f[D3Q27System::INV_BSW] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1);
      f[D3Q27System::INV_BSE] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BSE,x1,x2+1,x3+1);
      f[D3Q27System::INV_BNW] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BNW,x1+1,x2,x3+1);
      f[D3Q27System::INV_BNE] = (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BNE,x1,x2,x3+1);

      f[D3Q27System::ZERO] = (*this->zeroDistributionsFrom)(x1,x2,x3);
   }
}
//////////////////////////////////////////////////////////////////////////
inline void D3Q27ETFullDirectConnector2::distributeData(int x1, int x2, int x3)
{
   if(invStep)
   {
      (*this->localDistributionsTo)(D3Q27System::ET_E,x1,  x2,  x3) = f[D3Q27System::E];
      (*this->localDistributionsTo)(D3Q27System::ET_N,x1,  x2,  x3) = f[D3Q27System::N];
      (*this->localDistributionsTo)(D3Q27System::ET_T,x1,  x2,  x3) = f[D3Q27System::T];
      (*this->localDistributionsTo)(D3Q27System::ET_NE,x1,  x2,  x3) = f[D3Q27System::NE];
      (*this->localDistributionsTo)(D3Q27System::ET_NW,x1+1,x2,  x3) = f[D3Q27System::NW];
      (*this->localDistributionsTo)(D3Q27System::ET_TE,x1,  x2,  x3) = f[D3Q27System::TE];
      (*this->localDistributionsTo)(D3Q27System::ET_TW,x1+1,x2,  x3) = f[D3Q27System::TW];
      (*this->localDistributionsTo)(D3Q27System::ET_TN,x1,  x2,  x3) = f[D3Q27System::TN];
      (*this->localDistributionsTo)(D3Q27System::ET_TS,x1,  x2+1,x3) = f[D3Q27System::TS];
      (*this->localDistributionsTo)(D3Q27System::ET_TNE,x1,  x2,  x3) = f[D3Q27System::TNE];
      (*this->localDistributionsTo)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f[D3Q27System::TNW];
      (*this->localDistributionsTo)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f[D3Q27System::TSE];
      (*this->localDistributionsTo)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f[D3Q27System::TSW];

      (*this->nonLocalDistributionsTo)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f[D3Q27System::W ];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f[D3Q27System::S ];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f[D3Q27System::B ];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f[D3Q27System::SW];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f[D3Q27System::SE];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f[D3Q27System::BW];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f[D3Q27System::BE];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f[D3Q27System::BS];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f[D3Q27System::BN];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f[D3Q27System::BSW];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f[D3Q27System::BSE];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f[D3Q27System::BNW];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = f[D3Q27System::BNE];

      (*this->zeroDistributionsTo)(x1,x2,x3) = f[D3Q27System::ZERO];
   }
   else
   {
      (*this->localDistributionsTo)(D3Q27System::ET_E,x1,  x2,  x3) = f[D3Q27System::INV_E];
      (*this->localDistributionsTo)(D3Q27System::ET_N,x1,  x2,  x3) = f[D3Q27System::INV_N];
      (*this->localDistributionsTo)(D3Q27System::ET_T,x1,  x2,  x3) = f[D3Q27System::INV_T];
      (*this->localDistributionsTo)(D3Q27System::ET_NE,x1,  x2,  x3) = f[D3Q27System::INV_NE];
      (*this->localDistributionsTo)(D3Q27System::ET_NW,x1+1,x2,  x3) = f[D3Q27System::INV_NW];
      (*this->localDistributionsTo)(D3Q27System::ET_TE,x1,  x2,  x3) = f[D3Q27System::INV_TE];
      (*this->localDistributionsTo)(D3Q27System::ET_TW,x1+1,x2,  x3) = f[D3Q27System::INV_TW];
      (*this->localDistributionsTo)(D3Q27System::ET_TN,x1,  x2,  x3) = f[D3Q27System::INV_TN];
      (*this->localDistributionsTo)(D3Q27System::ET_TS,x1,  x2+1,x3) = f[D3Q27System::INV_TS];
      (*this->localDistributionsTo)(D3Q27System::ET_TNE,x1,  x2,  x3) = f[D3Q27System::INV_TNE];
      (*this->localDistributionsTo)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f[D3Q27System::INV_TNW];
      (*this->localDistributionsTo)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f[D3Q27System::INV_TSE];
      (*this->localDistributionsTo)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f[D3Q27System::INV_TSW];

      (*this->nonLocalDistributionsTo)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f[D3Q27System::INV_W ];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f[D3Q27System::INV_S ];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f[D3Q27System::INV_B ];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f[D3Q27System::INV_SW];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f[D3Q27System::INV_SE];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f[D3Q27System::INV_BW];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f[D3Q27System::INV_BE];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f[D3Q27System::INV_BS];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f[D3Q27System::INV_BN];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f[D3Q27System::INV_BSW];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f[D3Q27System::INV_BSE];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f[D3Q27System::INV_BNW];
      (*this->nonLocalDistributionsTo)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = f[D3Q27System::INV_BNE];

      (*this->zeroDistributionsTo)(x1,x2,x3) = f[D3Q27System::ZERO];
   }
}
#endif 

