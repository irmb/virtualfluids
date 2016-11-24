#include "D3Q27EsoTwist3DSplittedVector.h"
#include "EsoTwistD3Q27System.h"

D3Q27EsoTwist3DSplittedVector::D3Q27EsoTwist3DSplittedVector()
{
}
//////////////////////////////////////////////////////////////////////////
D3Q27EsoTwist3DSplittedVector::D3Q27EsoTwist3DSplittedVector( size_t nx1, size_t nx2, size_t nx3, LBMReal value )
{
   this->NX1 = nx1;
   this->NX2 = nx2;
   this->NX3 = nx3;

   this->localDistributions    = CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal,IndexerX4X3X2X1>(13, nx1+1, nx2+1, nx3+1, value));
   this->nonLocalDistributions = CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal,IndexerX4X3X2X1>(13, nx1+1, nx2+1, nx3+1, value));

   this->zeroDistributions     = CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal,IndexerX3X2X1>(nx1, nx2, nx3, value));
}
//////////////////////////////////////////////////////////////////////////
D3Q27EsoTwist3DSplittedVector::~D3Q27EsoTwist3DSplittedVector()
{

}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::swap()
{
   std::swap( this->localDistributions, this->nonLocalDistributions );
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::getDistribution(LBMReal* const f, size_t x1, size_t x2, size_t x3)
{
   f[D3Q27System::E] = (*this->localDistributions)(D3Q27System::ET_E, x1,x2,x3);
   f[D3Q27System::N] = (*this->localDistributions)(D3Q27System::ET_N,x1,x2,x3);  
   f[D3Q27System::T] = (*this->localDistributions)(D3Q27System::ET_T,x1,x2,x3);
   f[D3Q27System::NE] = (*this->localDistributions)(D3Q27System::ET_NE,x1,x2,x3);
   f[D3Q27System::NW] = (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,x3);
   f[D3Q27System::TE] = (*this->localDistributions)(D3Q27System::ET_TE,x1,x2,x3);
   f[D3Q27System::TW] = (*this->localDistributions)(D3Q27System::ET_TW, x1+1,x2,x3);
   f[D3Q27System::TN] = (*this->localDistributions)(D3Q27System::ET_TN,x1,x2,x3);
   f[D3Q27System::TS] = (*this->localDistributions)(D3Q27System::ET_TS,x1,x2+1,x3);
   f[D3Q27System::TNE] = (*this->localDistributions)(D3Q27System::ET_TNE,x1,x2,x3);
   f[D3Q27System::TNW] = (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,x3);
   f[D3Q27System::TSE] = (*this->localDistributions)(D3Q27System::ET_TSE,x1,x2+1,x3);
   f[D3Q27System::TSW] = (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3);

   f[D3Q27System::W ] = (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,x3  );
   f[D3Q27System::S ] = (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,x2+1,x3  );
   f[D3Q27System::B ] = (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,x2,x3+1  );
   f[D3Q27System::SW] = (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3 );
   f[D3Q27System::SE] = (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,x2+1,x3 );
   f[D3Q27System::BW] = (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,x3+1 );
   f[D3Q27System::BE] = (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,x2,x3+1 );
   f[D3Q27System::BS] = (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,x2+1,x3+1 );
   f[D3Q27System::BN] = (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,x2,x3+1 );
   f[D3Q27System::BSW] = (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1);
   f[D3Q27System::BSE] = (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,x2+1,x3+1);
   f[D3Q27System::BNW] = (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,x3+1);
   f[D3Q27System::BNE] = (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,x2,x3+1);

   f[D3Q27System::ZERO] = (*this->zeroDistributions)(x1,x2,x3);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setDistribution(const LBMReal* const f, size_t x1, size_t x2, size_t x3)
{
   (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3)   = f[D3Q27System::INV_E];
   (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3)   = f[D3Q27System::INV_N];
   (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3)   = f[D3Q27System::INV_T];
   (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3)  = f[D3Q27System::INV_NE];
   (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3)  = f[D3Q27System::INV_NW];
   (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3)  = f[D3Q27System::INV_TE];
   (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3)  = f[D3Q27System::INV_TW];
   (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3)  = f[D3Q27System::INV_TN];
   (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3)  = f[D3Q27System::INV_TS];
   (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f[D3Q27System::INV_TNE];
   (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f[D3Q27System::INV_TNW];
   (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f[D3Q27System::INV_TSE];
   (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f[D3Q27System::INV_TSW];

   (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f[D3Q27System::INV_W ];
   (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f[D3Q27System::INV_S ];
   (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f[D3Q27System::INV_B ];
   (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f[D3Q27System::INV_SW];
   (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f[D3Q27System::INV_SE];
   (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f[D3Q27System::INV_BW];
   (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f[D3Q27System::INV_BE];
   (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f[D3Q27System::INV_BS];
   (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f[D3Q27System::INV_BN];
   (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f[D3Q27System::INV_BSW];
   (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f[D3Q27System::INV_BSE];
   (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f[D3Q27System::INV_BNW];
   (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = f[D3Q27System::INV_BNE];

   (*this->zeroDistributions)(x1,x2,x3) = f[D3Q27System::ZERO];
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::getDistributionInv(LBMReal* const f, size_t x1, size_t x2, size_t x3)
{
   f[D3Q27System::INV_E] = (*this->localDistributions)(D3Q27System::ET_E, x1,x2,x3);
   f[D3Q27System::INV_N] = (*this->localDistributions)(D3Q27System::ET_N,x1,x2,x3);  
   f[D3Q27System::INV_T] = (*this->localDistributions)(D3Q27System::ET_T,x1,x2,x3);
   f[D3Q27System::INV_NE] = (*this->localDistributions)(D3Q27System::ET_NE,x1,x2,x3);
   f[D3Q27System::INV_NW] = (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,x3);
   f[D3Q27System::INV_TE] = (*this->localDistributions)(D3Q27System::ET_TE,x1,x2,x3);
   f[D3Q27System::INV_TW] = (*this->localDistributions)(D3Q27System::ET_TW, x1+1,x2,x3);
   f[D3Q27System::INV_TN] = (*this->localDistributions)(D3Q27System::ET_TN,x1,x2,x3);
   f[D3Q27System::INV_TS] = (*this->localDistributions)(D3Q27System::ET_TS,x1,x2+1,x3);
   f[D3Q27System::INV_TNE] = (*this->localDistributions)(D3Q27System::ET_TNE,x1,x2,x3);
   f[D3Q27System::INV_TNW] = (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,x3);
   f[D3Q27System::INV_TSE] = (*this->localDistributions)(D3Q27System::ET_TSE,x1,x2+1,x3);
   f[D3Q27System::INV_TSW] = (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3);

   f[D3Q27System::INV_W ] = (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,x3  );
   f[D3Q27System::INV_S ] = (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,x2+1,x3  );
   f[D3Q27System::INV_B ] = (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,x2,x3+1  );
   f[D3Q27System::INV_SW] = (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3 );
   f[D3Q27System::INV_SE] = (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,x2+1,x3 );
   f[D3Q27System::INV_BW] = (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,x3+1 );
   f[D3Q27System::INV_BE] = (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,x2,x3+1 );
   f[D3Q27System::INV_BS] = (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,x2+1,x3+1 );
   f[D3Q27System::INV_BN] = (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,x2,x3+1 );
   f[D3Q27System::INV_BSW] = (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1);
   f[D3Q27System::INV_BSE] = (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,x2+1,x3+1);
   f[D3Q27System::INV_BNW] = (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,x3+1);
   f[D3Q27System::INV_BNE] = (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,x2,x3+1);

   f[D3Q27System::ZERO] = (*this->zeroDistributions)(x1,x2,x3);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setDistributionInv(const LBMReal* const f, size_t x1, size_t x2, size_t x3)
{
   (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f[D3Q27System::E];
   (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f[D3Q27System::N];
   (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f[D3Q27System::T];
   (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f[D3Q27System::NE];
   (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3) = f[D3Q27System::NW];
   (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f[D3Q27System::TE];
   (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3) = f[D3Q27System::TW];
   (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f[D3Q27System::TN];
   (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3) = f[D3Q27System::TS];
   (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f[D3Q27System::TNE];
   (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f[D3Q27System::TNW];
   (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f[D3Q27System::TSE];
   (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f[D3Q27System::TSW];

   (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f[D3Q27System::W ];
   (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f[D3Q27System::S ];
   (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f[D3Q27System::B ];
   (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f[D3Q27System::SW];
   (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f[D3Q27System::SE];
   (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f[D3Q27System::BW];
   (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f[D3Q27System::BE];
   (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f[D3Q27System::BS];
   (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f[D3Q27System::BN];
   (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f[D3Q27System::BSW];
   (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f[D3Q27System::BSE];
   (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f[D3Q27System::BNW];
   (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = f[D3Q27System::BNE];

   (*this->zeroDistributions)(x1,x2,x3) = f[D3Q27System::ZERO];
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setDistributionForDirection(const LBMReal* const f, size_t x1, size_t x2, size_t x3, unsigned long int direction)
{
   bool directionFlag = false;
   if ((direction & EsoTwistD3Q27System::etE) == EsoTwistD3Q27System::etE)
      (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f[D3Q27System::E]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etW) == EsoTwistD3Q27System::etW)
      (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f[D3Q27System::W]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etS) == EsoTwistD3Q27System::etS)
      (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f[D3Q27System::S]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etN) == EsoTwistD3Q27System::etN)
      (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f[D3Q27System::N]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etB) == EsoTwistD3Q27System::etB)
      (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f[D3Q27System::B]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etT) == EsoTwistD3Q27System::etT)
      (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f[D3Q27System::T]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etSW) == EsoTwistD3Q27System::etSW)
      (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f[D3Q27System::SW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etNE) == EsoTwistD3Q27System::etNE)
      (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f[D3Q27System::NE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etNW) == EsoTwistD3Q27System::etNW)
      (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f[D3Q27System::NW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etSE) == EsoTwistD3Q27System::etSE)
      (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3) = f[D3Q27System::SE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBW) == EsoTwistD3Q27System::etBW)
      (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f[D3Q27System::BW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTE) == EsoTwistD3Q27System::etTE)
      (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f[D3Q27System::TE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTW) == EsoTwistD3Q27System::etTW)
      (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f[D3Q27System::TW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBE) == EsoTwistD3Q27System::etBE)
      (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3) = f[D3Q27System::BE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBS) == EsoTwistD3Q27System::etBS)
      (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f[D3Q27System::BS]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTN) == EsoTwistD3Q27System::etTN)
      (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f[D3Q27System::TN]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTS) == EsoTwistD3Q27System::etTS)
      (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f[D3Q27System::TS]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBN) == EsoTwistD3Q27System::etBN)
      (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3) = f[D3Q27System::BN]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBSW) == EsoTwistD3Q27System::etBSW)
      (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f[D3Q27System::BSW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTNE) == EsoTwistD3Q27System::etTNE)
      (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f[D3Q27System::TNE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBSE) == EsoTwistD3Q27System::etBSE)
      (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f[D3Q27System::BSE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTNW) == EsoTwistD3Q27System::etTNW)
      (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f[D3Q27System::TNW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBNW) == EsoTwistD3Q27System::etBNW)
      (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f[D3Q27System::BNW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTSE) == EsoTwistD3Q27System::etTSE)
      (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f[D3Q27System::TSE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBNE) == EsoTwistD3Q27System::etBNE)
      (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f[D3Q27System::BNE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTSW) == EsoTwistD3Q27System::etTSW)
      (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = f[D3Q27System::TSW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::ZERO) == EsoTwistD3Q27System::ZERO)
      (*this->zeroDistributions)(x1,x2,x3) = f[D3Q27System::ZERO]; directionFlag=true;
#ifdef _DEBUG
   if(!directionFlag)UB_THROW( UbException(UB_EXARGS, "Direction didn't find") );
#endif //DEBUG
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setDistributionForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, int direction)
{
   switch (direction)
   {
   case D3Q27System::E :
      (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f;
      break;
   case D3Q27System::W :
      (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f;
      break;
   case D3Q27System::S :
      (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f;
      break;
   case D3Q27System::N :
      (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f;
      break;
   case D3Q27System::B :
      (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f;
      break;
   case D3Q27System::T :
      (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f;
      break;
   case D3Q27System::SW :
      (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f;
      break;
   case D3Q27System::NE :
      (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f;
      break;
   case D3Q27System::NW :
      (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f;
      break;
   case D3Q27System::SE :
      (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3) = f;
      break;
   case D3Q27System::BW :
      (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f;
      break;
   case D3Q27System::TE :
      (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f;
      break;
   case D3Q27System::TW :
      (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f;
      break;
   case D3Q27System::BE :
      (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3) = f;
      break;
   case D3Q27System::BS :
      (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f;
      break;
   case D3Q27System::TN :
      (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f;
      break;
   case D3Q27System::TS :
      (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f;
      break;
   case D3Q27System::BN :
      (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3) = f;
      break;
   case D3Q27System::BSW :
      (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f;
      break;
   case D3Q27System::TNE :
      (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f;
      break;
   case D3Q27System::BSE :
      (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f;
      break;
   case D3Q27System::TNW :
      (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f;
      break;
   case D3Q27System::BNW :
      (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f;
      break;
   case D3Q27System::TSE :
      (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f;
      break;
   case D3Q27System::BNE :
      (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f;
      break;
   case D3Q27System::TSW :
      (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = f;
      break;
   case D3Q27System::ZERO :
      (*this->zeroDistributions)(x1,x2,x3) = f;
      break;
   default:
      UB_THROW( UbException(UB_EXARGS, "Direction didn't find") );     
   }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setDistributionInvForDirection(const LBMReal* const f, size_t x1, size_t x2, size_t x3, unsigned long int direction)
{
   bool directionFlag = false;
   if ((direction & EsoTwistD3Q27System::etE) == EsoTwistD3Q27System::etE)
       (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f[D3Q27System::E]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etW) == EsoTwistD3Q27System::etW)
      (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f[D3Q27System::W]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etS) == EsoTwistD3Q27System::etS)
       (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f[D3Q27System::S]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etN) == EsoTwistD3Q27System::etN)
      (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f[D3Q27System::N]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etB) == EsoTwistD3Q27System::etB)
       (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f[D3Q27System::B]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etT) == EsoTwistD3Q27System::etT)
      (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f[D3Q27System::T]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etSW) == EsoTwistD3Q27System::etSW)
       (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f[D3Q27System::SW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etNE) == EsoTwistD3Q27System::etNE)
      (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f[D3Q27System::NE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etNW) == EsoTwistD3Q27System::etNW)
       (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3) = f[D3Q27System::NW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etSE) == EsoTwistD3Q27System::etSE)
      (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f[D3Q27System::SE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBW) == EsoTwistD3Q27System::etBW)
       (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f[D3Q27System::BW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTE) == EsoTwistD3Q27System::etTE)
      (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f[D3Q27System::TE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTW) == EsoTwistD3Q27System::etTW)
       (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3) = f[D3Q27System::TW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBE) == EsoTwistD3Q27System::etBE)
      (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f[D3Q27System::BE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBS) == EsoTwistD3Q27System::etBS)
       (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f[D3Q27System::BS]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTN) == EsoTwistD3Q27System::etTN)
      (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f[D3Q27System::TN]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTS) == EsoTwistD3Q27System::etTS)
       (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3) = f[D3Q27System::TS]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBN) == EsoTwistD3Q27System::etBN)
      (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f[D3Q27System::BN]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBSW) == EsoTwistD3Q27System::etBSW)
       (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f[D3Q27System::BSW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTNE) == EsoTwistD3Q27System::etTNE)
      (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f[D3Q27System::TNE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBSE) == EsoTwistD3Q27System::etBSE)
       (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f[D3Q27System::BSE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTNW) == EsoTwistD3Q27System::etTNW)
      (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f[D3Q27System::TNW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBNW) == EsoTwistD3Q27System::etBNW)
       (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f[D3Q27System::BNW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTSE) == EsoTwistD3Q27System::etTSE)
      (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f[D3Q27System::TSE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBNE) == EsoTwistD3Q27System::etBNE)
       (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1)= f[D3Q27System::BNE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTSW) == EsoTwistD3Q27System::etTSW)
      (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f[D3Q27System::TSW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::ZERO) == EsoTwistD3Q27System::ZERO)
      (*this->zeroDistributions)(x1,x2,x3) = f[D3Q27System::ZERO]; directionFlag=true;
#ifdef _DEBUG
   if(!directionFlag)UB_THROW( UbException(UB_EXARGS, "Direction didn't find") );
#endif //DEBUG
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setDistributionInvForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, unsigned long int direction)
{
   switch (direction)
   {
   case D3Q27System::E :
      (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f;
      break;
   case D3Q27System::W :
      (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f;
      break;
   case D3Q27System::S :
      (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f;
      break;
   case D3Q27System::N :
      (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f;
      break;
   case D3Q27System::B :
      (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f;
      break;
   case D3Q27System::T :
      (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f;
      break;
   case D3Q27System::SW :
      (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f;
      break;
   case D3Q27System::NE :
      (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f;
      break;
   case D3Q27System::NW :
      (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3) = f;
      break;
   case D3Q27System::SE :
      (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f;
      break;
   case D3Q27System::BW :
      (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f;
      break;
   case D3Q27System::TE :
      (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f;
      break;
   case D3Q27System::TW :
      (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3) = f;
      break;
   case D3Q27System::BE :
      (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f;
      break;
   case D3Q27System::BS :
      (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f;
      break;
   case D3Q27System::TN :
      (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f;
      break;
   case D3Q27System::TS :
      (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3) = f;
      break;
   case D3Q27System::BN :
      (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f;
      break;
   case D3Q27System::BSW :
      (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f;
      break;
   case D3Q27System::TNE :
      (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f;
      break;
   case D3Q27System::BSE :
      (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f;
      break;
   case D3Q27System::TNW :
      (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f;
      break;
   case D3Q27System::BNW :
      (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f;
      break;
   case D3Q27System::TSE :
      (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f;
      break;
   case D3Q27System::BNE :
      (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = f;
      break;
   case D3Q27System::TSW :
      (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f;
      break;
   case D3Q27System::ZERO :
      (*this->zeroDistributions)(x1,x2,x3) = f;
      break;
   default:
      UB_THROW( UbException(UB_EXARGS, "Direction didn't find") );     
   }
}
//////////////////////////////////////////////////////////////////////////
LBMReal D3Q27EsoTwist3DSplittedVector::getDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction)
{
   switch (direction)
   {
   case D3Q27System::W :
      return (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    );
   case D3Q27System::E :
      return (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3);
   case D3Q27System::N :
      return (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3);
   case D3Q27System::S :
      return (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    );
   case D3Q27System::T :
      return (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3);
   case D3Q27System::B :
      return (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  );
   case D3Q27System::NE :
      return (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3);
   case D3Q27System::SW :
      return (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   );
   case D3Q27System::SE :
      return (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   );
   case D3Q27System::NW :
      return (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3);
   case D3Q27System::TE :
      return (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3);
   case D3Q27System::BW :
      return (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 );
   case D3Q27System::BE :
      return (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 );
   case D3Q27System::TW :
      return (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3);
   case D3Q27System::TN :
      return (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3);
   case D3Q27System::BS :
      return (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 );
   case D3Q27System::BN :
      return (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 );
   case D3Q27System::TS :
      return (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3);
   case D3Q27System::TNE :
      return (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3);
   case D3Q27System::BSW :
      return (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1);
   case D3Q27System::TNW :
      return (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3);
   case D3Q27System::BSE :
      return (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1);
   case D3Q27System::TSE :
      return (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3);
   case D3Q27System::BNW :
      return (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1);
   case D3Q27System::TSW :
      return (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3);
   case D3Q27System::BNE :
      return (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1);
   case D3Q27System::ZERO :
      return (*this->zeroDistributions)(x1,x2,x3);
   default:
      UB_THROW( UbException(UB_EXARGS, "Direction didn't find") );     
   }
}
//////////////////////////////////////////////////////////////////////////
LBMReal D3Q27EsoTwist3DSplittedVector::getDistributionInvForDirection(size_t x1, size_t x2, size_t x3, int direction)
{
   switch (direction)
   {
   case D3Q27System::E :
      return (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    );
   case D3Q27System::W :
      return (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3);
   case D3Q27System::S :
      return (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3);
   case D3Q27System::N :
      return (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    );
   case D3Q27System::B :
      return (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3);
   case D3Q27System::T :
      return (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  );
   case D3Q27System::SW :
      return (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3);
   case D3Q27System::NE :
      return (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   );
   case D3Q27System::NW :
      return (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   );
   case D3Q27System::SE :
      return (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3);
   case D3Q27System::BW :
      return (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3);
   case D3Q27System::TE :
      return (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 );
   case D3Q27System::TW :
      return (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 );
   case D3Q27System::BE :
      return (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3);
   case D3Q27System::BS :
      return (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3);
   case D3Q27System::TN :
      return (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 );
   case D3Q27System::TS :
      return (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 );
   case D3Q27System::BN :
      return (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3);
   case D3Q27System::BSW :
      return (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3);
   case D3Q27System::TNE :
      return (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1);
   case D3Q27System::BSE :
      return (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3);
   case D3Q27System::TNW :
      return (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1);
   case D3Q27System::BNW :
      return (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3);
   case D3Q27System::TSE :
      return (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1);
   case D3Q27System::BNE :
      return (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3);
   case D3Q27System::TSW :
      return (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1);
   case D3Q27System::ZERO :
      return (*this->zeroDistributions)(x1,x2,x3);
   default:
      UB_THROW( UbException(UB_EXARGS, "Direction didn't find") );     
   }
}
//////////////////////////////////////////////////////////////////////////
size_t D3Q27EsoTwist3DSplittedVector::getNX1() const
{
   return NX1;
}
//////////////////////////////////////////////////////////////////////////
size_t D3Q27EsoTwist3DSplittedVector::getNX2() const
{
   return NX2;
}
//////////////////////////////////////////////////////////////////////////
size_t D3Q27EsoTwist3DSplittedVector::getNX3() const
{
   return NX3;
}
//////////////////////////////////////////////////////////////////////////
CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr D3Q27EsoTwist3DSplittedVector::getLocalDistributions()
{
   return this->localDistributions;
}
//////////////////////////////////////////////////////////////////////////
CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr D3Q27EsoTwist3DSplittedVector::getNonLocalDistributions()
{
   return this->nonLocalDistributions;
}
//////////////////////////////////////////////////////////////////////////
CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr D3Q27EsoTwist3DSplittedVector::getZeroDistributions()
{
   return this->zeroDistributions;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::getDistributionAfterLastStep( LBMReal* const f, size_t x1, size_t x2, size_t x3 )
{

}

//////////////////////////////////////////////////////////////////////////

