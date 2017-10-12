#include "EsoTwistD3Q27SparseData.h"
#include "EsoTwistD3Q27System.h"

EsoTwistD3Q27SparseData::SparseData EsoTwistD3Q27SparseData::localDistributions;
EsoTwistD3Q27SparseData::SparseData EsoTwistD3Q27SparseData::nonLocalDistributions;
EsoTwistD3Q27SparseData::SparseData EsoTwistD3Q27SparseData::zeroDistributions;

size_t EsoTwistD3Q27SparseData::nx1=0;
size_t EsoTwistD3Q27SparseData::nx2=0;
size_t EsoTwistD3Q27SparseData::nx3=0;
size_t EsoTwistD3Q27SparseData::nx4=0;
size_t EsoTwistD3Q27SparseData::nx5=0;
//////////////////////////////////////////////////////////////////////////
EsoTwistD3Q27SparseData::EsoTwistD3Q27SparseData()
{
}
//////////////////////////////////////////////////////////////////////////
EsoTwistD3Q27SparseData::EsoTwistD3Q27SparseData( size_t ibx[3], size_t nx[3], size_t level, double value )
{
   this->NX1 = nx[0];
   this->NX2 = nx[1];
   this->NX3 = nx[2];

   this->ox1 = ibx[0]*nx[0]; 
   this->ox2 = ibx[1]*nx[1];
   this->ox3 = ibx[2]*nx[2];
   this->level = level;

   ld = SparseDataPtr(&localDistributions);
   nld = SparseDataPtr(&nonLocalDistributions);
   zd = SparseDataPtr(&zeroDistributions);

   for(int x3 = 0; x3 < NX3; x3++)
   {
      for(int x2 = 0; x2 < NX2; x2++)
      {
         for(int x1 = 0; x1 < NX1; x1++)
         {
            for(int f = 0; f < 13; f++)
            {
               ld->insert(std::make_pair(index4D(f,x1,x2,x3),value));
               nld->insert(std::make_pair(index4D(f,x1,x2,x3),value));
            }
            zd->insert(std::make_pair(index3D(x1,x2,x3),value));
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
EsoTwistD3Q27SparseData::~EsoTwistD3Q27SparseData()
{
   ld = SparseDataPtr();
   nld = SparseDataPtr();
   zd = SparseDataPtr();

}
//////////////////////////////////////////////////////////////////////////
void EsoTwistD3Q27SparseData::swap()
{
   std::swap( this->ld, this->nld );
}
//////////////////////////////////////////////////////////////////////////
void EsoTwistD3Q27SparseData::getDistribution(LBMReal* const f, size_t x1, size_t x2, size_t x3)
{
   size_t x1plus = x1+1;
   size_t x2plus = x2+1;
   size_t x3plus = x3+1;

   f[D3Q27System::E] = (*ld)[index4D(D3Q27System::ET_E, x1,x2,x3)];
   f[D3Q27System::N] = (*ld)[index4D(D3Q27System::ET_N,x1,x2,x3)];  
   f[D3Q27System::T] = (*ld)[index4D(D3Q27System::ET_T,x1,x2,x3)];
   f[D3Q27System::NE] = (*ld)[index4D(D3Q27System::ET_NE,x1,x2,x3)];
   f[D3Q27System::NW] = (*ld)[index4D(D3Q27System::ET_NW,x1plus,x2,x3)];
   f[D3Q27System::TE] = (*ld)[index4D(D3Q27System::ET_TE,x1,x2,x3)];
   f[D3Q27System::TW] = (*ld)[index4D(D3Q27System::ET_TW, x1plus,x2,x3)];
   f[D3Q27System::TN] = (*ld)[index4D(D3Q27System::ET_TN,x1,x2,x3)];
   f[D3Q27System::TS] = (*ld)[index4D(D3Q27System::ET_TS,x1,x2plus,x3)];
   f[D3Q27System::TNE] = (*ld)[index4D(D3Q27System::ET_TNE,x1,x2,x3)];
   f[D3Q27System::TNW] = (*ld)[index4D(D3Q27System::ET_TNW,x1plus,x2,x3)];
   f[D3Q27System::TSE] = (*ld)[index4D(D3Q27System::ET_TSE,x1,x2plus,x3)];
   f[D3Q27System::TSW] = (*ld)[index4D(D3Q27System::ET_TSW,x1plus,x2plus,x3)];

   f[D3Q27System::W ] = (*nld)[index4D(D3Q27System::ET_W,x1plus,x2,x3  )];
   f[D3Q27System::S ] = (*nld)[index4D(D3Q27System::ET_S,x1,x2plus,x3  )];
   f[D3Q27System::B ] = (*nld)[index4D(D3Q27System::ET_B,x1,x2,x3plus  )];
   f[D3Q27System::SW] = (*nld)[index4D(D3Q27System::ET_SW,x1plus,x2plus,x3 )];
   f[D3Q27System::SE] = (*nld)[index4D(D3Q27System::ET_SE,x1,x2plus,x3 )];
   f[D3Q27System::BW] = (*nld)[index4D(D3Q27System::ET_BW,x1plus,x2,x3plus )];
   f[D3Q27System::BE] = (*nld)[index4D(D3Q27System::ET_BE,x1,x2,x3plus )];
   f[D3Q27System::BS] = (*nld)[index4D(D3Q27System::ET_BS,x1,x2plus,x3plus )];
   f[D3Q27System::BN] = (*nld)[index4D(D3Q27System::ET_BN,x1,x2,x3plus )];
   f[D3Q27System::BSW] = (*nld)[index4D(D3Q27System::ET_BSW,x1plus,x2plus,x3plus)];
   f[D3Q27System::BSE] = (*nld)[index4D(D3Q27System::ET_BSE,x1,x2plus,x3plus)];
   f[D3Q27System::BNW] = (*nld)[index4D(D3Q27System::ET_BNW,x1plus,x2,x3plus)];
   f[D3Q27System::BNE] = (*nld)[index4D(D3Q27System::ET_BNE,x1,x2,x3plus)];

   f[D3Q27System::ZERO] = (*zd)[index3D(x1,x2,x3)];
}
////////////////////////////////////////////////////////////////////////
void EsoTwistD3Q27SparseData::setDistribution(const LBMReal* const f, size_t x1, size_t x2, size_t x3)
{
   size_t x1plus = x1+1;
   size_t x2plus = x2+1;
   size_t x3plus = x3+1;

   (*ld)[index4D(D3Q27System::ET_E,x1,  x2,  x3)] = f[D3Q27System::INV_E];
   (*ld)[index4D(D3Q27System::ET_N,x1,  x2,  x3)] = f[D3Q27System::INV_N];
   (*ld)[index4D(D3Q27System::ET_T,x1,  x2,  x3)] = f[D3Q27System::INV_T];
   (*ld)[index4D(D3Q27System::ET_NE,x1,  x2,  x3)] = f[D3Q27System::INV_NE];
   (*ld)[index4D(D3Q27System::ET_NW,x1plus,x2,  x3)] = f[D3Q27System::INV_NW];
   (*ld)[index4D(D3Q27System::ET_TE,x1,  x2,  x3)] = f[D3Q27System::INV_TE];
   (*ld)[index4D(D3Q27System::ET_TW,x1plus,x2,  x3)] = f[D3Q27System::INV_TW];
   (*ld)[index4D(D3Q27System::ET_TN,x1,  x2,  x3)] = f[D3Q27System::INV_TN];
   (*ld)[index4D(D3Q27System::ET_TS,x1,  x2plus,x3)] = f[D3Q27System::INV_TS];
   (*ld)[index4D(D3Q27System::ET_TNE,x1,  x2,  x3)] = f[D3Q27System::INV_TNE];
   (*ld)[index4D(D3Q27System::ET_TNW,x1plus,x2,  x3)] = f[D3Q27System::INV_TNW];
   (*ld)[index4D(D3Q27System::ET_TSE,x1,  x2plus,x3)] = f[D3Q27System::INV_TSE];
   (*ld)[index4D(D3Q27System::ET_TSW,x1plus,x2plus,x3)] = f[D3Q27System::INV_TSW];

   (*nld)[index4D(D3Q27System::ET_W,x1plus,x2,  x3    )] = f[D3Q27System::INV_W ];
   (*nld)[index4D(D3Q27System::ET_S,x1,  x2plus,x3    )] = f[D3Q27System::INV_S ];
   (*nld)[index4D(D3Q27System::ET_B,x1,  x2,  x3plus  )] = f[D3Q27System::INV_B ];
   (*nld)[index4D(D3Q27System::ET_SW,x1plus,x2plus,x3   )] = f[D3Q27System::INV_SW];
   (*nld)[index4D(D3Q27System::ET_SE,x1,  x2plus,x3   )] = f[D3Q27System::INV_SE];
   (*nld)[index4D(D3Q27System::ET_BW,x1plus,x2,  x3plus )] = f[D3Q27System::INV_BW];
   (*nld)[index4D(D3Q27System::ET_BE,x1,  x2,  x3plus )] = f[D3Q27System::INV_BE];
   (*nld)[index4D(D3Q27System::ET_BS,x1,  x2plus,x3plus )] = f[D3Q27System::INV_BS];
   (*nld)[index4D(D3Q27System::ET_BN,x1,  x2,  x3plus )] = f[D3Q27System::INV_BN];
   (*nld)[index4D(D3Q27System::ET_BSW,x1plus,x2plus,x3plus)] = f[D3Q27System::INV_BSW];
   (*nld)[index4D(D3Q27System::ET_BSE,x1,  x2plus,x3plus)] = f[D3Q27System::INV_BSE];
   (*nld)[index4D(D3Q27System::ET_BNW,x1plus,x2,  x3plus)] = f[D3Q27System::INV_BNW];
   (*nld)[index4D(D3Q27System::ET_BNE,x1,  x2,  x3plus)] = f[D3Q27System::INV_BNE];

   (*zd)[index3D(x1,x2,x3)] = f[D3Q27System::ZERO];
}
////////////////////////////////////////////////////////////////////////
void EsoTwistD3Q27SparseData::getDistributionInv(LBMReal* const f, size_t x1, size_t x2, size_t x3)
{
   size_t x1plus = x1+1;
   size_t x2plus = x2+1;
   size_t x3plus = x3+1;

   f[D3Q27System::INV_E] = (*ld)[index4D(D3Q27System::ET_E, x1,x2,x3)];
   f[D3Q27System::INV_N] = (*ld)[index4D(D3Q27System::ET_N,x1,x2,x3)];  
   f[D3Q27System::INV_T] = (*ld)[index4D(D3Q27System::ET_T,x1,x2,x3)];
   f[D3Q27System::INV_NE] = (*ld)[index4D(D3Q27System::ET_NE,x1,x2,x3)];
   f[D3Q27System::INV_NW] = (*ld)[index4D(D3Q27System::ET_NW,x1plus,x2,x3)];
   f[D3Q27System::INV_TE] = (*ld)[index4D(D3Q27System::ET_TE,x1,x2,x3)];
   f[D3Q27System::INV_TW] = (*ld)[index4D(D3Q27System::ET_TW, x1plus,x2,x3)];
   f[D3Q27System::INV_TN] = (*ld)[index4D(D3Q27System::ET_TN,x1,x2,x3)];
   f[D3Q27System::INV_TS] = (*ld)[index4D(D3Q27System::ET_TS,x1,x2plus,x3)];
   f[D3Q27System::INV_TNE] = (*ld)[index4D(D3Q27System::ET_TNE,x1,x2,x3)];
   f[D3Q27System::INV_TNW] = (*ld)[index4D(D3Q27System::ET_TNW,x1plus,x2,x3)];
   f[D3Q27System::INV_TSE] = (*ld)[index4D(D3Q27System::ET_TSE,x1,x2plus,x3)];
   f[D3Q27System::INV_TSW] = (*ld)[index4D(D3Q27System::ET_TSW,x1plus,x2plus,x3)];

   f[D3Q27System::INV_W ] = (*nld)[index4D(D3Q27System::ET_W,x1plus,x2,x3  )];
   f[D3Q27System::INV_S ] = (*nld)[index4D(D3Q27System::ET_S,x1,x2plus,x3  )];
   f[D3Q27System::INV_B ] = (*nld)[index4D(D3Q27System::ET_B,x1,x2,x3plus  )];
   f[D3Q27System::INV_SW] = (*nld)[index4D(D3Q27System::ET_SW,x1plus,x2plus,x3 )];
   f[D3Q27System::INV_SE] = (*nld)[index4D(D3Q27System::ET_SE,x1,x2plus,x3 )];
   f[D3Q27System::INV_BW] = (*nld)[index4D(D3Q27System::ET_BW,x1plus,x2,x3plus )];
   f[D3Q27System::INV_BE] = (*nld)[index4D(D3Q27System::ET_BE,x1,x2,x3plus )];
   f[D3Q27System::INV_BS] = (*nld)[index4D(D3Q27System::ET_BS,x1,x2plus,x3plus )];
   f[D3Q27System::INV_BN] = (*nld)[index4D(D3Q27System::ET_BN,x1,x2,x3plus )];
   f[D3Q27System::INV_BSW] = (*nld)[index4D(D3Q27System::ET_BSW,x1plus,x2plus,x3plus)];
   f[D3Q27System::INV_BSE] = (*nld)[index4D(D3Q27System::ET_BSE,x1,x2plus,x3plus)];
   f[D3Q27System::INV_BNW] = (*nld)[index4D(D3Q27System::ET_BNW,x1plus,x2,x3plus)];
   f[D3Q27System::INV_BNE] = (*nld)[index4D(D3Q27System::ET_BNE,x1,x2,x3plus)];

   f[D3Q27System::ZERO] = (*zd)[index3D(x1,x2,x3)];
}
//////////////////////////////////////////////////////////////////////////
void EsoTwistD3Q27SparseData::setDistributionInv(const LBMReal* const f, size_t x1, size_t x2, size_t x3)
{
   size_t x1plus = x1+1;
   size_t x2plus = x2+1;
   size_t x3plus = x3+1;

   (*ld)[index4D(D3Q27System::ET_E,x1,  x2,  x3)] = f[D3Q27System::E];
   (*ld)[index4D(D3Q27System::ET_N,x1,  x2,  x3)] = f[D3Q27System::N];
   (*ld)[index4D(D3Q27System::ET_T,x1,  x2,  x3)] = f[D3Q27System::T];
   (*ld)[index4D(D3Q27System::ET_NE,x1,  x2,  x3)] = f[D3Q27System::NE];
   (*ld)[index4D(D3Q27System::ET_NW,x1plus,x2,  x3)] = f[D3Q27System::NW];
   (*ld)[index4D(D3Q27System::ET_TE,x1,  x2,  x3)] = f[D3Q27System::TE];
   (*ld)[index4D(D3Q27System::ET_TW,x1plus,x2,  x3)] = f[D3Q27System::TW];
   (*ld)[index4D(D3Q27System::ET_TN,x1,  x2,  x3)] = f[D3Q27System::TN];
   (*ld)[index4D(D3Q27System::ET_TS,x1,  x2plus,x3)] = f[D3Q27System::TS];
   (*ld)[index4D(D3Q27System::ET_TNE,x1,  x2,  x3)] = f[D3Q27System::TNE];
   (*ld)[index4D(D3Q27System::ET_TNW,x1plus,x2,  x3)] = f[D3Q27System::TNW];
   (*ld)[index4D(D3Q27System::ET_TSE,x1,  x2plus,x3)] = f[D3Q27System::TSE];
   (*ld)[index4D(D3Q27System::ET_TSW,x1plus,x2plus,x3)] = f[D3Q27System::TSW];

   (*nld)[index4D(D3Q27System::ET_W,x1plus,x2,  x3    )] = f[D3Q27System::W ];
   (*nld)[index4D(D3Q27System::ET_S,x1,  x2plus,x3    )] = f[D3Q27System::S ];
   (*nld)[index4D(D3Q27System::ET_B,x1,  x2,  x3plus  )] = f[D3Q27System::B ];
   (*nld)[index4D(D3Q27System::ET_SW,x1plus,x2plus,x3   )] = f[D3Q27System::SW];
   (*nld)[index4D(D3Q27System::ET_SE,x1,  x2plus,x3   )] = f[D3Q27System::SE];
   (*nld)[index4D(D3Q27System::ET_BW,x1plus,x2,  x3plus )] = f[D3Q27System::BW];
   (*nld)[index4D(D3Q27System::ET_BE,x1,  x2,  x3plus )] = f[D3Q27System::BE];
   (*nld)[index4D(D3Q27System::ET_BS,x1,  x2plus,x3plus )] = f[D3Q27System::BS];
   (*nld)[index4D(D3Q27System::ET_BN,x1,  x2,  x3plus )] = f[D3Q27System::BN];
   (*nld)[index4D(D3Q27System::ET_BSW,x1plus,x2plus,x3plus)] = f[D3Q27System::BSW];
   (*nld)[index4D(D3Q27System::ET_BSE,x1,  x2plus,x3plus)] = f[D3Q27System::BSE];
   (*nld)[index4D(D3Q27System::ET_BNW,x1plus,x2,  x3plus)] = f[D3Q27System::BNW];
   (*nld)[index4D(D3Q27System::ET_BNE,x1,  x2,  x3plus)] = f[D3Q27System::BNE];

   (*zd)[index3D(x1,x2,x3)] = f[D3Q27System::ZERO];
}
//////////////////////////////////////////////////////////////////////////
void EsoTwistD3Q27SparseData::setDistributionForDirection(const LBMReal* const f, size_t x1, size_t x2, size_t x3, unsigned long int direction)
{
   size_t x1plus = x1+1;
   size_t x2plus = x2+1;
   size_t x3plus = x3+1;

   bool directionFlag = false;
   if ((direction & EsoTwistD3Q27System::etE) == EsoTwistD3Q27System::etE)
      (*nld)[index4D(D3Q27System::ET_W,x1plus,x2,  x3    )] = f[D3Q27System::E]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etW) == EsoTwistD3Q27System::etW)
      (*ld)[index4D(D3Q27System::ET_E,x1,  x2,  x3)] = f[D3Q27System::W]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etS) == EsoTwistD3Q27System::etS)
      (*ld)[index4D(D3Q27System::ET_N,x1,  x2,  x3)] = f[D3Q27System::S]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etN) == EsoTwistD3Q27System::etN)
      (*nld)[index4D(D3Q27System::ET_S,x1,  x2plus,x3    )] = f[D3Q27System::N]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etB) == EsoTwistD3Q27System::etB)
      (*ld)[index4D(D3Q27System::ET_T,x1,  x2,  x3)] = f[D3Q27System::B]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etT) == EsoTwistD3Q27System::etT)
      (*nld)[index4D(D3Q27System::ET_B,x1,  x2,  x3plus  )] = f[D3Q27System::T]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etSW) == EsoTwistD3Q27System::etSW)
      (*ld)[index4D(D3Q27System::ET_NE,x1,  x2,  x3)] = f[D3Q27System::SW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etNE) == EsoTwistD3Q27System::etNE)
      (*nld)[index4D(D3Q27System::ET_SW,x1plus,x2plus,x3   )] = f[D3Q27System::NE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etNW) == EsoTwistD3Q27System::etNW)
      (*nld)[index4D(D3Q27System::ET_SE,x1,  x2plus,x3   )] = f[D3Q27System::NW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etSE) == EsoTwistD3Q27System::etSE)
      (*ld)[index4D(D3Q27System::ET_NW,x1plus,x2,  x3)] = f[D3Q27System::SE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBW) == EsoTwistD3Q27System::etBW)
      (*ld)[index4D(D3Q27System::ET_TE,x1,  x2,  x3)] = f[D3Q27System::BW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTE) == EsoTwistD3Q27System::etTE)
      (*nld)[index4D(D3Q27System::ET_BW,x1plus,x2,  x3plus )] = f[D3Q27System::TE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTW) == EsoTwistD3Q27System::etTW)
      (*nld)[index4D(D3Q27System::ET_BE,x1,  x2,  x3plus )] = f[D3Q27System::TW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBE) == EsoTwistD3Q27System::etBE)
      (*ld)[index4D(D3Q27System::ET_TW,x1plus,x2,  x3)] = f[D3Q27System::BE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBS) == EsoTwistD3Q27System::etBS)
      (*ld)[index4D(D3Q27System::ET_TN,x1,  x2,  x3)] = f[D3Q27System::BS]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTN) == EsoTwistD3Q27System::etTN)
      (*nld)[index4D(D3Q27System::ET_BS,x1,  x2plus,x3plus )] = f[D3Q27System::TN]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTS) == EsoTwistD3Q27System::etTS)
      (*nld)[index4D(D3Q27System::ET_BN,x1,  x2,  x3plus )] = f[D3Q27System::TS]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBN) == EsoTwistD3Q27System::etBN)
      (*ld)[index4D(D3Q27System::ET_TS,x1,  x2plus,x3)] = f[D3Q27System::BN]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBSW) == EsoTwistD3Q27System::etBSW)
      (*ld)[index4D(D3Q27System::ET_TNE,x1,  x2,  x3)] = f[D3Q27System::BSW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTNE) == EsoTwistD3Q27System::etTNE)
      (*nld)[index4D(D3Q27System::ET_BSW,x1plus,x2plus,x3plus)] = f[D3Q27System::TNE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBSE) == EsoTwistD3Q27System::etBSE)
      (*ld)[index4D(D3Q27System::ET_TNW,x1plus,x2,  x3)] = f[D3Q27System::BSE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTNW) == EsoTwistD3Q27System::etTNW)
      (*nld)[index4D(D3Q27System::ET_BSE,x1,  x2plus,x3plus)] = f[D3Q27System::TNW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBNW) == EsoTwistD3Q27System::etBNW)
      (*ld)[index4D(D3Q27System::ET_TSE,x1,  x2plus,x3)] = f[D3Q27System::BNW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTSE) == EsoTwistD3Q27System::etTSE)
      (*nld)[index4D(D3Q27System::ET_BNW,x1plus,x2,  x3plus)] = f[D3Q27System::TSE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBNE) == EsoTwistD3Q27System::etBNE)
      (*ld)[index4D(D3Q27System::ET_TSW,x1plus,x2plus,x3)] = f[D3Q27System::BNE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTSW) == EsoTwistD3Q27System::etTSW)
      (*nld)[index4D(D3Q27System::ET_BNE,x1,  x2,  x3plus)] = f[D3Q27System::TSW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::ZERO) == EsoTwistD3Q27System::ZERO)
      (*zd)[index3D(x1,x2,x3)] = f[D3Q27System::ZERO]; directionFlag=true;
#ifdef _DEBUG
   if(!directionFlag)UB_THROW( UbException(UB_EXARGS, "Direction didn't find") );
#endif //DEBUG
}
//////////////////////////////////////////////////////////////////////////
void EsoTwistD3Q27SparseData::setDistributionForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, int direction)
{
   switch (direction)
   {
   case D3Q27System::E :
      (*nld)[index4D(D3Q27System::ET_W,x1+1,x2,  x3    )] = f;
      break;
   case D3Q27System::W :
      (*ld)[index4D(D3Q27System::ET_E,x1,  x2,  x3)] = f;
      break;
   case D3Q27System::S :
      (*ld)[index4D(D3Q27System::ET_N,x1,  x2,  x3)] = f;
      break;
   case D3Q27System::N :
      (*nld)[index4D(D3Q27System::ET_S,x1,  x2+1,x3    )] = f;
      break;
   case D3Q27System::B :
      (*ld)[index4D(D3Q27System::ET_T,x1,  x2,  x3)] = f;
      break;
   case D3Q27System::T :
      (*nld)[index4D(D3Q27System::ET_B,x1,  x2,  x3+1  )] = f;
      break;
   case D3Q27System::SW :
      (*ld)[index4D(D3Q27System::ET_NE,x1,  x2,  x3)] = f;
      break;
   case D3Q27System::NE :
      (*nld)[index4D(D3Q27System::ET_SW,x1+1,x2+1,x3   )] = f;
      break;
   case D3Q27System::NW :
      (*nld)[index4D(D3Q27System::ET_SE,x1,  x2+1,x3   )] = f;
      break;
   case D3Q27System::SE :
      (*ld)[index4D(D3Q27System::ET_NW,x1+1,x2,  x3)] = f;
      break;
   case D3Q27System::BW :
      (*ld)[index4D(D3Q27System::ET_TE,x1,  x2,  x3)] = f;
      break;
   case D3Q27System::TE :
      (*nld)[index4D(D3Q27System::ET_BW,x1+1,x2,  x3+1 )] = f;
      break;
   case D3Q27System::TW :
      (*nld)[index4D(D3Q27System::ET_BE,x1,  x2,  x3+1 )] = f;
      break;
   case D3Q27System::BE :
      (*ld)[index4D(D3Q27System::ET_TW,x1+1,x2,  x3)] = f;
      break;
   case D3Q27System::BS :
      (*ld)[index4D(D3Q27System::ET_TN,x1,  x2,  x3)] = f;
      break;
   case D3Q27System::TN :
      (*nld)[index4D(D3Q27System::ET_BS,x1,  x2+1,x3+1 )] = f;
      break;
   case D3Q27System::TS :
      (*nld)[index4D(D3Q27System::ET_BN,x1,  x2,  x3+1 )] = f;
      break;
   case D3Q27System::BN :
      (*ld)[index4D(D3Q27System::ET_TS,x1,  x2+1,x3)] = f;
      break;
   case D3Q27System::BSW :
      (*ld)[index4D(D3Q27System::ET_TNE,x1,  x2,  x3)] = f;
      break;
   case D3Q27System::TNE :
      (*nld)[index4D(D3Q27System::ET_BSW,x1+1,x2+1,x3+1)] = f;
      break;
   case D3Q27System::BSE :
      (*ld)[index4D(D3Q27System::ET_TNW,x1+1,x2,  x3)] = f;
      break;
   case D3Q27System::TNW :
      (*nld)[index4D(D3Q27System::ET_BSE,x1,  x2+1,x3+1)] = f;
      break;
   case D3Q27System::BNW :
      (*ld)[index4D(D3Q27System::ET_TSE,x1,  x2+1,x3)] = f;
      break;
   case D3Q27System::TSE :
      (*nld)[index4D(D3Q27System::ET_BNW,x1+1,x2,  x3+1)] = f;
      break;
   case D3Q27System::BNE :
      (*ld)[index4D(D3Q27System::ET_TSW,x1+1,x2+1,x3)] = f;
      break;
   case D3Q27System::TSW :
      (*nld)[index4D(D3Q27System::ET_BNE,x1,  x2,  x3+1)] = f;
      break;
   case D3Q27System::ZERO :
      (*zd)[index3D(x1,x2,x3)] = f;
      break;
   default:
      UB_THROW( UbException(UB_EXARGS, "Direction didn't find"));     
   }
}
//////////////////////////////////////////////////////////////////////////
void EsoTwistD3Q27SparseData::setDistributionInvForDirection(const LBMReal* const f, size_t x1, size_t x2, size_t x3, unsigned long int direction)
{
   size_t x1plus = x1+1;
   size_t x2plus = x2+1;
   size_t x3plus = x3+1;

   bool directionFlag = false;
   if ((direction & EsoTwistD3Q27System::etE) == EsoTwistD3Q27System::etE)
      (*ld)[index4D(D3Q27System::ET_E,x1,  x2,  x3)] = f[D3Q27System::E]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etW) == EsoTwistD3Q27System::etW)
      (*nld)[index4D(D3Q27System::ET_W,x1plus,x2,  x3    )] = f[D3Q27System::W]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etS) == EsoTwistD3Q27System::etS)
      (*nld)[index4D(D3Q27System::ET_S,x1,  x2plus,x3    )] = f[D3Q27System::S]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etN) == EsoTwistD3Q27System::etN)
      (*ld)[index4D(D3Q27System::ET_N,x1,  x2,  x3)] = f[D3Q27System::N]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etB) == EsoTwistD3Q27System::etB)
      (*nld)[index4D(D3Q27System::ET_B,x1,  x2,  x3plus  )] = f[D3Q27System::B]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etT) == EsoTwistD3Q27System::etT)
      (*ld)[index4D(D3Q27System::ET_T,x1,  x2,  x3)] = f[D3Q27System::T]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etSW) == EsoTwistD3Q27System::etSW)
      (*nld)[index4D(D3Q27System::ET_SW,x1plus,x2plus,x3   )] = f[D3Q27System::SW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etNE) == EsoTwistD3Q27System::etNE)
      (*ld)[index4D(D3Q27System::ET_NE,x1,  x2,  x3)] = f[D3Q27System::NE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etNW) == EsoTwistD3Q27System::etNW)
      (*ld)[index4D(D3Q27System::ET_NW,x1plus,x2,  x3)] = f[D3Q27System::NW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etSE) == EsoTwistD3Q27System::etSE)
      (*nld)[index4D(D3Q27System::ET_SE,x1,  x2plus,x3   )] = f[D3Q27System::SE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBW) == EsoTwistD3Q27System::etBW)
      (*nld)[index4D(D3Q27System::ET_BW,x1plus,x2,  x3plus )] = f[D3Q27System::BW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTE) == EsoTwistD3Q27System::etTE)
      (*ld)[index4D(D3Q27System::ET_TE,x1,  x2,  x3)] = f[D3Q27System::TE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTW) == EsoTwistD3Q27System::etTW)
      (*ld)[index4D(D3Q27System::ET_TW,x1plus,x2,  x3)] = f[D3Q27System::TW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBE) == EsoTwistD3Q27System::etBE)
      (*nld)[index4D(D3Q27System::ET_BE,x1,  x2,  x3plus )] = f[D3Q27System::BE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBS) == EsoTwistD3Q27System::etBS)
      (*nld)[index4D(D3Q27System::ET_BS,x1,  x2plus,x3plus )] = f[D3Q27System::BS]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTN) == EsoTwistD3Q27System::etTN)
      (*ld)[index4D(D3Q27System::ET_TN,x1,  x2,  x3)] = f[D3Q27System::TN]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTS) == EsoTwistD3Q27System::etTS)
      (*ld)[index4D(D3Q27System::ET_TS,x1,  x2plus,x3)] = f[D3Q27System::TS]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBN) == EsoTwistD3Q27System::etBN)
      (*nld)[index4D(D3Q27System::ET_BN,x1,  x2,  x3plus )] = f[D3Q27System::BN]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBSW) == EsoTwistD3Q27System::etBSW)
      (*nld)[index4D(D3Q27System::ET_BSW,x1plus,x2plus,x3plus)] = f[D3Q27System::BSW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTNE) == EsoTwistD3Q27System::etTNE)
      (*ld)[index4D(D3Q27System::ET_TNE,x1,  x2,  x3)] = f[D3Q27System::TNE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBSE) == EsoTwistD3Q27System::etBSE)
      (*nld)[index4D(D3Q27System::ET_BSE,x1,  x2plus,x3plus)] = f[D3Q27System::BSE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTNW) == EsoTwistD3Q27System::etTNW)
      (*ld)[index4D(D3Q27System::ET_TNW,x1plus,x2,  x3)] = f[D3Q27System::TNW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBNW) == EsoTwistD3Q27System::etBNW)
      (*nld)[index4D(D3Q27System::ET_BNW,x1plus,x2,  x3plus)] = f[D3Q27System::BNW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTSE) == EsoTwistD3Q27System::etTSE)
      (*ld)[index4D(D3Q27System::ET_TSE,x1,  x2plus,x3)] = f[D3Q27System::TSE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etBNE) == EsoTwistD3Q27System::etBNE)
      (*nld)[index4D(D3Q27System::ET_BNE,x1,  x2,  x3plus)] = f[D3Q27System::BNE]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::etTSW) == EsoTwistD3Q27System::etTSW)
      (*ld)[index4D(D3Q27System::ET_TSW,x1plus,x2plus,x3)] = f[D3Q27System::TSW]; directionFlag=true;
   if ((direction & EsoTwistD3Q27System::ZERO) == EsoTwistD3Q27System::ZERO)
      (*zd)[index3D(x1,x2,x3)] = f[D3Q27System::ZERO]; directionFlag=true;
#ifdef _DEBUG
   if(!directionFlag)UB_THROW( UbException(UB_EXARGS, "Direction didn't find") );
#endif //DEBUG
}
//////////////////////////////////////////////////////////////////////////
void EsoTwistD3Q27SparseData::setDistributionInvForDirection(LBMReal f, size_t x1, size_t x2, size_t x3, unsigned long int direction)
{
   switch (direction)
   {
   case D3Q27System::E :
      (*ld)[index4D(D3Q27System::ET_E,x1,  x2,  x3)] = f;
      break;
   case D3Q27System::W :
      (*nld)[index4D(D3Q27System::ET_W,x1+1,x2,  x3    )] = f;
      break;
   case D3Q27System::S :
      (*nld)[index4D(D3Q27System::ET_S,x1,  x2+1,x3    )] = f;
      break;
   case D3Q27System::N :
      (*ld)[index4D(D3Q27System::ET_N,x1,  x2,  x3)] = f;
      break;
   case D3Q27System::B :
      (*nld)[index4D(D3Q27System::ET_B,x1,  x2,  x3+1  )] = f;
      break;
   case D3Q27System::T :
      (*ld)[index4D(D3Q27System::ET_T,x1,  x2,  x3)] = f;
      break;
   case D3Q27System::SW :
      (*nld)[index4D(D3Q27System::ET_SW,x1+1,x2+1,x3   )] = f;
      break;
   case D3Q27System::NE :
      (*ld)[index4D(D3Q27System::ET_NE,x1,  x2,  x3)] = f;
      break;
   case D3Q27System::NW :
      (*ld)[index4D(D3Q27System::ET_NW,x1+1,x2,  x3)] = f;
      break;
   case D3Q27System::SE :
      (*nld)[index4D(D3Q27System::ET_SE,x1,  x2+1,x3   )] = f;
      break;
   case D3Q27System::BW :
      (*nld)[index4D(D3Q27System::ET_BW,x1+1,x2,  x3+1 )] = f;
      break;
   case D3Q27System::TE :
      (*ld)[index4D(D3Q27System::ET_TE,x1,  x2,  x3)] = f;
      break;
   case D3Q27System::TW :
      (*ld)[index4D(D3Q27System::ET_TW,x1+1,x2,  x3)] = f;
      break;
   case D3Q27System::BE :
      (*nld)[index4D(D3Q27System::ET_BE,x1,  x2,  x3+1 )] = f;
      break;
   case D3Q27System::BS :
      (*nld)[index4D(D3Q27System::ET_BS,x1,  x2+1,x3+1 )] = f;
      break;
   case D3Q27System::TN :
      (*ld)[index4D(D3Q27System::ET_TN,x1,  x2,  x3)] = f;
      break;
   case D3Q27System::TS :
      (*ld)[index4D(D3Q27System::ET_TS,x1,  x2+1,x3)] = f;
      break;
   case D3Q27System::BN :
      (*nld)[index4D(D3Q27System::ET_BN,x1,  x2,  x3+1 )] = f;
      break;
   case D3Q27System::BSW :
      (*nld)[index4D(D3Q27System::ET_BSW,x1+1,x2+1,x3+1)] = f;
      break;
   case D3Q27System::TNE :
      (*ld)[index4D(D3Q27System::ET_TNE,x1,  x2,  x3)] = f;
      break;
   case D3Q27System::BSE :
      (*nld)[index4D(D3Q27System::ET_BSE,x1,  x2+1,x3+1)] = f;
      break;
   case D3Q27System::TNW :
      (*ld)[index4D(D3Q27System::ET_TNW,x1+1,x2,  x3)] = f;
      break;
   case D3Q27System::BNW :
      (*nld)[index4D(D3Q27System::ET_BNW,x1+1,x2,  x3+1)] = f;
      break;
   case D3Q27System::TSE :
      (*ld)[index4D(D3Q27System::ET_TSE,x1,  x2+1,x3)] = f;
      break;
   case D3Q27System::BNE :
      (*nld)[index4D(D3Q27System::ET_BNE,x1,  x2,  x3+1)] = f;
      break;
   case D3Q27System::TSW :
      (*ld)[index4D(D3Q27System::ET_TSW,x1+1,x2+1,x3)] = f;
      break;
   case D3Q27System::ZERO :
      (*zd)[index3D(x1,x2,x3)] = f;
      break;
   default:
      UB_THROW( UbException(UB_EXARGS, "Direction didn't find"));     
   }
}
//////////////////////////////////////////////////////////////////////////
LBMReal EsoTwistD3Q27SparseData::getDistributionInvForDirection(size_t x1, size_t x2, size_t x3, int direction)
{
   switch (direction)
   {
   case D3Q27System::E :
      return (*nld)[index4D(D3Q27System::ET_W,x1+1,x2,  x3    )];
   case D3Q27System::W :
      return (*ld)[index4D(D3Q27System::ET_E,x1,  x2,  x3)];
   case D3Q27System::S :
      return (*ld)[index4D(D3Q27System::ET_N,x1,  x2,  x3)];
   case D3Q27System::N :
      return (*nld)[index4D(D3Q27System::ET_S,x1,  x2+1,x3    )];
   case D3Q27System::B :
      return (*ld)[index4D(D3Q27System::ET_T,x1,  x2,  x3)];
   case D3Q27System::T :
      return (*nld)[index4D(D3Q27System::ET_B,x1,  x2,  x3+1  )];
   case D3Q27System::SW :
      return (*ld)[index4D(D3Q27System::ET_NE,x1,  x2,  x3)];
   case D3Q27System::NE :
      return (*nld)[index4D(D3Q27System::ET_SW,x1+1,x2+1,x3   )];
   case D3Q27System::NW :
      return (*nld)[index4D(D3Q27System::ET_SE,x1,  x2+1,x3   )];
   case D3Q27System::SE :
      return (*ld)[index4D(D3Q27System::ET_NW,x1+1,x2,  x3)];
   case D3Q27System::BW :
      return (*ld)[index4D(D3Q27System::ET_TE,x1,  x2,  x3)];
   case D3Q27System::TE :
      return (*nld)[index4D(D3Q27System::ET_BW,x1+1,x2,  x3+1 )];
   case D3Q27System::TW :
      return (*nld)[index4D(D3Q27System::ET_BE,x1,  x2,  x3+1 )];
   case D3Q27System::BE :
      return (*ld)[index4D(D3Q27System::ET_TW,x1+1,x2,  x3)];
   case D3Q27System::BS :
      return (*ld)[index4D(D3Q27System::ET_TN,x1,  x2,  x3)];
   case D3Q27System::TN :
      return (*nld)[index4D(D3Q27System::ET_BS,x1,  x2+1,x3+1 )];
   case D3Q27System::TS :
      return (*nld)[index4D(D3Q27System::ET_BN,x1,  x2,  x3+1 )];
   case D3Q27System::BN :
      return (*ld)[index4D(D3Q27System::ET_TS,x1,  x2+1,x3)];
   case D3Q27System::BSW :
      return (*ld)[index4D(D3Q27System::ET_TNE,x1,  x2,  x3)];
   case D3Q27System::TNE :
      return (*nld)[index4D(D3Q27System::ET_BSW,x1+1,x2+1,x3+1)];
   case D3Q27System::BSE :
      return (*ld)[index4D(D3Q27System::ET_TNW,x1+1,x2,  x3)];
   case D3Q27System::TNW :
      return (*nld)[index4D(D3Q27System::ET_BSE,x1,  x2+1,x3+1)];
   case D3Q27System::BNW :
      return (*ld)[index4D(D3Q27System::ET_TSE,x1,  x2+1,x3)];
   case D3Q27System::TSE :
      return (*nld)[index4D(D3Q27System::ET_BNW,x1+1,x2,  x3+1)];
   case D3Q27System::BNE :
      return (*ld)[index4D(D3Q27System::ET_TSW,x1+1,x2+1,x3)];
   case D3Q27System::TSW :
      return (*nld)[index4D(D3Q27System::ET_BNE,x1,  x2,  x3+1)];
   case D3Q27System::ZERO :
      return (*zd)[index3D(x1,x2,x3)];
   default:
      UB_THROW( UbException(UB_EXARGS, "Direction didn't find") );     
   }
}
//////////////////////////////////////////////////////////////////////////
size_t EsoTwistD3Q27SparseData::getNX1() const
{
   return NX1;
}
//////////////////////////////////////////////////////////////////////////
size_t EsoTwistD3Q27SparseData::getNX2() const
{
   return NX2;
}
//////////////////////////////////////////////////////////////////////////
size_t EsoTwistD3Q27SparseData::getNX3() const
{
   return NX3;
}
//////////////////////////////////////////////////////////////////////////
EsoTwistD3Q27SparseData::SparseDataPtr EsoTwistD3Q27SparseData::getLocalDistributions()
{
   return ld;
}
//////////////////////////////////////////////////////////////////////////
EsoTwistD3Q27SparseData::SparseDataPtr EsoTwistD3Q27SparseData::getNonLocalDistributions()
{
   return nld;
}
//////////////////////////////////////////////////////////////////////////
EsoTwistD3Q27SparseData::SparseDataPtr EsoTwistD3Q27SparseData::getZeroDistributions()
{
   return zd;
}
//////////////////////////////////////////////////////////////////////////
void EsoTwistD3Q27SparseData::setSize( int nx[4] )
{
   nx1 = 13;
   nx2 = nx[0];
   nx3 = nx[1];
   nx4 = nx[2];
   nx5 = nx[3];
}
//////////////////////////////////////////////////////////////////////////




