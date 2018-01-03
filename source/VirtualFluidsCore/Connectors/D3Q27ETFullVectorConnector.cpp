#include "D3Q27ETFullVectorConnector.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "LBMKernel.h"
#include "DataSet3D.h"
//////////////////////////////////////////////////////////////////////////
D3Q27ETFullVectorConnector::D3Q27ETFullVectorConnector(Block3DPtr block
   , VectorTransmitterPtr sender
   , VectorTransmitterPtr receiver
   , int sendDir)
   : RemoteBlock3DConnector(block, sender, receiver, sendDir)
{
   if (!block || !sender || !receiver)
      UB_THROW(UbException(UB_EXARGS, "sender or receiver == NULL!!"));

}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETFullVectorConnector::init()
{
   maxX1 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX1() - 1;
   maxX2 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX2() - 1;
   maxX3 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX3() - 1;

   fDis = std::dynamic_pointer_cast<EsoTwist3D>(block.lock()->getKernel()->getDataSet()->getFdistributions());

   int anz = 27;
   switch (sendDir)
   {
   case D3Q27System::ZERO: UB_THROW(UbException(UB_EXARGS, "ZERO not allowed")); break;
   case D3Q27System::E:
   case D3Q27System::W: sender->getData().resize(maxX2*maxX3*anz, 0.0);   break;
   case D3Q27System::N:
   case D3Q27System::S: sender->getData().resize(maxX1*maxX3*anz, 0.0);   break;
   case D3Q27System::T:
   case D3Q27System::B: sender->getData().resize(maxX1*maxX2*anz, 0.0);   break;

   case D3Q27System::NE:
   case D3Q27System::SW:
   case D3Q27System::SE:
   case D3Q27System::NW:  sender->getData().resize(maxX3*anz, 0.0);   break;

   case D3Q27System::TE:
   case D3Q27System::BW:
   case D3Q27System::BE:
   case D3Q27System::TW:  sender->getData().resize(maxX2*anz, 0.0);   break;

   case D3Q27System::TN:
   case D3Q27System::BS:
   case D3Q27System::BN:
   case D3Q27System::TS:  sender->getData().resize(maxX1*anz, 0.0);   break;

   case D3Q27System::TNE:
   case D3Q27System::BSW:
   case D3Q27System::BNE:
   case D3Q27System::TSW:
   case D3Q27System::TSE:
   case D3Q27System::BNW:
   case D3Q27System::BSE:
   case D3Q27System::TNW:  sender->getData().resize(anz, 0.0);   break;

   default: UB_THROW(UbException(UB_EXARGS, "unknown sendDir"));
   }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETFullVectorConnector::fillSendVectors()
{
   localDistributions = std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->fDis)->getLocalDistributions();
   nonLocalDistributions = std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->fDis)->getNonLocalDistributions();
   zeroDistributions = std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->fDis)->getZeroDistributions();

   vector_type& sdata = sender->getData();

   int index = 0;
   //EAST
   if (sendDir == D3Q27System::E)
   {
      for (int x3 = 1; x3 < maxX3; x3++)
      {
         for (int x2 = 1; x2 < maxX2; x2++)
         {
            fillData(sdata, index, maxX1 - 1, x2, x3);
         }
      }
   }
   //WEST
   else if (sendDir == D3Q27System::W)
   {
      for (int x3 = 1; x3 < maxX3; x3++)
      {
         for (int x2 = 1; x2 < maxX2; x2++)
         {
            fillData(sdata, index, 1, x2, x3);
         }
      }
   }
   //NORTH
   else if (sendDir == D3Q27System::N)
   {
      for (int x3 = 1; x3 < maxX3; x3++)
      {
         for (int x1 = 1; x1 < maxX1; x1++)
         {
            fillData(sdata, index, x1, maxX2 - 1, x3);
         }
      }
   }
   //SOUTH
   else if (sendDir == D3Q27System::S)
   {
      for (int x3 = 1; x3 < maxX3; x3++)
      {
         for (int x1 = 1; x1 < maxX1; x1++)
         {
            fillData(sdata, index, x1, 1, x3);
         }
      }
   }
   //TOP
   else if (sendDir == D3Q27System::T)
   {
      for (int x2 = 1; x2 < maxX2; x2++)
      {
         for (int x1 = 1; x1 < maxX1; x1++)
         {
            fillData(sdata, index, x1, x2, maxX3 - 1);
         }
      }
   }
   //BOTTOM
   else if (sendDir == D3Q27System::B)
   {
      for (int x2 = 1; x2 < maxX2; x2++)
      {
         for (int x1 = 1; x1 < maxX1; x1++)
         {
            fillData(sdata, index, x1, x2, 1);
         }
      }
   }
   //NE NW SW SE
   else if (sendDir == D3Q27System::NE || sendDir == D3Q27System::NW || sendDir == D3Q27System::SW || sendDir == D3Q27System::SE)
   {
      int x1 = 0;
      int x2 = 0;
      switch (sendDir)
      {
      case D3Q27System::NE:
         x1 = maxX1 - 1;
         x2 = maxX2 - 1;
         break;
      case D3Q27System::NW:
         x1 = 1;
         x2 = maxX2 - 1;
         break;
      case D3Q27System::SW:
         x1 = 1;
         x2 = 1;
         break;
      case D3Q27System::SE:
         x1 = maxX1 - 1;
         x2 = 1;
         break;
      }
      for (int x3 = 1; x3 < maxX3; x3++)
      {
         fillData(sdata, index, x1, x2, x3);
      }
   }
   //TE TW BW BE
   else if (sendDir == D3Q27System::TE || sendDir == D3Q27System::TW || sendDir == D3Q27System::BW || sendDir == D3Q27System::BE)
   {
      int x1 = 0;
      int x3 = 0;
      switch (sendDir)
      {
      case D3Q27System::TE:
         x1 = maxX1 - 1;
         x3 = maxX3 - 1;
         break;
      case D3Q27System::TW:
         x1 = 1;
         x3 = maxX3 - 1;
         break;
      case D3Q27System::BW:
         x1 = 1;
         x3 = 1;
         break;
      case D3Q27System::BE:
         x1 = maxX1 - 1;
         x3 = 1;
         break;
      }
      for (int x2 = 1; x2 < maxX2; x2++)
      {
         fillData(sdata, index, x1, x2, x3);
      }
   }
   //TN BN BS TS
   else if (sendDir == D3Q27System::TN || sendDir == D3Q27System::BN || sendDir == D3Q27System::BS || sendDir == D3Q27System::TS)
   {
      int x2 = 0;
      int x3 = 0;
      switch (sendDir)
      {
      case D3Q27System::TN:
         x3 = maxX3 - 1;
         x2 = maxX2 - 1;
         break;
      case D3Q27System::BN:
         x3 = 1;
         x2 = maxX2 - 1;
         break;
      case D3Q27System::BS:
         x3 = 1;
         x2 = 1;
         break;
      case D3Q27System::TS:
         x3 = maxX3 - 1;
         x2 = 1;
         break;
      }
      for (int x1 = 1; x1 < maxX1; x1++)
      {
         fillData(sdata, index, x1, x2, x3);
      }
   }
   //TNE TNW TSW TSE BNE BNW BSW BSE
   else if (sendDir == D3Q27System::TNE || sendDir == D3Q27System::TNW || sendDir == D3Q27System::TSW || sendDir == D3Q27System::TSE
      || sendDir == D3Q27System::BNE || sendDir == D3Q27System::BNW || sendDir == D3Q27System::BSW || sendDir == D3Q27System::BSE)
   {
      int x1 = 0;
      int x2 = 0;
      int x3 = 0;
      switch (sendDir)
      {
      case D3Q27System::TNE:   x1 = maxX1 - 1; x2 = maxX2 - 1; x3 = maxX3 - 1; break;
      case D3Q27System::TNW:   x1 = 1;       x2 = maxX2 - 1; x3 = maxX3 - 1; break;
      case D3Q27System::TSW:   x1 = 1;       x2 = 1;       x3 = maxX3 - 1; break;
      case D3Q27System::TSE:   x1 = maxX1 - 1; x2 = 1;       x3 = maxX3 - 1; break;
      case D3Q27System::BNE:   x1 = maxX1 - 1; x2 = maxX2 - 1; x3 = 1;       break;
      case D3Q27System::BNW:   x1 = 1;       x2 = maxX2 - 1; x3 = 1;       break;
      case D3Q27System::BSW:   x1 = 1;       x2 = 1;       x3 = 1;       break;
      case D3Q27System::BSE:   x1 = maxX1 - 1; x2 = 1;       x3 = 1;       break;
      }
      fillData(sdata, index, x1, x2, x3);
   }
   else UB_THROW(UbException(UB_EXARGS, "unknown dir"));
}
////////////////////////////////////////////////////////////////////////
void D3Q27ETFullVectorConnector::distributeReceiveVectors()
{
   /*e.g. connector sendet nach EAST --> empfaengt daten aus WEST ;-)*/

   localDistributions = std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->fDis)->getLocalDistributions();
   nonLocalDistributions = std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->fDis)->getNonLocalDistributions();
   zeroDistributions = std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->fDis)->getZeroDistributions();

   vector_type& rdata = receiver->getData();

   int index = 0;

   if (sendDir == D3Q27System::W)
   {
      for (int x3 = 1; x3 < maxX3; x3++)
      {
         for (int x2 = 1; x2 < maxX2; x2++)
         {
            distributeData(rdata, index, 0, x2, x3);
         }
      }
   }
   else if (sendDir == D3Q27System::E)
   {
      for (int x3 = 1; x3 < maxX3; x3++)
      {
         for (int x2 = 1; x2 < maxX2; x2++)
         {
            distributeData(rdata, index, maxX1, x2, x3);
         }
      }
   }
   else if (sendDir == D3Q27System::S)
   {
      for (int x3 = 1; x3 < maxX3; x3++)
      {
         for (int x1 = 1; x1 < maxX1; x1++)
         {
            distributeData(rdata, index, x1, 0, x3);
         }
      }
   }
   else if (sendDir == D3Q27System::N)
   {
      for (int x3 = 1; x3 < maxX3; x3++)
      {
         for (int x1 = 1; x1 < maxX1; x1++)
         {
            distributeData(rdata, index, x1, maxX2, x3);
         }
      }
   }
   else if (sendDir == D3Q27System::B)
   {
      for (int x2 = 1; x2 < maxX2; x2++)
      {
         for (int x1 = 1; x1 < maxX1; x1++)
         {
            distributeData(rdata, index, x1, x2, 0);
         }
      }
   }
   else if (sendDir == D3Q27System::T)
   {
      for (int x2 = 1; x2 < maxX2; x2++)
      {
         for (int x1 = 1; x1 < maxX1; x1++)
         {
            distributeData(rdata, index, x1, x2, maxX3);
         }
      }
   }
   //NE NW SW SE
   else if (sendDir == D3Q27System::NE || sendDir == D3Q27System::NW || sendDir == D3Q27System::SW || sendDir == D3Q27System::SE)
   {
      int x1 = 0;
      int x2 = 0;
      switch (sendDir)  //wenn sendir NE dann kommen werte von SW
      {
      case D3Q27System::NE:
         x1 = maxX1;
         x2 = maxX2;
         break;
      case D3Q27System::NW:
         x1 = 0;
         x2 = maxX2;
         break;
      case D3Q27System::SW:
         x1 = 0;
         x2 = 0;
         break;
      case D3Q27System::SE:
         x1 = maxX1;
         x2 = 0;
         break;
      }
      for (int x3 = 1; x3 < maxX3; x3++)
      {
         distributeData(rdata, index, x1, x2, x3);
      }

   }
   //TE TW BW BE
   else if (sendDir == D3Q27System::TE || sendDir == D3Q27System::TW || sendDir == D3Q27System::BW || sendDir == D3Q27System::BE)

   {
      int x1 = 0;
      int x3 = 0;
      switch (sendDir)  //wenn sendir NE dann kommen werte von SW
      {
      case D3Q27System::TE:
         x1 = maxX1;
         x3 = maxX3;
         break;
      case D3Q27System::TW:
         x1 = 0;
         x3 = maxX3;
         break;
      case D3Q27System::BW:
         x1 = 0;
         x3 = 0;
         break;
      case D3Q27System::BE:
         x1 = maxX1;
         x3 = 0;
         break;
      }
      for (int x2 = 1; x2 < maxX2; x2++)
      {
         distributeData(rdata, index, x1, x2, x3);
      }
   }
   //TN BN BS TS
   else if (sendDir == D3Q27System::TN || sendDir == D3Q27System::BN || sendDir == D3Q27System::BS || sendDir == D3Q27System::TS)
   {
      int x2 = 0;
      int x3 = 0;
      switch (sendDir)
      {
      case D3Q27System::TN:
         x3 = maxX3;
         x2 = maxX2;
         break;
      case D3Q27System::BN:
         x3 = 0;
         x2 = maxX2;
         break;
      case D3Q27System::BS:
         x3 = 0;
         x2 = 0;
         break;
      case D3Q27System::TS:
         x3 = maxX3;
         x2 = 0;
         break;

      }
      for (int x1 = 1; x1 < maxX1; x1++)
      {
         distributeData(rdata, index, x1, x2, x3);
      }
   }
   //TNE TNW TSW TSE BNE BNW BSW BSE
   else if (sendDir == D3Q27System::TNE || sendDir == D3Q27System::TNW || sendDir == D3Q27System::TSW || sendDir == D3Q27System::TSE
      || sendDir == D3Q27System::BNE || sendDir == D3Q27System::BNW || sendDir == D3Q27System::BSW || sendDir == D3Q27System::BSE)
   {
      int x1 = 0;
      int x2 = 0;
      int x3 = 0;

      switch (sendDir)
      {
      case D3Q27System::TNE:
         x1 = maxX1;
         x2 = maxX2;
         x3 = maxX3;
         break;
      case D3Q27System::TNW:
         x1 = 0;
         x2 = maxX2;
         x3 = maxX3;
         break;
      case D3Q27System::TSW:
         x1 = 0;
         x2 = 0;
         x3 = maxX3;
         break;
      case D3Q27System::TSE:
         x1 = maxX1;
         x2 = 0;
         x3 = maxX3;
         break;
      case D3Q27System::BNE:
         x1 = maxX1;
         x2 = maxX2;
         x3 = 0;
         break;
      case D3Q27System::BNW:
         x1 = 0;
         x2 = maxX2;
         x3 = 0;
         break;
      case D3Q27System::BSW:
         x1 = 0;
         x2 = 0;
         x3 = 0;
         break;
      case D3Q27System::BSE:
         x1 = maxX1;
         x2 = 0;
         x3 = 0;
         break;
      }
      distributeData(rdata, index, x1, x2, x3);
   }
   else UB_THROW(UbException(UB_EXARGS, "unknown dir"));
}
//////////////////////////////////////////////////////////////////////////


