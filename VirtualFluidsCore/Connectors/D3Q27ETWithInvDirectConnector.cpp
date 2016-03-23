#include "D3Q27ETWithInvDirectConnector.h"
#include "LBMKernelETD3Q27.h"
#include "EsoTwistD3Q27System.h"


using namespace std;

//Im Anschluss sind die Bulkbereiche synchron 

//*==========================================================*/
void D3Q27ETWithInvDirectConnector::sendVectors()
{
   EsoTwist3DPtr  fFrom =  boost::dynamic_pointer_cast<EsoTwist3D>(from.lock()->getKernel()->getDataSet()->getFdistributions());
   EsoTwist3DPtr  fTo   = boost::dynamic_pointer_cast<EsoTwist3D>(to.lock()->getKernel()->getDataSet()->getFdistributions());

   int maxX1 = (int)fFrom->getNX1()-1;
   int maxX2 = (int)fFrom->getNX2()-1;
   int maxX3 = (int)fFrom->getNX3()-1;

   LBMReal f[D3Q27System::ENDF+1];

   //EAST
   if(sendDir==D3Q27System::E)
   {
      for(int x3=1; x3<maxX3; x3++)   
      {
         for(int x2=1; x2<maxX2; x2++)   
         {
            if(invStep) fFrom->getDistribution(f,maxX1-1,x2,x3);
            else fFrom->getDistributionInv(f,maxX1-1,x2,x3);
            if(invStep) fTo->setDistributionInvForDirection(f,0,x2,x3,EsoTwistD3Q27System::etSE | EsoTwistD3Q27System::etE 
                                                      | EsoTwistD3Q27System::etNE | EsoTwistD3Q27System::etTE 
                                                      | EsoTwistD3Q27System::etBE | EsoTwistD3Q27System::etTNE
                                                      | EsoTwistD3Q27System::etTSE | EsoTwistD3Q27System::etBNE
                                                      | EsoTwistD3Q27System::etBSE);
            else fTo->setDistributionForDirection(f,0,x2,x3,EsoTwistD3Q27System::etSE | EsoTwistD3Q27System::etE 
                                                      | EsoTwistD3Q27System::etNE | EsoTwistD3Q27System::etTE 
                                                      | EsoTwistD3Q27System::etBE | EsoTwistD3Q27System::etTNE
                                                      | EsoTwistD3Q27System::etTSE | EsoTwistD3Q27System::etBNE
                                                      | EsoTwistD3Q27System::etBSE);
         }
      }
   }
   //WEST
   else if(sendDir==D3Q27System::W)
   {
      for(int x3=1; x3<maxX3; x3++)   
      {
         for(int x2=1; x2<maxX2; x2++)   
         {
            if(invStep) fFrom->getDistribution(f,1,x2,x3);
            else fFrom->getDistributionInv(f,1,x2,x3);
            if(invStep) fTo->setDistributionInvForDirection(f,maxX1,x2,x3,EsoTwistD3Q27System::etSW | EsoTwistD3Q27System::etW 
                                                      | EsoTwistD3Q27System::etNW | EsoTwistD3Q27System::etTW 
                                                      | EsoTwistD3Q27System::etBW | EsoTwistD3Q27System::etTNW
                                                      | EsoTwistD3Q27System::etTSW | EsoTwistD3Q27System::etBNW
                                                      | EsoTwistD3Q27System::etBSW);
            else fTo->setDistributionForDirection(f,maxX1,x2,x3,EsoTwistD3Q27System::etSW | EsoTwistD3Q27System::etW 
                                                        | EsoTwistD3Q27System::etNW | EsoTwistD3Q27System::etTW 
                                                        | EsoTwistD3Q27System::etBW | EsoTwistD3Q27System::etTNW
                                                        | EsoTwistD3Q27System::etTSW | EsoTwistD3Q27System::etBNW
                                                        | EsoTwistD3Q27System::etBSW);
         }
      }
   }
   //NORTH
   else if(sendDir==D3Q27System::N)
   {
      for(int x3=1; x3<maxX3; x3++)   
      {
         for(int x1=1; x1<maxX1; x1++)        
         {                                    
            if(invStep) fFrom->getDistribution(f,x1,maxX2-1,x3);
            else fFrom->getDistributionInv(f,x1,maxX2-1,x3);
            if(invStep) fTo->setDistributionInvForDirection(f,x1,0,x3,EsoTwistD3Q27System::etNW | EsoTwistD3Q27System::etN 
                                                      | EsoTwistD3Q27System::etNE | EsoTwistD3Q27System::etTN 
                                                      | EsoTwistD3Q27System::etBN | EsoTwistD3Q27System::etTNW
                                                      | EsoTwistD3Q27System::etTNE | EsoTwistD3Q27System::etBNW
                                                      | EsoTwistD3Q27System::etBNE);
            else fTo->setDistributionForDirection(f,x1,0,x3,EsoTwistD3Q27System::etNW | EsoTwistD3Q27System::etN 
                                                    | EsoTwistD3Q27System::etNE | EsoTwistD3Q27System::etTN 
                                                    | EsoTwistD3Q27System::etBN | EsoTwistD3Q27System::etTNW
                                                    | EsoTwistD3Q27System::etTNE | EsoTwistD3Q27System::etBNW
                                                    | EsoTwistD3Q27System::etBNE);
         }
      }
   }
   //SOUTH
   else if(sendDir==D3Q27System::S)
   {
      for(int x3=1; x3<maxX3; x3++)   
      {
         for(int x1=1; x1<maxX1; x1++)            
         {                                        
            if(invStep) fFrom->getDistribution(f,x1,1,x3);
            else fFrom->getDistributionInv(f,x1,1,x3);
            if(invStep) fTo->setDistributionInvForDirection(f,x1,maxX2,x3,EsoTwistD3Q27System::etSW | EsoTwistD3Q27System::etS 
                                                         | EsoTwistD3Q27System::etSE | EsoTwistD3Q27System::etTS
                                                         | EsoTwistD3Q27System::etBS | EsoTwistD3Q27System::etTSW
                                                         | EsoTwistD3Q27System::etTSE | EsoTwistD3Q27System::etBSW
                                                         | EsoTwistD3Q27System::etBSE);
            else fTo->setDistributionForDirection(f,x1,maxX2,x3,EsoTwistD3Q27System::etSW | EsoTwistD3Q27System::etS 
                                                        | EsoTwistD3Q27System::etSE | EsoTwistD3Q27System::etTS
                                                        | EsoTwistD3Q27System::etBS | EsoTwistD3Q27System::etTSW
                                                        | EsoTwistD3Q27System::etTSE | EsoTwistD3Q27System::etBSW
                                                        | EsoTwistD3Q27System::etBSE);
         }
      }
   }

   //TOP
   else if(sendDir==D3Q27System::T)
   {
      for(int x2=1; x2<maxX2; x2++)   
      {
         for(int x1=1; x1<maxX1; x1++)            
         {                                        
            if(invStep) fFrom->getDistribution(f,x1,x2,maxX3-1);
            else fFrom->getDistributionInv(f,x1,x2,maxX3-1);
            if(invStep) fTo->setDistributionInvForDirection(f,x1,x2,0,EsoTwistD3Q27System::etTW | EsoTwistD3Q27System::etT 
                                                      | EsoTwistD3Q27System::etTE | EsoTwistD3Q27System::etTS
                                                      | EsoTwistD3Q27System::etTN | EsoTwistD3Q27System::etTNE
                                                      | EsoTwistD3Q27System::etTNW | EsoTwistD3Q27System::etTSW
                                                      | EsoTwistD3Q27System::etTSE);
            else fTo->setDistributionForDirection(f,x1,x2,0,EsoTwistD3Q27System::etTW | EsoTwistD3Q27System::etT 
                                                      | EsoTwistD3Q27System::etTE | EsoTwistD3Q27System::etTS
                                                      | EsoTwistD3Q27System::etTN | EsoTwistD3Q27System::etTNE
                                                      | EsoTwistD3Q27System::etTNW | EsoTwistD3Q27System::etTSW
                                                      | EsoTwistD3Q27System::etTSE);
         }
      }
   }
   //BOTTOM
   else if(sendDir==D3Q27System::B)
   {
      for(int x2=1; x2<maxX2; x2++)   
      {
         for(int x1=1; x1<maxX1; x1++)            
         {                                        
            if(invStep) fFrom->getDistribution(f,x1,x2,1);
            else fFrom->getDistributionInv(f,x1,x2,1);
            if(invStep) fTo->setDistributionInvForDirection(f,x1,x2,maxX3,EsoTwistD3Q27System::etBW | EsoTwistD3Q27System::etB 
                                                         | EsoTwistD3Q27System::etBE | EsoTwistD3Q27System::etBS
                                                         | EsoTwistD3Q27System::etBN | EsoTwistD3Q27System::etBNE
                                                         | EsoTwistD3Q27System::etBNW | EsoTwistD3Q27System::etBSW
                                                         | EsoTwistD3Q27System::etBSE);
            else fTo->setDistributionForDirection(f,x1,x2,maxX3,EsoTwistD3Q27System::etBW | EsoTwistD3Q27System::etB 
                                                         | EsoTwistD3Q27System::etBE | EsoTwistD3Q27System::etBS
                                                         | EsoTwistD3Q27System::etBN | EsoTwistD3Q27System::etBNE
                                                         | EsoTwistD3Q27System::etBNW | EsoTwistD3Q27System::etBSW
                                                         | EsoTwistD3Q27System::etBSE);
         }
      }
   }
   //NORTHEAST
   else if(sendDir==D3Q27System::NE)
   {   
      for(int x3=1; x3<maxX3; x3++)   
      {
         if(invStep) fFrom->getDistribution(f,maxX1-1,maxX2-1,x3);
         else fFrom->getDistributionInv(f,maxX1-1,maxX2-1,x3);
         if(invStep) fTo->setDistributionInvForDirection(f,0,0,x3,EsoTwistD3Q27System::etNE
                                                   | EsoTwistD3Q27System::etTNE
                                                   | EsoTwistD3Q27System::etBNE);
         else fTo->setDistributionForDirection(f,0,0,x3,EsoTwistD3Q27System::etNE
                                                  | EsoTwistD3Q27System::etTNE
                                                  | EsoTwistD3Q27System::etBNE);
      }
   }
   //NORTHWEST
   else if(sendDir==D3Q27System::NW)
   {   
      for(int x3=1; x3<maxX3; x3++)   
      {
         if(invStep)fFrom->getDistribution(f,1,maxX2-1,x3);
         else fFrom->getDistributionInv(f,1,maxX2-1,x3);
         if(invStep) fTo->setDistributionInvForDirection(f,maxX1,0,x3,EsoTwistD3Q27System::etNW
                                                   | EsoTwistD3Q27System::etTNW
                                                   | EsoTwistD3Q27System::etBNW);
         else fTo->setDistributionForDirection(f,maxX1,0,x3,EsoTwistD3Q27System::etNW
                                                   | EsoTwistD3Q27System::etTNW
                                                   | EsoTwistD3Q27System::etBNW);
      }
   }
   //SOUTHWEST
   else if(sendDir==D3Q27System::SW)
   {  
      for(int x3=1; x3<maxX3; x3++)   
      {
         if(invStep) fFrom->getDistribution(f,1,1,x3);
         else fFrom->getDistributionInv(f,1,1,x3);
         if(invStep) fTo->setDistributionInvForDirection(f,maxX1,maxX2,x3,EsoTwistD3Q27System::etSW
                                                         | EsoTwistD3Q27System::etTSW
                                                         | EsoTwistD3Q27System::etBSW);
         else fTo->setDistributionForDirection(f,maxX1,maxX2,x3,EsoTwistD3Q27System::etSW
                                                         | EsoTwistD3Q27System::etTSW
                                                         | EsoTwistD3Q27System::etBSW);
      }
   }
   //SOUTHEAST
   else if(sendDir==D3Q27System::SE)
   {   
      for(int x3=1; x3<maxX3; x3++)   
      {
         if(invStep) fFrom->getDistribution(f,maxX1-1,1,x3);
         else fFrom->getDistributionInv(f,maxX1-1,1,x3);
         if(invStep) fTo->setDistributionInvForDirection(f,0,maxX2,x3,EsoTwistD3Q27System::etSE
                                                      | EsoTwistD3Q27System::etTSE
                                                      | EsoTwistD3Q27System::etBSE);
         else fTo->setDistributionForDirection(f,0,maxX2,x3,EsoTwistD3Q27System::etSE
                                                      | EsoTwistD3Q27System::etTSE
                                                      | EsoTwistD3Q27System::etBSE);
      }
   }
   else if(sendDir==D3Q27System::TE)
      for(int x2=1; x2<maxX2; x2++)
      {
         if(invStep) fFrom->getDistribution(f,maxX1-1,x2,maxX3-1);
         else fFrom->getDistributionInv(f,maxX1-1,x2,maxX3-1);
         if(invStep) fTo->setDistributionInvForDirection(f,0,x2,0,EsoTwistD3Q27System::etTE
                                                | EsoTwistD3Q27System::etTNE
                                                | EsoTwistD3Q27System::etTSE);
         else fTo->setDistributionForDirection(f,0,x2,0,EsoTwistD3Q27System::etTE
                                                | EsoTwistD3Q27System::etTNE
                                                | EsoTwistD3Q27System::etTSE);
      }
   else if(sendDir==D3Q27System::BW)
      for(int x2=1; x2<maxX2; x2++)
      {
         if(invStep) fFrom->getDistribution(f,1,x2,1);
         else fFrom->getDistributionInv(f,1,x2,1);
         if(invStep) fTo->setDistributionInvForDirection(f,maxX1,x2,maxX3,EsoTwistD3Q27System::etBW
                                                         | EsoTwistD3Q27System::etBNW
                                                         | EsoTwistD3Q27System::etBSW);
         else fTo->setDistributionForDirection(f,maxX1,x2,maxX3,EsoTwistD3Q27System::etBW
                                                         | EsoTwistD3Q27System::etBNW
                                                         | EsoTwistD3Q27System::etBSW);
      }
   else if(sendDir==D3Q27System::BE)
      for(int x2=1; x2<maxX2; x2++)
      {
         if(invStep) fFrom->getDistribution(f,maxX1-1,x2,1);
         else fFrom->getDistributionInv(f,maxX1-1,x2,1);
         if(invStep) fTo->setDistributionInvForDirection(f,0,x2,maxX3,EsoTwistD3Q27System::etBE
                                                      | EsoTwistD3Q27System::etBNE
                                                      | EsoTwistD3Q27System::etBSE);
         else fTo->setDistributionForDirection(f,0,x2,maxX3,EsoTwistD3Q27System::etBE
                                                      | EsoTwistD3Q27System::etBNE
                                                      | EsoTwistD3Q27System::etBSE);
      }
   else if(sendDir==D3Q27System::TW)
      for(int x2=1; x2<maxX2; x2++)
      {
         if(invStep) fFrom->getDistribution(f,1,x2,maxX3-1);
         else fFrom->getDistributionInv(f,1,x2,maxX3-1);
         if(invStep) fTo->setDistributionInvForDirection(f,maxX1,x2,0,EsoTwistD3Q27System::etTW
                                                      | EsoTwistD3Q27System::etTNW
                                                      | EsoTwistD3Q27System::etTSW);
         else fTo->setDistributionForDirection(f,maxX1,x2,0,EsoTwistD3Q27System::etTW
                                                      | EsoTwistD3Q27System::etTNW
                                                      | EsoTwistD3Q27System::etTSW);
      }
   else if(sendDir==D3Q27System::TN)
      for(int x1=1; x1<maxX1; x1++)
      {
         if(invStep) fFrom->getDistribution(f,x1,maxX2-1,maxX3-1);
         else fFrom->getDistributionInv(f,x1,maxX2-1,maxX3-1);
         if(invStep) fTo->setDistributionInvForDirection(f,x1,0,0,EsoTwistD3Q27System::etTN
                                                | EsoTwistD3Q27System::etTNW
                                                | EsoTwistD3Q27System::etTNE);
         else fTo->setDistributionForDirection(f,x1,0,0,EsoTwistD3Q27System::etTN
                                                | EsoTwistD3Q27System::etTNW
                                                | EsoTwistD3Q27System::etTNE);
      }
   else if(sendDir==D3Q27System::BS)
      for(int x1=1; x1<maxX1; x1++)
      {
         fFrom->getDistributionInv(f,x1,1,1);
         fTo->setDistributionForDirection(f,x1,maxX2,maxX3,EsoTwistD3Q27System::etBS
                                                         | EsoTwistD3Q27System::etBSW
                                                         | EsoTwistD3Q27System::etBSE);
      }
   else if(sendDir==D3Q27System::BN)
      for(int x1=1; x1<maxX1; x1++)
      {
         if(invStep) fFrom->getDistribution(f,x1,maxX2-1,1);
         else fFrom->getDistributionInv(f,x1,maxX2-1,1);
         if(invStep) fTo->setDistributionInvForDirection(f,x1,0,maxX3,EsoTwistD3Q27System::etBN
                                                      | EsoTwistD3Q27System::etBNW
                                                      | EsoTwistD3Q27System::etBNE);
         else fTo->setDistributionForDirection(f,x1,0,maxX3,EsoTwistD3Q27System::etBN
                                                      | EsoTwistD3Q27System::etBNW
                                                      | EsoTwistD3Q27System::etBNE);
      }

   else if(sendDir==D3Q27System::TS)
      for(int x1=1; x1<maxX1; x1++)
      {
         if(invStep) fFrom->getDistribution(f,x1,1,maxX3-1);
         else fFrom->getDistributionInv(f,x1,1,maxX3-1);
         if(invStep) fTo->setDistributionInvForDirection(f,x1,maxX2,0,EsoTwistD3Q27System::etTS
                                                      | EsoTwistD3Q27System::etTSW
                                                      | EsoTwistD3Q27System::etTSE);
         else fTo->setDistributionForDirection(f,x1,maxX2,0,EsoTwistD3Q27System::etTS
                                                      | EsoTwistD3Q27System::etTSW
                                                      | EsoTwistD3Q27System::etTSE);
      }

   else if(sendDir==D3Q27System::TSW)
   {
      if(invStep) fFrom->getDistribution(f,1,1,maxX3-1);
      else fFrom->getDistributionInv(f,1,1,maxX3-1);
      if(invStep) fTo->setDistributionInvForDirection(f,maxX1,maxX2,0,EsoTwistD3Q27System::etTSW);
      else fTo->setDistributionForDirection(f,maxX1,maxX2,0,EsoTwistD3Q27System::etTSW);
   }
   else if(sendDir==D3Q27System::TSE)
   {
      if(invStep) fFrom->getDistribution(f,maxX1-1,1,maxX3-1);
      else fFrom->getDistributionInv(f,maxX1-1,1,maxX3-1);
      if(invStep) fTo->setDistributionInvForDirection(f,0,maxX2,0,EsoTwistD3Q27System::etTSE);
      else fTo->setDistributionForDirection(f,0,maxX2,0,EsoTwistD3Q27System::etTSE);
   }
   else if(sendDir==D3Q27System::TNW)
   {
      if(invStep) fFrom->getDistribution(f,1,maxX2-1,maxX3-1);
      else fFrom->getDistributionInv(f,1,maxX2-1,maxX3-1);
      if(invStep) fTo->setDistributionInvForDirection(f,maxX1,0,0,EsoTwistD3Q27System::etTNW);
      else fTo->setDistributionForDirection(f,maxX1,0,0,EsoTwistD3Q27System::etTNW);
   }
   else if(sendDir==D3Q27System::TNE)
   {
      if(invStep) fFrom->getDistribution(f,maxX1-1,maxX2-1,maxX3-1);
      else fFrom->getDistributionInv(f,maxX1-1,maxX2-1,maxX3-1);
      if(invStep) fTo->setDistributionInvForDirection(f,0,0,0,EsoTwistD3Q27System::etTNE);
      else fTo->setDistributionForDirection(f,0,0,0,EsoTwistD3Q27System::etTNE);
   }
   else if(sendDir==D3Q27System::BSW)
   {
      if(invStep) fFrom->getDistribution(f,1,1,1);
      else fFrom->getDistributionInv(f,1,1,1);
      if(invStep) fTo->setDistributionInvForDirection(f,maxX1,maxX2,maxX3,EsoTwistD3Q27System::etBSW);
      else fTo->setDistributionForDirection(f,maxX1,maxX2,maxX3,EsoTwistD3Q27System::etBSW);
   }
   else if(sendDir==D3Q27System::BSE)
   {
      if(invStep) fFrom->getDistribution(f,maxX1-1,1,1);
      else fFrom->getDistributionInv(f,maxX1-1,1,1);
      if(invStep) fTo->setDistributionInvForDirection(f,0,maxX2,maxX3,EsoTwistD3Q27System::etBSE);
      else fTo->setDistributionForDirection(f,0,maxX2,maxX3,EsoTwistD3Q27System::etBSE);
   }
   else if(sendDir==D3Q27System::BNW)
   {
      if(invStep) fFrom->getDistribution(f,1,maxX2-1,1);
      else fFrom->getDistributionInv(f,1,maxX2-1,1);
      if(invStep) fTo->setDistributionInvForDirection(f,maxX1,0,maxX3,EsoTwistD3Q27System::etBNW);
      else fTo->setDistributionForDirection(f,maxX1,0,maxX3,EsoTwistD3Q27System::etBNW);
   }
   else if(sendDir==D3Q27System::BNE)
   {
      if(invStep) fFrom->getDistribution(f,maxX1-1,maxX2-1,1);
      else fFrom->getDistributionInv(f,maxX1-1,maxX2-1,1);
      if(invStep) fTo->setDistributionInvForDirection(f,0,0,maxX3,EsoTwistD3Q27System::etBNE);
      else fTo->setDistributionForDirection(f,0,0,maxX3,EsoTwistD3Q27System::etBNE);
   }
   else UB_THROW( UbException(UB_EXARGS,"unknown dir") );
}
//////////////////////////////////////////////////////////////////////////
double D3Q27ETWithInvDirectConnector::getSendRecieveTime()
{
   return 0;
}
