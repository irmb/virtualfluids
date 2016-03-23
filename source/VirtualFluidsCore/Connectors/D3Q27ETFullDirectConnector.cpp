#include "D3Q27ETFullDirectConnector.h"
#include "LBMKernelETD3Q27.h"
#include "EsoTwistD3Q27System.h"


using namespace std;

D3Q27ETFullDirectConnector::D3Q27ETFullDirectConnector(Block3DPtr from, Block3DPtr to, int sendDir)
   : LocalBlock3DConnector(from, to, sendDir)

{

}
//////////////////////////////////////////////////////////////////////////
D3Q27ETFullDirectConnector::~D3Q27ETFullDirectConnector()
{

}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETFullDirectConnector::init()
{
   fFrom =  boost::dynamic_pointer_cast<EsoTwist3D>(from.lock()->getKernel()->getDataSet()->getFdistributions());
   fTo   = boost::dynamic_pointer_cast<EsoTwist3D>(to.lock()->getKernel()->getDataSet()->getFdistributions());

   maxX2 = (int)fFrom->getNX2()-1;
   maxX3 = (int)fFrom->getNX3()-1;
   maxX1 = (int)fFrom->getNX1()-1;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETFullDirectConnector::sendVectors()
{
   //EAST
   if(sendDir==D3Q27System::E)
   {
      for(int x3=1; x3<maxX3; x3++)   
      {
         for(int x2=1; x2<maxX2; x2++)   
         {
            fillData(fFrom,maxX1-1,x2,x3);
            distributeData(fTo,0,x2,x3);
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
            fillData(fFrom,1,x2,x3);
            distributeData(fTo,maxX1,x2,x3);
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
            fillData(fFrom,x1,maxX2-1,x3);
            distributeData(fTo,x1,0,x3);
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
            fillData(fFrom,x1,1,x3);
            distributeData(fTo,x1,maxX2,x3);
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
            fillData(fFrom,x1,x2,maxX3-1);
            distributeData(fTo,x1,x2,0);
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
            fillData(fFrom,x1,x2,1);
            distributeData(fTo,x1,x2,maxX3);
         }
      }
   }
   //NORTHEAST
   else if(sendDir==D3Q27System::NE)
   {   
      for(int x3=1; x3<maxX3; x3++)   
      {
         fillData(fFrom,maxX1-1,maxX2-1,x3);
         distributeData(fTo,0,0,x3);
      }
   }
   //NORTHWEST
   else if(sendDir==D3Q27System::NW)
   {   
      for(int x3=1; x3<maxX3; x3++)   
      {
         fillData(fFrom,1,maxX2-1,x3);
         distributeData(fTo,maxX1,0,x3);
      }
   }
   //SOUTHWEST
   else if(sendDir==D3Q27System::SW)
   {  
      for(int x3=1; x3<maxX3; x3++)   
      {
         fillData(fFrom,1,1,x3);
         distributeData(fTo,maxX1,maxX2,x3);
      }
   }
   //SOUTHEAST
   else if(sendDir==D3Q27System::SE)
   {   
      for(int x3=1; x3<maxX3; x3++)   
      {
         fillData(fFrom,maxX1-1,1,x3);
         distributeData(fTo,0,maxX2,x3);
      }
   }
   else if(sendDir==D3Q27System::TE)
      for(int x2=1; x2<maxX2; x2++)
      {
         fillData(fFrom,maxX1-1,x2,maxX3-1);
         distributeData(fTo,0,x2,0);
      }
   else if(sendDir==D3Q27System::BW)
      for(int x2=1; x2<maxX2; x2++)
      {
         fillData(fFrom,1,x2,1);
         distributeData(fTo,maxX1,x2,maxX3);
      }
   else if(sendDir==D3Q27System::BE)
      for(int x2=1; x2<maxX2; x2++)
      {
         fillData(fFrom,maxX1-1,x2,1);
         distributeData(fTo,0,x2,maxX3);
      }
   else if(sendDir==D3Q27System::TW)
      for(int x2=1; x2<maxX2; x2++)
      {
         fillData(fFrom,1,x2,maxX3-1);
         distributeData(fTo,maxX1,x2,0);
      }
   else if(sendDir==D3Q27System::TN)
      for(int x1=1; x1<maxX1; x1++)
      {
         fillData(fFrom,x1,maxX2-1,maxX3-1);
         distributeData(fTo,x1,0,0);
      }
   else if(sendDir==D3Q27System::BS)
      for(int x1=1; x1<maxX1; x1++)
      {
         fillData(fFrom,x1,1,1);
         distributeData(fTo,x1,maxX2,maxX3);
      }
   else if(sendDir==D3Q27System::BN)
      for(int x1=1; x1<maxX1; x1++)
      {
         fillData(fFrom,x1,maxX2-1,1);
         distributeData(fTo,x1,0,maxX3);
      }

   else if(sendDir==D3Q27System::TS)
      for(int x1=1; x1<maxX1; x1++)
      {
         fillData(fFrom,x1,1,maxX3-1);
         distributeData(fTo,x1,maxX2,0);
      }

   else if(sendDir==D3Q27System::TSW)
   {
      fillData(fFrom,1,1,maxX3-1);
      distributeData(fTo,maxX1,maxX2,0);
   }
   else if(sendDir==D3Q27System::TSE)
   {
      fillData(fFrom,maxX1-1,1,maxX3-1);
      distributeData(fTo,0,maxX2,0);
   }
   else if(sendDir==D3Q27System::TNW)
   {
      fillData(fFrom,1,maxX2-1,maxX3-1);
      distributeData(fTo,maxX1,0,0);
   }
   else if(sendDir==D3Q27System::TNE)
   {
      fillData(fFrom,maxX1-1,maxX2-1,maxX3-1);
      distributeData(fTo,0,0,0);
   }
   else if(sendDir==D3Q27System::BSW)
   {
      fillData(fFrom,1,1,1);
      distributeData(fTo,maxX1,maxX2,maxX3);
   }
   else if(sendDir==D3Q27System::BSE)
   {
      fillData(fFrom,maxX1-1,1,1);
      distributeData(fTo,0,maxX2,maxX3);
   }
   else if(sendDir==D3Q27System::BNW)
   {
      fillData(fFrom,1,maxX2-1,1);
      distributeData(fTo,maxX1,0,maxX3);
   }
   else if(sendDir==D3Q27System::BNE)
   {
      fillData(fFrom,maxX1-1,maxX2-1,1);
      distributeData(fTo,0,0,maxX3);
   }
   else UB_THROW( UbException(UB_EXARGS,"unknown dir") );
}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETFullDirectConnector::fillData(EsoTwist3DPtr  fFrom, const int& x1, const int& x2, const int& x3)
{
   if(invStep)
   {
      fFrom->getDistribution(f,x1,x2,x3);
   }
   else
   {
      fFrom->getDistributionInv(f,x1,x2,x3);
   }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETFullDirectConnector::distributeData(EsoTwist3DPtr  fTo, const int& x1, const int& x2, const int& x3)
{
   if(invStep)
   {
      fTo->setDistributionInv(f,x1,x2,x3);
   }
   else
   {
      fTo->setDistribution(f,x1,x2,x3);
   }
}
//////////////////////////////////////////////////////////////////////////
