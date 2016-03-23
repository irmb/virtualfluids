#include "D3Q27ETFullDirectConnector2.h"
#include "LBMKernelETD3Q27.h"
#include "D3Q27EsoTwist3DSplittedVector.h"


using namespace std;

D3Q27ETFullDirectConnector2::D3Q27ETFullDirectConnector2(Block3DPtr from, Block3DPtr to, int sendDir) 
   : LocalBlock3DConnector(from, to, sendDir)

{

}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETFullDirectConnector2::init()
{
   maxX1 = (int)this->from.lock()->getKernel()->getDataSet()->getFdistributions()->getNX1()-1;
   maxX2 = (int)this->from.lock()->getKernel()->getDataSet()->getFdistributions()->getNX2()-1;
   maxX3 = (int)this->from.lock()->getKernel()->getDataSet()->getFdistributions()->getNX3()-1;

   fFrom =  boost::dynamic_pointer_cast<EsoTwist3D>(from.lock()->getKernel()->getDataSet()->getFdistributions());
   fTo   = boost::dynamic_pointer_cast<EsoTwist3D>(to.lock()->getKernel()->getDataSet()->getFdistributions());

   //localDistributionsFrom = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->from.lock()->getKernel()->getDataSet()->getFdistributions())->getLocalDistributions();
   //nonLocalDistributionsFrom = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->from.lock()->getKernel()->getDataSet()->getFdistributions())->getNonLocalDistributions();
   //zeroDistributionsFrom = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->from.lock()->getKernel()->getDataSet()->getFdistributions())->getZeroDistributions();

   //localDistributionsTo = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->to.lock()->getKernel()->getDataSet()->getFdistributions())->getLocalDistributions();
   //nonLocalDistributionsTo = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->to.lock()->getKernel()->getDataSet()->getFdistributions())->getNonLocalDistributions();
   //zeroDistributionsTo = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->to.lock()->getKernel()->getDataSet()->getFdistributions())->getZeroDistributions();
}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETFullDirectConnector2::sendVectors()
{
   //EsoTwist3DPtr  fFrom =  boost::dynamic_pointer_cast<EsoTwist3D>(from.lock()->getKernel()->getDataSet()->getFdistributions());
   //EsoTwist3DPtr  fTo   = boost::dynamic_pointer_cast<EsoTwist3D>(to.lock()->getKernel()->getDataSet()->getFdistributions());

   //int maxX1 = (int)fFrom->getNX1()-1;
   //int maxX2 = (int)fFrom->getNX2()-1;
   //int maxX3 = (int)fFrom->getNX3()-1;

   //LBMReal f[D3Q27System::ENDF+1];

   //localDistributionsFrom    = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->from.lock()->getKernel()->getDataSet()->getFdistributions())->getLocalDistributions();
   //nonLocalDistributionsFrom = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->from.lock()->getKernel()->getDataSet()->getFdistributions())->getNonLocalDistributions();
   //zeroDistributionsFrom     = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->from.lock()->getKernel()->getDataSet()->getFdistributions())->getZeroDistributions();

   //localDistributionsTo      = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->to.lock()->getKernel()->getDataSet()->getFdistributions())->getLocalDistributions();
   //nonLocalDistributionsTo   = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->to.lock()->getKernel()->getDataSet()->getFdistributions())->getNonLocalDistributions();
   //zeroDistributionsTo       = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->to.lock()->getKernel()->getDataSet()->getFdistributions())->getZeroDistributions();

   localDistributionsFrom = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->fFrom)->getLocalDistributions();
   nonLocalDistributionsFrom = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->fFrom)->getNonLocalDistributions();
   zeroDistributionsFrom = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->fFrom)->getZeroDistributions();

   localDistributionsTo = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->fTo)->getLocalDistributions();
   nonLocalDistributionsTo = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->fTo)->getNonLocalDistributions();
   zeroDistributionsTo = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(this->fTo)->getZeroDistributions();

   //EAST
   if(sendDir==D3Q27System::E)
   {
      for(int x3=1; x3<maxX3; x3++)   
      {
         for(int x2=1; x2<maxX2; x2++)   
         {
            fillData(maxX1-1,x2,x3);
            distributeData(0,x2,x3);
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
            fillData(1,x2,x3);
            distributeData(maxX1,x2,x3);
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
            fillData(x1,maxX2-1,x3);
            distributeData(x1,0,x3);
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
            fillData(x1,1,x3);
            distributeData(x1,maxX2,x3);
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
            fillData(x1,x2,maxX3-1);
            distributeData(x1,x2,0);
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
            fillData(x1,x2,1);
            distributeData(x1,x2,maxX3);
         }
      }
   }
   //NORTHEAST
   else if(sendDir==D3Q27System::NE)
   {   
      for(int x3=1; x3<maxX3; x3++)   
      {
         fillData(maxX1-1,maxX2-1,x3);
         distributeData(0,0,x3);
      }
   }
   //NORTHWEST
   else if(sendDir==D3Q27System::NW)
   {   
      for(int x3=1; x3<maxX3; x3++)   
      {
         fillData(1,maxX2-1,x3);
         distributeData(maxX1,0,x3);
      }
   }
   //SOUTHWEST
   else if(sendDir==D3Q27System::SW)
   {  
      for(int x3=1; x3<maxX3; x3++)   
      {
         fillData(1,1,x3);
         distributeData(maxX1,maxX2,x3);
      }
   }
   //SOUTHEAST
   else if(sendDir==D3Q27System::SE)
   {   
      for(int x3=1; x3<maxX3; x3++)   
      {
         fillData(maxX1-1,1,x3);
         distributeData(0,maxX2,x3);
      }
   }
   else if(sendDir==D3Q27System::TE)
      for(int x2=1; x2<maxX2; x2++)
      {
         fillData(maxX1-1,x2,maxX3-1);
         distributeData(0,x2,0);
      }
   else if(sendDir==D3Q27System::BW)
      for(int x2=1; x2<maxX2; x2++)
      {
         fillData(1,x2,1);
         distributeData(maxX1,x2,maxX3);
      }
   else if(sendDir==D3Q27System::BE)
      for(int x2=1; x2<maxX2; x2++)
      {
         fillData(maxX1-1,x2,1);
         distributeData(0,x2,maxX3);
      }
   else if(sendDir==D3Q27System::TW)
      for(int x2=1; x2<maxX2; x2++)
      {
         fillData(1,x2,maxX3-1);
         distributeData(maxX1,x2,0);
      }
   else if(sendDir==D3Q27System::TN)
      for(int x1=1; x1<maxX1; x1++)
      {
         fillData(x1,maxX2-1,maxX3-1);
         distributeData(x1,0,0);
      }
   else if(sendDir==D3Q27System::BS)
      for(int x1=1; x1<maxX1; x1++)
      {
         fillData(x1,1,1);
         distributeData(x1,maxX2,maxX3);
      }
   else if(sendDir==D3Q27System::BN)
      for(int x1=1; x1<maxX1; x1++)
      {
         fillData(x1,maxX2-1,1);
         distributeData(x1,0,maxX3);
      }

   else if(sendDir==D3Q27System::TS)
      for(int x1=1; x1<maxX1; x1++)
      {
         fillData(x1,1,maxX3-1);
         distributeData(x1,maxX2,0);
      }

   else if(sendDir==D3Q27System::TSW)
   {
      fillData(1,1,maxX3-1);
      distributeData(maxX1,maxX2,0);
   }
   else if(sendDir==D3Q27System::TSE)
   {
      fillData(maxX1-1,1,maxX3-1);
      distributeData(0,maxX2,0);
   }
   else if(sendDir==D3Q27System::TNW)
   {
      fillData(1,maxX2-1,maxX3-1);
      distributeData(maxX1,0,0);
   }
   else if(sendDir==D3Q27System::TNE)
   {
      fillData(maxX1-1,maxX2-1,maxX3-1);
      distributeData(0,0,0);
   }
   else if(sendDir==D3Q27System::BSW)
   {
      fillData(1,1,1);
      distributeData(maxX1,maxX2,maxX3);
   }
   else if(sendDir==D3Q27System::BSE)
   {
      fillData(maxX1-1,1,1);
      distributeData(0,maxX2,maxX3);
   }
   else if(sendDir==D3Q27System::BNW)
   {
      fillData(1,maxX2-1,1);
      distributeData(maxX1,0,maxX3);
   }
   else if(sendDir==D3Q27System::BNE)
   {
      fillData(maxX1-1,maxX2-1,1);
      distributeData(0,0,maxX3);
   }
   else UB_THROW( UbException(UB_EXARGS,"unknown dir") );
}
//////////////////////////////////////////////////////////////////////////
//void D3Q27ETFullDirectConnector2::fillData(/*EsoTwist3DPtr  fFrom,*/ int x1, int x2, int x3)
//{
//   if(rfStep)
//   {
//      fFrom->getDistribution(f,x1,x2,x3);
//   }
//   else
//   {
//      fFrom->getDistributionInv(f,x1,x2,x3);
//   }
//}
////////////////////////////////////////////////////////////////////////////
//void D3Q27ETFullDirectConnector2::distributeData(/*EsoTwist3DPtr  fTo,*/ int x1, int x2, int x3)
//{
//   if(rfStep)
//   {
//      fTo->setDistributionInv(f,x1,x2,x3);
//   }
//   else
//   {
//      fTo->setDistribution(f,x1,x2,x3);
//   }
//}
//////////////////////////////////////////////////////////////////////////
//double D3Q27ETFullDirectConnector2::getSendRecieveTime()
//{
//   return 0;
//}
