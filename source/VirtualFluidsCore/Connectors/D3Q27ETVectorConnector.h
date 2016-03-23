#ifndef D3Q27VECTORCONNECTOR_H
#define D3Q27VECTORCONNECTOR_H

#include <vector>

#include "basics/transmitter/TbTransmitter.h"
#include "basics/container/CbVector.h"
#include "basics/utilities/UbTiming.h"
#include "Block3DConnector.h"
#include "D3Q27System.h"
#include "Block3D.h"
#include "LBMKernelETD3Q27.h"
#include "EsoTwistD3Q27System.h"

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

class Block3D;

//daten werden in einen vector (dieser befindet sich im transmitter) kopiert
//der vector wird via transmitter uebertragen
//transmitter kann ein lokal, MPI, RCG, CTL oder was auch immer fuer ein
//transmitter sein, der von Transmitter abgeleitet ist ;-)
template< typename VectorTransmitter > //TbTransmitter< CbVector<LBMReal> >
class D3Q27ETVectorConnector : public Block3DConnector
{
public:
   //typedef CbVector<LBMReal> DataType;
   //typedef boost::shared_ptr< TbTransmitter< DataType > > TransmitterPtr;
   typedef typename VectorTransmitter::value_type  vector_type;
   typedef boost::shared_ptr< VectorTransmitter >  VectorTransmitterPtr;
public:
   D3Q27ETVectorConnector(  Block3DPtr block
                       , VectorTransmitterPtr sender
                       , VectorTransmitterPtr receiver
                       , int sendDir); 

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
   bool isInterpolationConnectorFC() { return false; }

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

   VectorTransmitterPtr sender;
   VectorTransmitterPtr receiver;
   UbTiming timer;
};

//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
D3Q27ETVectorConnector< VectorTransmitter >::D3Q27ETVectorConnector(   Block3DPtr block
                                               , VectorTransmitterPtr sender
                                               , VectorTransmitterPtr receiver
                                               , int sendDir )
                                               :  Block3DConnector(sendDir)
                                               , block(block)
                                               , sender(sender)
                                               , receiver(receiver)
{
   if(!block || !sender || !receiver) 
      UB_THROW( UbException(UB_EXARGS,"sender or receiver == NULL!!") );
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
bool D3Q27ETVectorConnector< VectorTransmitter >::isLocalConnector()  
{ 
   return !this->isRemoteConnector(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
bool D3Q27ETVectorConnector< VectorTransmitter >::isRemoteConnector() 
{ 
   return (   ( sender && sender->isRemoteTransmitter() ) 
      || ( receiver && receiver->isRemoteTransmitter() ) );
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETVectorConnector< VectorTransmitter >::sendTransmitterDataSize()
{ 
   assert(sender  !=NULL); sender->sendDataSize();        
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETVectorConnector< VectorTransmitter >::receiveTransmitterDataSize()
{ 
   assert(receiver!=NULL); receiver->receiveDataSize();   
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETVectorConnector< VectorTransmitter >::prepareForSend()            
{ 
   assert(sender  !=NULL); sender->prepareForSend();      
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETVectorConnector< VectorTransmitter >::sendVectors()                
{ 
   assert(sender  !=NULL); sender->sendData();            
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETVectorConnector< VectorTransmitter >::prepareForReceive()          
{ 
   assert(receiver!=NULL); receiver->prepareForReceive(); 
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETVectorConnector< VectorTransmitter >::receiveVectors()             
{ 
   assert(receiver!=NULL); receiver->receiveData();       
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETVectorConnector< VectorTransmitter >::init()
{
   int maxX1 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX1()-1;
   int maxX2 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX2()-1;
   int maxX3 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX3()-1;

   int anz = 9;
   switch(sendDir)
   {		                  
   case D3Q27System::ZERO : UB_THROW( UbException(UB_EXARGS,"ZERO not allowed") ); break;
   case D3Q27System::E    : 
   case D3Q27System::W    : sender->getData().resize(maxX2*maxX3*anz, 0.0);   break; 
   case D3Q27System::N    :
   case D3Q27System::S    : sender->getData().resize(maxX1*maxX3*anz, 0.0);   break; 
   case D3Q27System::T    : 
   case D3Q27System::B    : sender->getData().resize(maxX1*maxX2*anz, 0.0);   break; 

   case D3Q27System::NE  :   
   case D3Q27System::SW  :   
   case D3Q27System::SE  :   
   case D3Q27System::NW  :  sender->getData().resize(maxX3*3, 0.0);   break; 

   case D3Q27System::TE  :   
   case D3Q27System::BW  :   
   case D3Q27System::BE  :   
   case D3Q27System::TW  :  sender->getData().resize(maxX2*3, 0.0);   break; 

   case D3Q27System::TN  :   
   case D3Q27System::BS  :   
   case D3Q27System::BN  :   
   case D3Q27System::TS  :  sender->getData().resize(maxX1*3, 0.0);   break;

   case D3Q27System::TNE : 
   case D3Q27System::BSW : 
   case D3Q27System::BNE : 
   case D3Q27System::TSW : 
   case D3Q27System::TSE : 
   case D3Q27System::BNW : 
   case D3Q27System::BSE : 
   case D3Q27System::TNW :  sender->getData().resize(1, 0.0);   break;

   default: UB_THROW( UbException(UB_EXARGS,"unknown sendDir") );
   }
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETVectorConnector< VectorTransmitter >::fillSendVectors() 
{ 
   DistributionArray3DPtr  fFrom = block.lock()->getKernel()->getDataSet()->getFdistributions();
   int maxX1 = (int)fFrom->getNX1()-1;
   int maxX2 = (int)fFrom->getNX2()-1;
   int maxX3 = (int)fFrom->getNX3()-1;

   LBMReal f[D3Q27System::ENDF+1];
   int sendDirForEdge[3];

   vector_type& data = sender->getData();

   int index = 0;
   //EAST
   if(sendDir==D3Q27System::E)
   {
      for(int x3=1; x3<maxX3; x3++)   
      {
         for(int x2=1; x2<maxX2; x2++)   
         {
            fFrom->getDistributionInv(f,maxX1-1,x2,x3);
            data[index++] = f[D3Q27System::E]; 
            data[index++] = f[D3Q27System::NE]; 
            data[index++] = f[D3Q27System::SE];
            data[index++] = f[D3Q27System::TE]; 
            data[index++] = f[D3Q27System::BE];
            data[index++] = f[D3Q27System::TNE];
            data[index++] = f[D3Q27System::TSE];
            data[index++] = f[D3Q27System::BNE];
            data[index++] = f[D3Q27System::BSE];
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
            fFrom->getDistributionInv(f,1,x2,x3);
            data[index++] = f[D3Q27System::W ]; 
            data[index++] = f[D3Q27System::NW]; 
            data[index++] = f[D3Q27System::SW];
            data[index++] = f[D3Q27System::TW]; 
            data[index++] = f[D3Q27System::BW];
            data[index++] = f[D3Q27System::TNW];
            data[index++] = f[D3Q27System::TSW];
            data[index++] = f[D3Q27System::BNW];
            data[index++] = f[D3Q27System::BSW];
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
            fFrom->getDistributionInv(f,x1,maxX2-1,x3);
            data[index++] = f[D3Q27System::N ]; 
            data[index++] = f[D3Q27System::NE]; 
            data[index++] = f[D3Q27System::NW];
            data[index++] = f[D3Q27System::TN]; 
            data[index++] = f[D3Q27System::BN];
            data[index++] = f[D3Q27System::TNE];
            data[index++] = f[D3Q27System::TNW];
            data[index++] = f[D3Q27System::BNE];
            data[index++] = f[D3Q27System::BNW];
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
            fFrom->getDistributionInv(f,x1,1,x3);
            data[index++] = f[D3Q27System::S ]; 
            data[index++] = f[D3Q27System::SE]; 
            data[index++] = f[D3Q27System::SW];
            data[index++] = f[D3Q27System::TS]; 
            data[index++] = f[D3Q27System::BS];
            data[index++] = f[D3Q27System::TSE];
            data[index++] = f[D3Q27System::TSW];
            data[index++] = f[D3Q27System::BSE];
            data[index++] = f[D3Q27System::BSW];
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
            fFrom->getDistributionInv(f,x1,x2,maxX3-1);
            data[index++] = f[D3Q27System::T ]; 
            data[index++] = f[D3Q27System::TE]; 
            data[index++] = f[D3Q27System::TW];
            data[index++] = f[D3Q27System::TN]; 
            data[index++] = f[D3Q27System::TS];
            data[index++] = f[D3Q27System::TNE];
            data[index++] = f[D3Q27System::TNW];
            data[index++] = f[D3Q27System::TSE];
            data[index++] = f[D3Q27System::TSW];
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
            fFrom->getDistributionInv(f,x1,x2,1);
            data[index++] = f[D3Q27System::B ]; 
            data[index++] = f[D3Q27System::BE]; 
            data[index++] = f[D3Q27System::BW];
            data[index++] = f[D3Q27System::BN]; 
            data[index++] = f[D3Q27System::BS];
            data[index++] = f[D3Q27System::BNE];
            data[index++] = f[D3Q27System::BNW];
            data[index++] = f[D3Q27System::BSE];
            data[index++] = f[D3Q27System::BSW];
         }
      }
   }
   //NE NW SW SE
   else if(sendDir==D3Q27System::NE || sendDir==D3Q27System::NW || sendDir==D3Q27System::SW || sendDir==D3Q27System::SE)
   {
      int x1 = 0;
      int x2 = 0;
      switch(sendDir)  
      {
      case D3Q27System::NE:   x1=maxX1-1; x2=maxX2-1; 
         sendDirForEdge[0]=D3Q27System::NE; 
         sendDirForEdge[1]=D3Q27System::TNE;
         sendDirForEdge[2]=D3Q27System::BNE; 
         break;
      case D3Q27System::NW:   x1=1;       x2=maxX2-1; 
         sendDirForEdge[0]=D3Q27System::NW; 
         sendDirForEdge[1]=D3Q27System::TNW;
         sendDirForEdge[2]=D3Q27System::BNW;          
         break;
      case D3Q27System::SW:   x1=1;       x2=1;       
         sendDirForEdge[0]=D3Q27System::SW; 
         sendDirForEdge[1]=D3Q27System::TSW;
         sendDirForEdge[2]=D3Q27System::BSW; 
         break;
      case D3Q27System::SE:   x1=maxX1-1; x2=1;       
         sendDirForEdge[0]=D3Q27System::SE; 
         sendDirForEdge[1]=D3Q27System::TSE;
         sendDirForEdge[2]=D3Q27System::BSE; 
         break;
      }
      for(int x3=1; x3<maxX3; x3++)
      {
         data[index++] = fFrom->getDistributionInvForDirection(x1,x2,x3,sendDirForEdge[0]); 
         data[index++] = fFrom->getDistributionInvForDirection(x1,x2,x3,sendDirForEdge[1]); 
         data[index++] = fFrom->getDistributionInvForDirection(x1,x2,x3,sendDirForEdge[2]); 
      }
   }
   //TE TW BW BE
   else if(sendDir==D3Q27System::TE || sendDir==D3Q27System::TW || sendDir==D3Q27System::BW || sendDir==D3Q27System::BE)
   {
      int x1 = 0;
      int x3 = 0;
      switch(sendDir)  
      {
      case D3Q27System::TE:   x1=maxX1-1; x3=maxX3-1; 
         sendDirForEdge[0]=D3Q27System::TE; 
         sendDirForEdge[1]=D3Q27System::TNE;
         sendDirForEdge[2]=D3Q27System::TSE; 
         break;
      case D3Q27System::TW:   x1=1;       x3=maxX3-1; 
         sendDirForEdge[0]=D3Q27System::TW; 
         sendDirForEdge[1]=D3Q27System::TNW;
         sendDirForEdge[2]=D3Q27System::TSW; 
         break;
      case D3Q27System::BW:   x1=1;       x3=1;       
         sendDirForEdge[0]=D3Q27System::BW; 
         sendDirForEdge[1]=D3Q27System::BNW;
         sendDirForEdge[2]=D3Q27System::BSW; 
         break;
      case D3Q27System::BE:   x1=maxX1-1; x3=1;      
         sendDirForEdge[0]=D3Q27System::BE; 
         sendDirForEdge[1]=D3Q27System::BNE;
         sendDirForEdge[2]=D3Q27System::BSE; 
         break;
      }
      for(int x2=1; x2<maxX2; x2++) 
      {
         data[index++] = fFrom->getDistributionInvForDirection(x1,x2,x3,sendDirForEdge[0]); 
         data[index++] = fFrom->getDistributionInvForDirection(x1,x2,x3,sendDirForEdge[1]); 
         data[index++] = fFrom->getDistributionInvForDirection(x1,x2,x3,sendDirForEdge[2]); 
      }
   }
   //TN BN BS TS
   else if(sendDir==D3Q27System::TN || sendDir==D3Q27System::BN || sendDir==D3Q27System::BS || sendDir==D3Q27System::TS)
   {
      int x2 = 0;
      int x3 = 0;
      switch(sendDir)  
      {
      case D3Q27System::TN:   x3=maxX3-1; x2=maxX2-1; 
         sendDirForEdge[0]=D3Q27System::TN; 
         sendDirForEdge[1]=D3Q27System::TNE;
         sendDirForEdge[2]=D3Q27System::TNW; 
         break;
      case D3Q27System::BN:   x3=1;       x2=maxX2-1; 
         sendDirForEdge[0]=D3Q27System::BN; 
         sendDirForEdge[1]=D3Q27System::BNE;
         sendDirForEdge[2]=D3Q27System::BNW; 
         break;
      case D3Q27System::BS:   x3=1;       x2=1;       
         sendDirForEdge[0]=D3Q27System::BS; 
         sendDirForEdge[1]=D3Q27System::BSE;
         sendDirForEdge[2]=D3Q27System::BSW; 
         break;
      case D3Q27System::TS:   x3=maxX3-1; x2=1;       
         sendDirForEdge[0]=D3Q27System::TS; 
         sendDirForEdge[1]=D3Q27System::TSE;
         sendDirForEdge[2]=D3Q27System::TSW; 
         break;
      }
      for(int x1=1; x1<maxX1; x1++)
      {
         data[index++] = fFrom->getDistributionInvForDirection(x1,x2,x3,sendDirForEdge[0]); 
         data[index++] = fFrom->getDistributionInvForDirection(x1,x2,x3,sendDirForEdge[1]); 
         data[index++] = fFrom->getDistributionInvForDirection(x1,x2,x3,sendDirForEdge[2]);
      }
   }
   //TNE TNW TSW TSE BNE BNW BSW BSE
   else if(sendDir==D3Q27System::TNE || sendDir==D3Q27System::TNW || sendDir==D3Q27System::TSW || sendDir==D3Q27System::TSE
      || sendDir==D3Q27System::BNE || sendDir==D3Q27System::BNW || sendDir==D3Q27System::BSW || sendDir==D3Q27System::BSE)
   {
      int x1 = 0;
      int x2 = 0;
      int x3 = 0;
      switch(sendDir) 
      {
      case D3Q27System::TNE:   x1=maxX1-1; x2=maxX2-1; x3=maxX3-1; break;
      case D3Q27System::TNW:   x1=1;       x2=maxX2-1; x3=maxX3-1; break;
      case D3Q27System::TSW:   x1=1;       x2=1;       x3=maxX3-1; break;
      case D3Q27System::TSE:   x1=maxX1-1; x2=1;       x3=maxX3-1; break;
      case D3Q27System::BNE:   x1=maxX1-1; x2=maxX2-1; x3=1;       break;
      case D3Q27System::BNW:   x1=1;       x2=maxX2-1; x3=1;       break;
      case D3Q27System::BSW:   x1=1;       x2=1;       x3=1;       break;
      case D3Q27System::BSE:   x1=maxX1-1; x2=1;       x3=1;       break;
      }
      data[index++] = fFrom->getDistributionInvForDirection(x1,x2,x3,sendDir);
   }
   else UB_THROW( UbException(UB_EXARGS,"unknown dir") );
}
////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter >
void D3Q27ETVectorConnector< VectorTransmitter >::distributeReceiveVectors()
{
   /*e.g. connector sendet nach EAST --> empfaengt daten aus WEST ;-)*/
   DistributionArray3DPtr  fTo = block.lock()->getKernel()->getDataSet()->getFdistributions();
   int maxX1 = (int)fTo->getNX1()-1;
   int maxX2 = (int)fTo->getNX2()-1;
   int maxX3 = (int)fTo->getNX3()-1;
   int index = 0;
   int sendDirForEdge[3];

   vector_type& data = receiver->getData();

   if(sendDir==D3Q27System::W)
   {
      for(int x3=1; x3<maxX3; x3++)  
      {
         for(int x2=1; x2<maxX2; x2++)   
         {
            fTo->setDistributionForDirection(data[index++],0,x2,x3,D3Q27System::E);
            fTo->setDistributionForDirection(data[index++],0,x2,x3,D3Q27System::NE); 
            fTo->setDistributionForDirection(data[index++],0,x2,x3,D3Q27System::SE);
            fTo->setDistributionForDirection(data[index++],0,x2,x3,D3Q27System::TE); 
            fTo->setDistributionForDirection(data[index++],0,x2,x3,D3Q27System::BE);
            fTo->setDistributionForDirection(data[index++],0,x2,x3,D3Q27System::TNE);
            fTo->setDistributionForDirection(data[index++],0,x2,x3,D3Q27System::TSE);
            fTo->setDistributionForDirection(data[index++],0,x2,x3,D3Q27System::BNE);
            fTo->setDistributionForDirection(data[index++],0,x2,x3,D3Q27System::BSE);
         }
      }
   }
   else if(sendDir==D3Q27System::E)
   {
      for(int x3=1; x3<maxX3; x3++)  
      {
         for(int x2=1; x2<maxX2; x2++)   
         {
            fTo->setDistributionForDirection(data[index++],maxX1,x2,x3,D3Q27System::W );
            fTo->setDistributionForDirection(data[index++],maxX1,x2,x3,D3Q27System::NW);
            fTo->setDistributionForDirection(data[index++],maxX1,x2,x3,D3Q27System::SW);
            fTo->setDistributionForDirection(data[index++],maxX1,x2,x3,D3Q27System::TW);
            fTo->setDistributionForDirection(data[index++],maxX1,x2,x3,D3Q27System::BW); 
            fTo->setDistributionForDirection(data[index++],maxX1,x2,x3,D3Q27System::TNW);
            fTo->setDistributionForDirection(data[index++],maxX1,x2,x3,D3Q27System::TSW);
            fTo->setDistributionForDirection(data[index++],maxX1,x2,x3,D3Q27System::BNW);
            fTo->setDistributionForDirection(data[index++],maxX1,x2,x3,D3Q27System::BSW);
         }
      }
   }
   else if(sendDir==D3Q27System::S)
   {
      for(int x3=1; x3<maxX3; x3++)  
      {
         for(int x1=1; x1<maxX1; x1++)   
         {
            fTo->setDistributionForDirection(data[index++],x1,0,x3,D3Q27System::N );
            fTo->setDistributionForDirection(data[index++],x1,0,x3,D3Q27System::NE);
            fTo->setDistributionForDirection(data[index++],x1,0,x3,D3Q27System::NW);
            fTo->setDistributionForDirection(data[index++],x1,0,x3,D3Q27System::TN);
            fTo->setDistributionForDirection(data[index++],x1,0,x3,D3Q27System::BN);
            fTo->setDistributionForDirection(data[index++],x1,0,x3,D3Q27System::TNE);
            fTo->setDistributionForDirection(data[index++],x1,0,x3,D3Q27System::TNW);
            fTo->setDistributionForDirection(data[index++],x1,0,x3,D3Q27System::BNE);
            fTo->setDistributionForDirection(data[index++],x1,0,x3,D3Q27System::BNW);
         }
      }
   }
   else if(sendDir==D3Q27System::N)
   {
      for(int x3=1; x3<maxX3; x3++)  
      {
         for(int x1=1; x1<maxX1; x1++)    
         {
            fTo->setDistributionForDirection(data[index++],x1,maxX2,x3,D3Q27System::S );
            fTo->setDistributionForDirection(data[index++],x1,maxX2,x3,D3Q27System::SE); 
            fTo->setDistributionForDirection(data[index++],x1,maxX2,x3,D3Q27System::SW); 
            fTo->setDistributionForDirection(data[index++],x1,maxX2,x3,D3Q27System::TS); 
            fTo->setDistributionForDirection(data[index++],x1,maxX2,x3,D3Q27System::BS);
            fTo->setDistributionForDirection(data[index++],x1,maxX2,x3,D3Q27System::TSE);
            fTo->setDistributionForDirection(data[index++],x1,maxX2,x3,D3Q27System::TSW);
            fTo->setDistributionForDirection(data[index++],x1,maxX2,x3,D3Q27System::BSE);
            fTo->setDistributionForDirection(data[index++],x1,maxX2,x3,D3Q27System::BSW);
         }
      }
   }
   else if(sendDir==D3Q27System::B)
   {
      for(int x2=1; x2<maxX2; x2++)  
      {
         for(int x1=1; x1<maxX1; x1++)   
         {
            fTo->setDistributionForDirection(data[index++],x1,x2,0,D3Q27System::T );
            fTo->setDistributionForDirection(data[index++],x1,x2,0,D3Q27System::TE);
            fTo->setDistributionForDirection(data[index++],x1,x2,0,D3Q27System::TW);
            fTo->setDistributionForDirection(data[index++],x1,x2,0,D3Q27System::TN); 
            fTo->setDistributionForDirection(data[index++],x1,x2,0,D3Q27System::TS);
            fTo->setDistributionForDirection(data[index++],x1,x2,0,D3Q27System::TNE);
            fTo->setDistributionForDirection(data[index++],x1,x2,0,D3Q27System::TNW);
            fTo->setDistributionForDirection(data[index++],x1,x2,0,D3Q27System::TSE);
            fTo->setDistributionForDirection(data[index++],x1,x2,0,D3Q27System::TSW);
         }
      }
   }
   else if(sendDir==D3Q27System::T)
   {
      for(int x2=1; x2<maxX2; x2++)  
      {
         for(int x1=1; x1<maxX1; x1++)   
         {
            fTo->setDistributionForDirection(data[index++],x1,x2,maxX3,D3Q27System::B );
            fTo->setDistributionForDirection(data[index++],x1,x2,maxX3,D3Q27System::BE);
            fTo->setDistributionForDirection(data[index++],x1,x2,maxX3,D3Q27System::BW);
            fTo->setDistributionForDirection(data[index++],x1,x2,maxX3,D3Q27System::BN);
            fTo->setDistributionForDirection(data[index++],x1,x2,maxX3,D3Q27System::BS);
            fTo->setDistributionForDirection(data[index++],x1,x2,maxX3,D3Q27System::BNE);
            fTo->setDistributionForDirection(data[index++],x1,x2,maxX3,D3Q27System::BNW);
            fTo->setDistributionForDirection(data[index++],x1,x2,maxX3,D3Q27System::BSE);
            fTo->setDistributionForDirection(data[index++],x1,x2,maxX3,D3Q27System::BSW);
         }
      }
   }
   //NE NW SW SE
   else if(sendDir==D3Q27System::NE || sendDir==D3Q27System::NW || sendDir==D3Q27System::SW || sendDir==D3Q27System::SE)
   {
      int inversDir = D3Q27System::getInvertDirection(sendDir);
      int x1 = 0;
      int x2 = 0;
      switch(sendDir)  //wenn sendir NE dann kommen werte von SW
      {
      case D3Q27System::NE:   x1=maxX1; x2=maxX2; 
         sendDirForEdge[0]=D3Q27System::SW; 
         sendDirForEdge[1]=D3Q27System::TSW;
         sendDirForEdge[2]=D3Q27System::BSW;
         break;
      case D3Q27System::NW:   x1=0;     x2=maxX2; 
         sendDirForEdge[0]=D3Q27System::SE; 
         sendDirForEdge[1]=D3Q27System::TSE;
         sendDirForEdge[2]=D3Q27System::BSE;        
         break;
      case D3Q27System::SW:   x1=0;     x2=0;       
         sendDirForEdge[0]=D3Q27System::NE; 
         sendDirForEdge[1]=D3Q27System::TNE;
         sendDirForEdge[2]=D3Q27System::BNE;  
         break;
      case D3Q27System::SE:   x1=maxX1; x2=0;  
         sendDirForEdge[0]=D3Q27System::NW; 
         sendDirForEdge[1]=D3Q27System::TNW;
         sendDirForEdge[2]=D3Q27System::BNW; 
         break;
      }
      for(int x3=1; x3<maxX3; x3++)
      {
         fTo->setDistributionForDirection(data[index++],x1, x2, x3, sendDirForEdge[0]);
         fTo->setDistributionForDirection(data[index++],x1, x2, x3, sendDirForEdge[1]);
         fTo->setDistributionForDirection(data[index++],x1, x2, x3, sendDirForEdge[2]);
      }

   }
   //TE TW BW BE
   else if(sendDir==D3Q27System::TE || sendDir==D3Q27System::TW || sendDir==D3Q27System::BW || sendDir==D3Q27System::BE)

   {
      int inversDir = D3Q27System::getInvertDirection(sendDir);
      int x1 = 0;
      int x3 = 0;
      switch(sendDir)  //wenn sendir NE dann kommen werte von SW
      {
      case D3Q27System::TE:   x1=maxX1; x3=maxX3; 
         sendDirForEdge[0]=D3Q27System::BW; 
         sendDirForEdge[1]=D3Q27System::BNW;
         sendDirForEdge[2]=D3Q27System::BSW;
         break;
      case D3Q27System::TW:   x1=0;     x3=maxX3; 
         sendDirForEdge[0]=D3Q27System::BE; 
         sendDirForEdge[1]=D3Q27System::BNE;
         sendDirForEdge[2]=D3Q27System::BSE; 
         break;
      case D3Q27System::BW:   x1=0;     x3=0;       
         sendDirForEdge[0]=D3Q27System::TE; 
         sendDirForEdge[1]=D3Q27System::TNE;
         sendDirForEdge[2]=D3Q27System::TSE; 
         break;
      case D3Q27System::BE:   x1=maxX1; x3=0;      
         sendDirForEdge[0]=D3Q27System::TW; 
         sendDirForEdge[1]=D3Q27System::TNW;
         sendDirForEdge[2]=D3Q27System::TSW; 
         break;
      }
      for(int x2=1; x2<maxX2; x2++) 
      {  
         fTo->setDistributionForDirection(data[index++],x1, x2, x3, sendDirForEdge[0]);
         fTo->setDistributionForDirection(data[index++],x1, x2, x3, sendDirForEdge[1]);
         fTo->setDistributionForDirection(data[index++],x1, x2, x3, sendDirForEdge[2]);
      }
   }
   //TN BN BS TS
   else if(sendDir==D3Q27System::TN || sendDir==D3Q27System::BN || sendDir==D3Q27System::BS || sendDir==D3Q27System::TS)
   {
      int inversDir = D3Q27System::getInvertDirection(sendDir);
      int x2 = 0;
      int x3 = 0;
      switch(sendDir)  
      {
      case D3Q27System::TN:   x3=maxX3; x2=maxX2; 
         sendDirForEdge[0]=D3Q27System::BS; 
         sendDirForEdge[1]=D3Q27System::BSE;
         sendDirForEdge[2]=D3Q27System::BSW;
         break;
      case D3Q27System::BN:   x3=0;       x2=maxX2; 
         sendDirForEdge[0]=D3Q27System::TS; 
         sendDirForEdge[1]=D3Q27System::TSE;
         sendDirForEdge[2]=D3Q27System::TSW; 
         break;
      case D3Q27System::BS:   x3=0;       x2=0;  
         sendDirForEdge[0]=D3Q27System::TN; 
         sendDirForEdge[1]=D3Q27System::TNE;
         sendDirForEdge[2]=D3Q27System::TNW; 
         break;
      case D3Q27System::TS:   x3=maxX3; x2=0;       
         sendDirForEdge[0]=D3Q27System::BN; 
         sendDirForEdge[1]=D3Q27System::BNE;
         sendDirForEdge[2]=D3Q27System::BNW; 
         break;

      }
      for(int x1=1; x1<maxX1; x1++) 
      {  
         fTo->setDistributionForDirection(data[index++],x1, x2, x3, sendDirForEdge[0]);
         fTo->setDistributionForDirection(data[index++],x1, x2, x3, sendDirForEdge[1]);
         fTo->setDistributionForDirection(data[index++],x1, x2, x3, sendDirForEdge[2]);
      }
   }
   //TNE TNW TSW TSE BNE BNW BSW BSE
   else if(sendDir==D3Q27System::TNE || sendDir==D3Q27System::TNW || sendDir==D3Q27System::TSW || sendDir==D3Q27System::TSE
      || sendDir==D3Q27System::BNE || sendDir==D3Q27System::BNW || sendDir==D3Q27System::BSW || sendDir==D3Q27System::BSE)
   {
      int inversDir = D3Q27System::getInvertDirection(sendDir);
      int x1 = 0;
      int x2 = 0;
      int x3 = 0;
      unsigned long int etDir = 0;
      switch(sendDir) 
      {
      case D3Q27System::TNE:   x1=maxX1; x2=maxX2; x3=maxX3; etDir=EsoTwistD3Q27System::etTNE; break;
      case D3Q27System::TNW:   x1=0; x2=maxX2; x3=maxX3; etDir=EsoTwistD3Q27System::etTNW; break;
      case D3Q27System::TSW:   x1=0; x2=0; x3=maxX3; etDir=EsoTwistD3Q27System::etTSW; break;
      case D3Q27System::TSE:   x1=maxX1; x2=0; x3=maxX3; etDir=EsoTwistD3Q27System::etTSE; break;
      case D3Q27System::BNE:   x1=maxX1; x2=maxX2; x3=0; etDir=EsoTwistD3Q27System::etBNE; break;
      case D3Q27System::BNW:   x1=0; x2=maxX2; x3=0; etDir=EsoTwistD3Q27System::etBNW; break;
      case D3Q27System::BSW:   x1=0; x2=0; x3=0; etDir=EsoTwistD3Q27System::etBSW; break;
      case D3Q27System::BSE:   x1=maxX1; x2=0; x3=0; etDir=EsoTwistD3Q27System::etBSE; break;
      }
      fTo->setDistributionForDirection(data[index++],x1, x2, x3, inversDir);
   }
   else UB_THROW( UbException(UB_EXARGS,"unknown dir") );
}
//////////////////////////////////////////////////////////////////////////
template< typename VectorTransmitter > 
double D3Q27ETVectorConnector<VectorTransmitter>::getSendRecieveTime()
{
   return 0;
}

#endif //D3Q27VECTORCONNECTOR_H
