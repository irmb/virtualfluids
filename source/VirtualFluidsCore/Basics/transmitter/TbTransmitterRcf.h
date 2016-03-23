//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef TBTRANSMITTERRCF_H
#define TBTRANSMITTERRCF_H

/*=========================================================================*/
/*  RCF Transmitter                                                        */
/*                                                                         */
/**
This Class provides the base for exception handling.
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 08.11.07
*/ 

/*
usage: ...
*/

#include <iostream>
#include <vector>
#include <map>

#include <basics/transmitter/TbTransmitter.h>
#include <basics/container/CbVector.h>
#include <basics/utilities/UbLogger.h>

#ifdef CAB_RCF
//////////////////////////////////////////////////////////////////////////
// RCF STUFF
//////////////////////////////////////////////////////////////////////////

#include <RCF/Idl.hpp>
#include <RCF/TcpEndpoint.hpp>
#include <boost/shared_ptr.hpp>

#include <3rdParty/rcf/RcfSerializationIncludes.h>
#include <3rdParty/rcf/RcfSystem.h>
#include <3rdParty/rcf/IRcfIpService.h>
#include <3rdParty/rcf/RcfConnection.h>

//zum ausstausch mittels RCF transmitter:
RCF_BEGIN(IRcfTransmitterReceiverService, "IRcfTransmitterReceiverService")
   RCF_METHOD_V2(void, receiveVectorForTransmitter,int /*tag*/, const std::vector<double>& /*data*/);
   RCF_METHOD_V2(void, receiveVectorForTransmitter,int /*tag*/, const CbVector<double>& /*data*/);
   RCF_METHOD_V2(void, receiveVectorForTransmitter,int /*tag*/, const std::vector<float>& /*data*/);
   RCF_METHOD_V2(void, receiveVectorForTransmitter,int /*tag*/, const CbVector<float>& /*data*/);
RCF_END(IRcfTransmitterReceiverService);

//////////////////////////////////////////////////////////////////////////
// TbVectorSenderRcf< Vector > 
//////////////////////////////////////////////////////////////////////////
template< typename Vector >
class TbVectorSenderRcf : public TbTransmitter< Vector >
{
public:
   typedef Vector value_type;   
   
   //static members
private:
   static std::vector< RcfClient<IRcfTransmitterReceiverService> >  recvServices;

   static void setRcfClients(const std::string& recvServiceID);

public:
   TbVectorSenderRcf(const std::string& recvServiceName,const RcfConnection& receiveProcess, const int& tag) 
      : TbTransmitter< Vector >()
   { 
      if( recvServices.empty() ) setRcfClients(recvServiceName);
      this->receiveRank	   = receiveProcess.getRank();
      this->tag            = tag;
      this->receiveProcess = receiveProcess;
   }

   bool isLocalTransmitter()  const { return false;                        }
   bool isRemoteTransmitter() const { return !this->isLocalTransmitter();  }

   void sendDataSize()    { this->sendData(); }
   void receiveDataSize() { UB_THROW( UbException(UB_EXARGS,"TbRcfVectorSender sends only") ); } 

   void sendData() 
   { 
      //remote prozess=client erhaelt daten
      recvServices[receiveRank].receiveVectorForTransmitter(tag, TbTransmitter< Vector >::getData());
   }
   void prepareForReceive() { UB_THROW( UbException(UB_EXARGS,"TbRcfVectorSender sends only") ); }
   Vector& receiveData()    { UB_THROW( UbException(UB_EXARGS,"TbRcfVectorSender sends only") ); }

   //info-section (usable for remote transmitter)
   int  getSendTbRank()    const { return  this->receiveRank; /*=der rank an den gesendet wird*/}
   int  getSendTbTag()     const { return  this->tag;                                           }
   int  getRecvFromRank()  const { UB_THROW( UbException(UB_EXARGS,"TbVectorSenderRcf sends only") ); }
   int  getRecvFromTag()   const { UB_THROW( UbException(UB_EXARGS,"TbVectorSenderRcf sends only") ); }

   std::string toString() const { return "TbVectorSenderRcf< "+(std::string)typeid(Vector).name()+"<"+(std::string)typeid(typename Vector::value_type).name() +"> > to rank (tag)"+UbSystem::toString(receiveRank)+"("+UbSystem::toString(tag)+") connection "+receiveProcess.toString(); }

private:
   RcfConnection receiveProcess;
   int tag;
   int receiveRank;
};

//////////////////////////////////////////////////////////////////////////
template< typename Vector >
std::vector< RcfClient<IRcfTransmitterReceiverService> >  TbVectorSenderRcf< Vector >::recvServices;

template< typename Vector >
void TbVectorSenderRcf< Vector >::setRcfClients(const std::string& recvServiceID)
{
   UBLOG(logINFO,"invoked setRcfClients");
   RcfConnection ipConnection = RcfSystem::getIpServiceConnection();
   if(!ipConnection) UB_THROW( UbException(UB_EXARGS,"ups, no IpServiceConnection") );
   RcfClient<IRcfIpService> ipService(RCF::TcpEndpoint(ipConnection.getIp(), ipConnection.getPort()));

   //Mit RecvServices verbinden
   std::vector<RcfConnection> connections = ipService.getServiceConnections(recvServiceID);
   if(connections.empty()) 
      UB_THROW( UbException(UB_EXARGS,"no existing RecvService with ID = "+recvServiceID) );
   std::sort(connections.begin(),connections.end(),RcfConnection::compareRank());

   for(std::size_t i=0; i<connections.size(); i++)
   {
      if( (int)i != connections[i].getRank() )
         UB_THROW( UbException(UB_EXARGS,"recvServices must have continougs ranks sarting from 0") );
      recvServices.push_back(RcfClient<IRcfTransmitterReceiverService>(RCF::TcpEndpoint(connections[i].getIp(), connections[i].getPort())) );
   }
}


//////////////////////////////////////////////////////////////////////////
// TbVectorReceiverRcf
//////////////////////////////////////////////////////////////////////////
template< typename Vector >
class TbVectorReceiverRcf : public TbTransmitter< Vector >
{
   typedef std::map<int, TbVectorReceiverRcf< Vector >* >  ReceiverMap;
   typedef typename ReceiverMap::iterator                  ReceiverMapIt;
   
public:
   typedef Vector value_type;

   //static members
private:
   static ReceiverMap  receiverMap;
   static boost::mutex staticReceiverMapMutex;

public:
   static void receiveVectorForTransmitter(int tag, const Vector& data)
   {
      TbVectorReceiverRcf< Vector >* receiver = NULL;
      
      //receiver ermitteln (map nicht thread-safe->lock! aber nur kurz, ansonsten ab einer gewissen anzahl an clients/blöcken->deadlock)
      //vermutlich brauch man den nicht mal... (denn das registrieren sollte zu diesem zeitpunkt abgeschlossen sein)
      //allerdings sollte man nicht gleichzeitig Suchoperationen in ner nich thread sicheren map durchführen!!!
      {
         boost::mutex::scoped_lock lock(staticReceiverMapMutex); //wenn man das ausserhalb macht, gibt es ab einer bestimmten anzahl an clients/bloecke einen deadlock
         
         ReceiverMapIt result = TbVectorReceiverRcf< Vector >::receiverMap.find(tag);
         if( result == TbVectorReceiverRcf< Vector >::receiverMap.end() )
            UB_THROW( UbException(UB_EXARGS,"receiver is not registered") );
         
         receiver = result->second;
         if(!receiver) 
            UB_THROW( UbException(UB_EXARGS,"receiver is NULL") );
      }
      
      boost::mutex::scoped_lock lock(receiver->bufferMutex); 

      // wait until buffer is empty 
      while(receiver->isBufferFull) 
         receiver->bufferEmptiedCondition.wait(lock); 

      //////////////////////////////////////////////////////////////////////////
      // put data in buffer 
      //old: receiver->buffer = data;
      
      //new:
      if( receiver->buffer.size() != data.size() )
      {
         receiver->buffer.resize( data.size() );
      }
      memcpy( (char*)&receiver->buffer[0], (char*)&data[0],  data.size()*sizeof(typename Vector::value_type) ); 
      
      //////////////////////////////////////////////////////////////////////////
      //evtl wartende clients benachrichtigen
      //notify "buffer filled" waiters 
      receiver->isBufferFull = true; 
      receiver->bufferFilledCondition.notify_all(); // notify_one muesste eigentlich reichen 
   }

public:
   TbVectorReceiverRcf(const int& tag, const int& recvFromRank/*just for info section*/) 
      : TbTransmitter< Vector >(), tag(tag), recvFromRank(recvFromRank)
   {
      {
         //receiver registrieren
         boost::mutex::scoped_lock lock( TbVectorReceiverRcf< Vector >::staticReceiverMapMutex ); 

         std::pair< ReceiverMapIt, bool> result = receiverMap.insert(std::make_pair(tag, this));
         if( !result.second )
         {
            ReceiverMapIt existingReceiver = TbVectorReceiverRcf< Vector >::receiverMap.find(tag);
            TbVectorReceiverRcf< Vector >::receiverMap.erase(existingReceiver);
            TbVectorReceiverRcf< Vector >::receiverMap.insert(std::pair<int, TbVectorReceiverRcf< Vector >* >(tag, this));
         }
      }
      isBufferFull = false; 
   }

   ~TbVectorReceiverRcf()
   {
      ReceiverMapIt existingReceiver = receiverMap.find(tag);
      TbVectorReceiverRcf< Vector >::receiverMap.erase(existingReceiver);
   }

   bool isLocalTransmitter()  const { return false;                        }
   bool isRemoteTransmitter() const { return !this->isLocalTransmitter();  }

   void sendDataSize()    { UB_THROW( UbException(UB_EXARGS,"TbRcfVectorReceiver receives only") ); }
   void receiveDataSize() { this->receiveData(); } 

   void sendData()        { UB_THROW( UbException(UB_EXARGS,"TbRcfVectorReceiver receives only") ); }

   Vector& receiveData()
   {
      boost::mutex::scoped_lock lock(bufferMutex); 

      // wait until buffer is full 
      while(!isBufferFull) 
         bufferFilledCondition.wait(lock); 

      // get data from buffer 
      //std::size_t dataSize = this->buffer.size();
      //if( this->getData().size()!=dataSize ) this->getData().resize(dataSize);
      //for(std::size_t i=0; i<dataSize; ++i)
      //   this->getData()[i]=buffer[i]; 
      
      //folgende assert schlaegt bei receiveDataSize(), das receiveData() aufgerufen wird fehl...
      //daher ausdokumentiert
      //assert( this->getData().size() == this->buffer.size() );
      
      this->getData().swap(this->buffer);
      
      isBufferFull = false; 

      // notify "buffer emptied" waiters 
      bufferEmptiedCondition.notify_all(); // notify_one sollte auch hier reichen 
      return this->getData(); 
   }

   //info-section (usable for remote transmitter)
   int  getSendTbRank()    const { UB_THROW( UbException(UB_EXARGS,"getSendTbRank  sends only") ); }
   int  getSendTbTag()     const { UB_THROW( UbException(UB_EXARGS,"getSendTbTag  sends only") ); }
   int  getRecvFromRank()  const { return  this->recvFromRank; /*=der rank an den gesendet wird*/ }
   int  getRecvFromTag()   const { return  this->tag;                                             }

   std::string toString() const { return "TbVectorReceiverRcf< "+(std::string)typeid(Vector).name()+"<"+(std::string)typeid(typename Vector::value_type).name() +"> > to rank (tag)"+UbSystem::toString(recvFromRank)+"("+UbSystem::toString(tag)+")"; }

private:
   int tag;
   int recvFromRank; //just for info-section

   boost::mutex     bufferMutex; 
   boost::condition bufferFilledCondition; 
   boost::condition bufferEmptiedCondition; 
   bool isBufferFull; 

   Vector buffer; 
};

//////////////////////////////////////////////////////////////////////////
// template< typename Vector >
// std::map<int, TbVectorReceiverRcf< Vector >* > TbVectorReceiverRcf< Vector >::receiverMap;
template< typename Vector >
typename TbVectorReceiverRcf< Vector >::ReceiverMap TbVectorReceiverRcf< Vector >::receiverMap;

template< typename Vector >
boost::mutex TbVectorReceiverRcf< Vector >::staticReceiverMapMutex;


// 
// 
// 
// 
// //derzeit funzt es nur mit vector<double> dafuer gibt es weiter unten eine spezialisierung
// //Grund: man muss  fuer jeden datentyp eine RCF methode beim service registrieren
// //        und derzeit ist eben nur eine fuer vector<double> vorhanden
// template< typename Vector >
// class TbVectorSenderRcf : public TbTransmitter< Vector >
// {
// public:
//    TbVectorSenderRcf(const std::string& calcServiceName, const RcfConnection& receiveProcess, const int& tag) 
//       : TbTransmitter< Vector >()
//    {  
//       UB_THROW( UbException("TbVectorSenderRcf::TbVectorSenderRcf() - TbRcfVectorSender not implmeneted for that type ") );    
//    }
//    void sendDataSize()      { UB_THROW( UbException("TbVectorSenderRcf::sendDataSize() - not defined for that type") );         }
//    void receiveDataSize()   { UB_THROW( UbException("TbVectorSenderRcf::receiveDataSize() - not defined for that type") );      } 
//    void sendData()          { UB_THROW( UbException("TbVectorSenderRcf::sendData() - not defined for that type") );             }
//    void prepareForReceive() { UB_THROW( UbException("TbVectorSenderRcf::prepareForReceive() - TbRcfVectorSender sends only") ); }
//    Vector& receiveData()         { UB_THROW( UbException("TbVectorSenderRcf::receiveData() - TbRcfVectorSender sends only") );       } 
//    std::string toString()   { return "undefined TbVectorSenderRcf"; }
// };
// 
// template< typename Vector  >
// class TbVectorReceiverRcf : public TbTransmitter< Vector >
// {
// public:
//    TbVectorReceiverRcf(const int& tag, const int& receiveRank) : TbTransmitter< Vector >()
//    {
//       UB_THROW( UbException("TbVectorReceiverRcf::TbVectorReceiverRcf() - not defined for that type") );  
//    }
//    void sendDataSize()    { UB_THROW( UbException("TbVectorReceiverRcf::sendDataSize() - not defined for that type") );     }    
//    void receiveDataSize() { UB_THROW( UbException("TbVectorReceiverRcf::receiveDataSize() - not defined for that type") );  } 
//    void sendData()        { UB_THROW( UbException("TbVectorReceiverRcf::sendData() - TbRcfVectorReceiver receives only") ); }
//    Vector& receiveData()       { UB_THROW( UbException("TbVectorReceiverRcf::receiveData() - not defined for that type") );      }
//    std::string toString() { return "undefined TbVectorReceiverRcf"; }
// };
// 
#endif // CAB_RCF

#endif //TBTRANSMITTERRCF_H
