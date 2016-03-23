//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef TBCBVECTORRCFPOOL_H
#define TBCBVECTORRCFPOOL_H

#ifdef CAB_RCF

/*=========================================================================*/
/*  ToCbVectorRcfPool                                                      */
/*                                                                         */
/**
This Class provides the base for exception handling.
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 09.10.08
*/ 

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include <cassert>

#include <RCF/Idl.hpp>
#include <RCF/TcpEndpoint.hpp>
#include <RCF/ByteBuffer.hpp>

#include <boost/shared_ptr.hpp>

#include <3rdParty/rcf/RcfSerializationIncludes.h>
#include <3rdParty/rcf/RcfSystem.h>
#include <3rdParty/rcf/IRcfIpService.h>
#include <3rdParty/rcf/RcfConnection.h>

#include <basics/transmitter/ToTransmitter.h>
#include <basics/container/CbVector.h>
#include <basics/container/CbVectorPool.h>
#include <basics/utilities/UbLogger.h>
#include <basics/utilities/UbLogger.h>


// zum ausstausch mittels RCF transmitter:
RCF_BEGIN(IRcfTransmitterReceiverPoolService, "IRcfTransmitterReceiverPoolService")
RCF_METHOD_V3(void, receiveRcfPoolData     , const float&  /*dummy*/, const int& /*poolKey*/, const RCF::ByteBuffer&         /*databuffer*/ )
RCF_METHOD_V3(void, receiveRcfPoolDataOrder, const float&  /*dummy*/, const int& /*poolKey*/, const std::vector< unsigned >& /*dataOrder*/  )
RCF_METHOD_V3(void, receiveRcfPoolData     , const double& /*dummy*/, const int& /*poolKey*/, const RCF::ByteBuffer&         /*databuffer*/ )
RCF_METHOD_V3(void, receiveRcfPoolDataOrder, const double& /*dummy*/, const int& /*poolKey*/, const std::vector< unsigned >& /*dataOrder*/  )
RCF_END(IRcfTransmitterReceiverPoolService);

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//ToRcfPoolSender/Receiver
//diese verschicken immer einen VectorPool. Letztlich einen langen vector,
//der eigentlich aus vielen kleinen besteht
//jeder RcfPoolVector hat einen pointer auf die startadresse in diesem vector
//die informationen werden im ToRcfVectorPool verwaltet
//RcfPoolVector verhaelt sich nach aussen hin mit einschraenkungen wie ein std::vector
//und kann somit bei den vector connector verwendet werden
//man kann die klassen theoretisch verallgemeinern.

template<typename T> class ToCbVectorSenderRcfPool;
template<typename T> class ToCbVectorReceiverRcfPool;

/*==================================================================*/
template<typename T>
class ToCbVectorRcfPool : public CbVectorPool<T>
{
public:
   typedef boost::shared_ptr< ToCbVectorRcfPool< T > > RcfPoolPtr;

   //////////////////////////////////////////////////////////////////////////
   typedef std::map< int, RcfPoolPtr >      RcfPoolPtrMap;
   typedef typename RcfPoolPtrMap::iterator RcfPoolPtrMapIter;

   //da BasisKlasse templateKlasse ist MUSS man hier die typedefs nochmal wiederholen!
   typedef typename CbVector<T>::value_type value_type;
   typedef typename CbVector<T>::size_type  size_type;
   typedef std::vector< value_type >        Pool;

   typedef unsigned CbVectorKey;
   typedef std::map< CbVectorKey, CbVector< value_type >* /*ptrVector*/  > CbVectorMap;
   typedef typename CbVectorMap::iterator CbVectorMapIter;

   //////////////////////////////////////////////////////////////////////////
   friend class ToCbVectorSenderRcfPool< T >; 
   friend class ToCbVectorReceiverRcfPool< T >; 

private:
   //////////////////////////////////////////////////////////////////////////
   static RcfPoolPtrMap poolMap;

//    //wenn pool als recevier fungiert
//    static ReceiverPoolMap  receiverPoolMap;
   static boost::mutex     staticPoolMapMutex;

   //wenn pool als sender fungiert
   static std::vector< RcfConnection >  recvConnections;   
   static std::vector< RcfClient<IRcfTransmitterReceiverPoolService> >  recvServices;
   static void setRcfSendToClients(const std::string& receiveServiceID, const int& rcfSendToRank);

public:
   //////////////////////////////////////////////////////////////////////////
   //STATIC MEMBERS
   //////////////////////////////////////////////////////////////////////////
   //createToCbVectorRcfPool:
   // poolKey      : Schluessel fuer eindeutige Indizierung in Map
   //              : (dieser schluessel ist eindeutig und muss bei send und recevicePool uebereinstimmen (ersetzt Tag)
   // rcfRemoteRank: rcf-rank des Empfaengers/Senders
   // rcfTag       : rcf-tag mit dem empfangen/gesendet wird
   static RcfPoolPtr createToCbVectorRcfSendPool(const std::string& rcfSendToServiceName, const int& poolKey, const int& rcfSendToRank, const size_type& startPoolSize = 20000 ) //startPoolSize*sizeof(T)/1024/1024 [MB]
   {
      if( poolMap.find(poolKey)!=ToCbVectorRcfPool< value_type >::poolMap.end() )
      {
         UB_THROW( UbException(UB_EXARGS,"es ist bereits ein Pool mit dem key vorhanden!!!") );
      }
      //pool erstellen
      RcfPoolPtr rcfPool(new ToCbVectorRcfPool<T>(rcfSendToServiceName, poolKey, rcfSendToRank, startPoolSize) ); 
      
      //pool "speichern"
      ToCbVectorRcfPool< value_type >::poolMap[poolKey] = rcfPool;

      ToCbVectorRcfPool::setRcfSendToClients(rcfSendToServiceName, rcfSendToRank);
      
      return rcfPool; 
   }
   /*==================================================================*/
   static RcfPoolPtr createToCbVectorRcfRecvPool(const int& poolKey, const int& rcfRecvRank, const size_type& startPoolSize = 20000 ) //startPoolSize*sizeof(T)/1024/1024 [MB]
   {
      if( poolMap.find(poolKey)!=poolMap.end() )
      {
         UB_THROW( UbException(UB_EXARGS,"es ist bereits ein Pool mit dem key vorhanden!!!") );
      }
      //pool erstellen
      RcfPoolPtr rcfPool(new ToCbVectorRcfPool<T>( "", poolKey, rcfRecvRank, startPoolSize ) ); 
                                                  
      //pool "speichern"
      ToCbVectorRcfPool< value_type >::poolMap[poolKey] = rcfPool;

      return rcfPool; 
   }
   /*==================================================================*/
   static void deleteToCbVectorRcfPool(const int& poolKey)
   {
      RcfPoolPtrMapIter it = ToCbVectorRcfPool< value_type >::poolMap.find(poolKey);
      if( it==poolMap.end() )
      {
         UB_THROW( UbException(UB_EXARGS,"kein Pool mit dem key vorhanden") );
      }
      ToCbVectorRcfPool< value_type >::poolMap.erase(it);
   }
   /*==================================================================*/
   static RcfPoolPtr getToCbVectorRcfPool(const int& poolKey)
   {
      RcfPoolPtrMapIter it;
      if( (it=ToCbVectorRcfPool< T >::poolMap.find(poolKey))!=ToCbVectorRcfPool< T >::poolMap.end() ) 
      {
         return it->second;
      }
      return NULL;
   }
   /*==================================================================*/
   static std::string getInfoString()
   {
      std::stringstream out;  
      out<<"ToCbVectorRcfPool<"<< typeid( T ).name()  << ") - Info:"<<std::endl;
      for(RcfPoolPtrMapIter it=poolMap.begin(); it!=poolMap.end(); ++it)
         out<<"pool with key("            <<std::setw(15)<<it->first<<") "
             <<"stores "                  <<std::setw(12)<<it->second->getNofStoredVectors() <<" vectors " 
             <<", elements to transfer = "<<std::setw(15)<<it->second->getPoolSize() 
             <<" ( "<< it->second->getPoolSize()*sizeof( T ) / ( 1024.0 * 1024.0 ) << " MB )" <<std::endl;
      return out.str();
   }
   /*==================================================================*/
   // checks if all vectors have one to one pool-entries
   static bool consistencyCheck()
   {
      for(RcfPoolPtrMapIter it=poolMap.begin(); it!=poolMap.end(); ++it)
      {
         if( !it->second-> CbVectorPool<T>::consistencyCheck() ) 
         {
            return false;         
         }
      }
      return true;
   }
   /*==================================================================*/
   static void receiveDataOrder(int poolKey, const std::vector< unsigned >& dataOrder)
   {
      RcfPoolPtr receiverPool;

      //receiver ermitteln (map nicht thread-safe->lock! aber nur kurz, ansonsten ab einer gewissen anzahl an clients/blöcken->deadlock)
      //vermutlich brauch man den nicht mal... (denn das registrieren sollte zu diesem zeitpunkt abgeschlossen sein)
      //allerdings sollte man nicht gleichzeitig Suchoperationen in ner nicht thread sicheren map durchführen!!!
      {
         boost::mutex::scoped_lock lock(staticPoolMapMutex); //wenn man das ausserhalb macht, gibt es ab einer bestimmten anzahl an clients/bloecke einen deadlock
         receiverPool = getToCbVectorRcfPool(poolKey);
         if(!receiverPool) UB_THROW( UbException(UB_EXARGS,"kein pool mit poolKey="+UbSystem::toString(poolKey)) );
      }

      boost::mutex::scoped_lock lock(receiverPool->receiveDataOrderMutex); 
      
      // wait until buffer is empty 
      while(receiverPool->receivedDataOrderBufferIsFull) 
         receiverPool->dataOrderVectorEmptiedCondition.wait(lock); 
      
      receiverPool->recvOrderVec = dataOrder;
                
      //////////////////////////////////////////////////////////////////////////
      //evtl wartende clients benachrichtigen
      //notify "buffer filled" waiters 
     
      receiverPool->receivedDataOrderBufferIsFull = true; 
      receiverPool->dataOrderVectorFilledCondition.notify_all(); 
   }
   /*==================================================================*/
   static void receivePoolData(const int& poolKey, const RCF::ByteBuffer& byteBuffer) //const typename CbVectorPool< T >::Pool& data)
   {
      RcfPoolPtr receiverPool;

      //receiver ermitteln (map nicht thread-safe->lock! aber nur kurz, ansonsten ab einer gewissen anzahl an clients/blöcken->deadlock)
      //vermutlich brauch man den nicht mal... (denn das registrieren sollte zu diesem zeitpunkt abgeschlossen sein)
      //allerdings sollte man nicht gleichzeitig Suchoperationen in ner nich thread sicheren map durchführen!!!
      {
         boost::mutex::scoped_lock lock(staticPoolMapMutex); //wenn man das ausserhalb macht, gibt es ab einer bestimmten anzahl an clients/bloecke einen deadlock
         receiverPool = getToCbVectorRcfPool(poolKey);
         if(!receiverPool) UB_THROW( UbException(UB_EXARGS,"kein pool mit poolKey="+UbSystem::toString(poolKey)) );

         //std::cout<<"poolMap.size()="<<poolMap.size()<<std::endl;
      }

      boost::mutex::scoped_lock lock(receiverPool->receiveMutex); 
      //UBLOG(logDEBUG5,"receivePoolVector - entered, pool");
      // wait until buffer is full 
      while(!receiverPool->performPoolUpdate)
         receiverPool->performPoolUpdateCond.wait(lock); 
      //UBLOG(logDEBUG5,"receivePoolVector - es kann losgehen buggercopy");

      //ACHTUNG! nie einen pool.swap(data) machen -> startadressen der "kleinen" vektoren passen sonst nimmer!
      if( receiverPool->pool.size()*sizeof( T ) != byteBuffer.getLength() )
      {
         UB_THROW( UbException(UB_EXARGS,"pool.size()!=byteBuffer.size()") );
      }
      memcpy( (char*)&receiverPool->pool[0], byteBuffer.getPtr(),  byteBuffer.getLength() ); 
//      memcpy( (char*)&receiverPool->pool[0], (char*)&data[0],  data.size()*sizeof( T ) ); 
      //UBLOG(logDEBUG5,"receivePoolVector - nach memcopy");

      receiverPool->poolWasUpdated    = true;
      receiverPool->performPoolUpdate = false;
      receiverPool->waitForPoolUpdateCond.notify_one(); // notify_one sollte auch hier reichen 
   }

protected:
   //////////////////////////////////////////////////////////////////////////
   ToCbVectorRcfPool(const std::string& rcfReceiveServiceName, const int& poolKey, const int& rcfRank, const size_type& startPoolSize  )
      :    CbVectorPool< value_type >( startPoolSize ) 
         , poolKey(poolKey)                           
         , nofStoredVectors(0) //=Anzahl an Vectoren im Pool, wird bei send/receiveDataOrder gesetzt
         , counterPrepareReceiveDataOrder(0)          
         , counterSendDataOrder(0)                    
         , counterReceiveDataOrder(0)                 
         , counterPrepareForReceive(0)                
         , counterReceive(0)                          
         , counterSend(0)                             
         , rcfReceiveServiceName(rcfReceiveServiceName)
         , rcfRank(rcfRank)                                 
         , receivedDataOrderBufferIsFull(false)
         , performPoolUpdate(false)
         , poolWasUpdated(false)
   {
   }

public:
   //returns key of Pool in RcfPoolMap
   int  getPoolKey()    const { return  this->poolKey;  }
   /*==================================================================*/
   //returns rank of process pool data will be send to/received from
   int  getRemoteRank() const { return  this->rcfRank;  }
   /*==================================================================*/
   //returns tag of process pool data will be send to/received from
   int  getRemoteTag()  const { return  this->rcfRank;  }

protected:
   /*==================================================================*/
   void sendDataOrder()
   {
      counterSendDataOrder++;
      if(counterSendDataOrder==this->cbVectorMap.size())
      {
         unsigned nofElements = (unsigned)this->cbVectorMap.size()*3+1;
         std::vector< unsigned >localSendOrderVec(nofElements); 
         unsigned index = 0;
         localSendOrderVec[index++] = (unsigned)this->pool.size(); //= laenge des vectors
         if(this->nextCbVectorStartIndexInPool != this->pool.size())  UB_THROW( UbException(UB_EXARGS,"an dieser Stelle sollten nextStartIndex und pool.size() identisch sein!!!") );
         
         for(CbVectorMapIter it = this->cbVectorMap.begin(); it!=this->cbVectorMap.end(); ++it)
         {
            CbVectorKey vectorKey=0;
            size_type   dataSize=0, startIndexInPool=0;
            this->getCbVectorData(*it->second/*vec*/, vectorKey, startIndexInPool, dataSize );
            if(it->first != vectorKey) UB_THROW( UbException(UB_EXARGS,"key mismatch!") );
            
            localSendOrderVec[index++] = (unsigned)vectorKey;         //vectorKey == allocator.getAllocatorKey()
            localSendOrderVec[index++] = (unsigned)startIndexInPool;  //startIndex in poolVector
            localSendOrderVec[index++] = (unsigned)dataSize;          //dataSize
         }
         
         //remote prozess=client erhaelt daten
         recvServices[this->rcfRank].receiveRcfPoolDataOrder(T(), this->poolKey, localSendOrderVec);
         
         counterSendDataOrder=0;

         nofStoredVectors = this->cbVectorMap.size();
      }
   }
   /*==================================================================*/
   void receiveDataOrder()
   { 
      counterReceiveDataOrder++;
      if(counterReceiveDataOrder==this->cbVectorMap.size())
      {
         boost::mutex::scoped_lock lock(receiveDataOrderMutex); 
         
         // wait until buffer is full 
         while(!receivedDataOrderBufferIsFull) 
            dataOrderVectorFilledCondition.wait(lock); //wird in receivePoolVectorForTransmitter freigegeben :)
         
         //////////////////////////////////////////////////////////////////////////
         unsigned nofElements = (unsigned)this->cbVectorMap.size()*3+1; //map MUSS auf beiden seiten gleich gross sein, sonst hat man ein grundsaetzliches problem ;-)

         if(nofElements!=(unsigned)recvOrderVec.size())
            UB_THROW( UbException(UB_EXARGS,"error... vec size stimmt nicht") );

         unsigned index = 0;
         this->nextCbVectorStartIndexInPool = recvOrderVec[index++]; //= laenge des vectors
         this->pool.resize(this->nextCbVectorStartIndexInPool);
         CbVectorMapIter it = this->cbVectorMap.begin();
         for(/*index*/; index<nofElements; index+=3, ++it)
         {
            CbVectorKey vectorKey        = (CbVectorKey)recvOrderVec.at(index  );
            size_type   startIndexInPool = (size_type)recvOrderVec.at(index+1);
            size_type   dataSize         = (size_type)recvOrderVec.at(index+2);

            if(it==this->cbVectorMap.end() || it->first != vectorKey ) 
               UB_THROW( UbException(UB_EXARGS,"entweder hat map nicht die gleiche reihenfolge oder vectorKey nicht vorhanden") );

            this->setCbVectorData(*it->second/*vec*/, vectorKey, startIndexInPool, dataSize );
         }
         if(it!=this->cbVectorMap.end())
            UB_THROW( UbException(UB_EXARGS,"error... in der map sind scheinbar noch weiter elemente vorhanden, die es auf der send seite nicht gibt...") );

         recvOrderVec.resize(0);

         // notify "buffer emptied" waiters 
         this->receivedDataOrderBufferIsFull = false; //->somit kann wieder neue reihenfolge empfangen werden
         dataOrderVectorEmptiedCondition.notify_all(); // notify_one sollte auch hier reichen 

         counterReceiveDataOrder = 0;
         nofStoredVectors = this->cbVectorMap.size();
      }
   }
   /*==================================================================*/
   void prepareForSendData() {}
   /*==================================================================*/
   void sendData()
   {
      counterSend++;
      if( counterSend == nofStoredVectors )
      {
         //remote prozess=client erhaelt daten
         //T() -> auf der empfangsseite wird automatisch die methode für den entsprechenden ToCbVectorRcfPool< T >
         //aufgerufen
         RCF::ByteBuffer byteBuffer( (char*)&this->pool[0], this->pool.size()*sizeof( T ), true );
         recvServices[this->rcfRank].receiveRcfPoolData( T(), this->poolKey, byteBuffer );
         counterSend=0;
      }
   }                              
   /*==================================================================*/
   void prepareForReceiveData() 
   {
      counterPrepareForReceive++;
      if( counterPrepareForReceive == this->nofStoredVectors )
      {
         boost::mutex::scoped_lock lock(receiveMutex); 
         //UBLOG(logDEBUG5,"prepareForReceiveData - entered -> notfifiziere performPoolUpdateCond");

         counterPrepareForReceive = 0;
         this->performPoolUpdate = true;
         performPoolUpdateCond.notify_one();
      }
   }
   /*==================================================================*/
   void receiveData()
   {
      if( counterReceive == 0 )
      {
         boost::mutex::scoped_lock lock(receiveMutex); 
         //UBLOG(logDEBUG5,"receiveData - wait for pool update");

         while(!this->poolWasUpdated)
            waitForPoolUpdateCond.wait(lock);
         this->poolWasUpdated    = false;
         //UBLOG(logDEBUG5,"receiveData - pool update seems to be finished");
      }

      counterReceive++;
      if( counterReceive == this->nofStoredVectors ) //alle receiver waren hier
      {
         counterReceive=0;
      }
   }

protected:
   int       poolKey; //eindeutiger schluessel fuer pool
   size_type nofStoredVectors;

   size_type counterPrepareReceiveDataOrder;
   size_type counterSendDataOrder;
   size_type counterReceiveDataOrder;
   size_type counterPrepareForReceive;
   size_type counterReceive;
   size_type counterSend;

   std::string rcfReceiveServiceName; //nur als SENDER wichtig!!
   int         rcfRank; 

   bool             receivedDataOrderBufferIsFull; 
   boost::mutex     receiveDataOrderMutex; 
   boost::condition dataOrderVectorFilledCondition; 
   boost::condition dataOrderVectorEmptiedCondition;
   std::vector< unsigned > recvOrderVec;

   bool             performPoolUpdate;
   boost::mutex     receiveMutex;
   bool             poolWasUpdated;
   boost::condition waitForPoolUpdateCond;
   boost::condition performPoolUpdateCond;
};

//////////////////////////////////////////////////////////////////////////
template< typename T >
std::vector< RcfClient<IRcfTransmitterReceiverPoolService> >  ToCbVectorRcfPool< T >::recvServices;

template< typename T >
std::vector< RcfConnection >  ToCbVectorRcfPool< T >::recvConnections;

template< typename T >                                              
void ToCbVectorRcfPool< T >::setRcfSendToClients(const std::string& recvServiceID, const int& rcfSendToRank)
{
   UBLOG(logINFO,"ToCbVectorRcfPool< T >::setRcfSendToClients - invoked setRcfClients");
   RcfConnection ipConnection = RcfSystem::getIpServiceConnection();
   if(!ipConnection) UB_THROW( UbException(UB_EXARGS,"ups, no IpServiceConnection") );
   RcfClient<IRcfIpService> ipService(RCF::TcpEndpoint(ipConnection.getIp(), ipConnection.getPort()));

   //////////////////////////////////////////////////////////////////////////
   //CalcService Verbindungsdaten holen und nach rank sortiere
   std::vector<RcfConnection> connections = ipService.getServiceConnections(recvServiceID);
   if(connections.empty()) UB_THROW( UbException(UB_EXARGS,"no existing RecvService with ID = "+recvServiceID) );
   std::sort(connections.begin(),connections.end(),RcfConnection::compareRank());

   //////////////////////////////////////////////////////////////////////////
   //CalcServiceClient für rcfSendToRank übernehmen
   assert( recvConnections.size() == recvServices.size() );

   if( rcfSendToRank >= (int)recvConnections.size() ) 
   {
      recvConnections.resize(rcfSendToRank+1);
      recvServices.resize( rcfSendToRank+1 );
   }
   
   //Anm.: nur, wenn nicht schon vorhanden (hierfür merkt man sich zusätzlich die connection!
   if( recvConnections[rcfSendToRank] != connections[rcfSendToRank] )
   {
      if( connections[rcfSendToRank].getRank() != rcfSendToRank )
         UB_THROW( UbException(UB_EXARGS,"error - ranks ranks anscheinend nicht kontinierlich von [0..n]") );
      
      recvConnections[rcfSendToRank] = connections[rcfSendToRank];
      recvServices[rcfSendToRank] = RcfClient<IRcfTransmitterReceiverPoolService>(RCF::TcpEndpoint(connections[rcfSendToRank].getIp(), connections[rcfSendToRank].getPort()));
      UBLOG(logINFO,"ToCbVectorRcfPool< T >::setRcfSendToClients - rank="<<rcfSendToRank<<" : set RcfClient with connection = " << connections[rcfSendToRank] );
   }
   else
   {
       UBLOG(logINFO,"ToCbVectorRcfPool< T >::setRcfSendToClients - rank="<<rcfSendToRank<<" : RcfClient already exists with connection = " << connections[rcfSendToRank] );
   }
   //for(std::size_t i=0; i<connections.size(); i++)
   //{
   //   if( (int)i != connections[i].getRank() )
   //      UB_THROW( UbException(UB_EXARGS,"recvServices must have continous ranks sarting from 0") );
   //   recvServices[i] = RcfClient<IRcfTransmitterReceiverPoolService>(RCF::TcpEndpoint(connections[i].getIp(), connections[i].getPort()));
   //   UBLOG(logINFO,"ToCbVectorRcfPool< T >::setRcfSendToClients - pos="<<i<<" : set RcfClient with connection = " << connections[i] );
   //}
}

template<typename T>
typename ToCbVectorRcfPool<T>::RcfPoolPtrMap ToCbVectorRcfPool<T>::poolMap;

template< typename T >
boost::mutex ToCbVectorRcfPool< T >::staticPoolMapMutex;


//////////////////////////////////////////////////////////////////////////
//  ToSenderRcfPool
//////////////////////////////////////////////////////////////////////////
template<typename T>
class ToCbVectorSenderRcfPool : public ToTransmitter< CbVector< T >  >
{
public:
   typedef CbVector< T > value_type;   

public:
   ToCbVectorSenderRcfPool(const unsigned int& cbVectorKey, ToCbVectorRcfPool< T >* rcfVectorPool)
      : rcfVectorPool(rcfVectorPool)
   { 
      this->getData().setAllocator( new CbVectorAllocatorPool<T>(cbVectorKey,this->rcfVectorPool) );
   }
   ~ToCbVectorSenderRcfPool()
   {
      if( this->rcfVectorPool->getNofStoredVectors()==1 ) //last entry!
      {
         ToCbVectorRcfPool< T >::deleteToCbVectorRcfPool(this->rcfVectorPool->getPoolKey());  
      }
   }

   bool isLocalTransmitter()  const { return false;                        }
   bool isRemoteTransmitter() const { return !this->isLocalTransmitter();  }

   void sendDataSize()          { this->rcfVectorPool->sendDataOrder(); }
   void receiveDataSize()       { UB_THROW( UbException(UB_EXARGS,"ToRcfPoolSender sends only") );  }   
   CbVector< T >& receiveData() { UB_THROW( UbException(UB_EXARGS,"ToRcfPoolSender sends only") );  }
   void prepareForSend()        { this->rcfVectorPool->prepareForSendData(); }
   void sendData()              { this->rcfVectorPool->sendData(); }

   //info-section (usable for remote transmitter)
   int  getSendToRank()   const { return  this->rcfVectorPool->getRemoteRank(); }
   int  getSendToTag()    const { return  this->rcfVectorPool->getRemoteTag();  }
   int  getRecvFromRank() const { UB_THROW( UbException(UB_EXARGS,"ToCbVectorSenderRcfPool sends only") ); }
   int  getRecvFromTag()  const { UB_THROW( UbException(UB_EXARGS,"ToCbVectorSenderRcfPool sends only") ); }

   std::string toString() const { return "ToCbVectorSenderRcfPool<"+(std::string)typeid(T).name()+" to rank (tag)"+UbSystem::toString(getSendToRank())+"("+UbSystem::toString(getSendToTag())+")"; }

protected:
   ToCbVectorRcfPool<T>* rcfVectorPool;
};

/*==================================================================*/
template<typename T>
class ToCbVectorReceiverRcfPool : public ToTransmitter< CbVector< T >  >
{
public:
   typedef CbVector< T > value_type;   

public:
   ToCbVectorReceiverRcfPool(const unsigned int& cbVectorKey, ToCbVectorRcfPool< T >* rcfVectorPool)
      : rcfVectorPool(rcfVectorPool)
   { 
      this->getData().setAllocator( new CbVectorAllocatorPool<T>(cbVectorKey, this->rcfVectorPool) );
   }
   ~ToCbVectorReceiverRcfPool()
   {
      if( this->rcfVectorPool->getNofStoredVectors()==1 ) //last entry!
      {
         UBLOG(logINFO,"ToCbVectorReceiverRcfPool - loesche map - poolKey "<<this->rcfVectorPool->getPoolKey());
         ToCbVectorRcfPool< T >::deleteToCbVectorRcfPool(this->rcfVectorPool->getPoolKey());  
      }
   }
   bool isLocalTransmitter()  const { return false;                         }
   bool isRemoteTransmitter() const { return !this->isLocalTransmitter();   }

   void sendDataSize()      { UB_THROW( UbException(UB_EXARGS,"ToCbVectorReceiverRcfPool receives only") ); }   
   void receiveDataSize()   { this->rcfVectorPool->receiveDataOrder();   }  
   void sendData()          { UB_THROW( UbException(UB_EXARGS,"ToCbVectorReceiverRcfPool receives only") ); }
   void prepareForReceive() { this->rcfVectorPool->prepareForReceiveData(); }
   CbVector< T >& receiveData()   { this->rcfVectorPool->receiveData(); return this->getData();  }

   //info-section (usable for remote transmitter)
   int  getSendToRank()   const { UB_THROW( UbException(UB_EXARGS,"ToCbVectorReceiverRcfPool receives only") ); }
   int  getSendToTag()    const { UB_THROW( UbException(UB_EXARGS,"ToCbVectorReceiverRcfPool receives only") ); }
   int  getRecvFromRank() const { return  this->rcfVectorPool->getRemoteRank();  }
   int  getRecvFromTag()  const { return  this->rcfVectorPool->getRemoteTag();  }

   std::string toString() const { return "ToCbVectorReceiverRcfPool<"+(std::string)typeid(T).name()+" to rank (tag)"+UbSystem::toString(getRecvFromRank())+"("+UbSystem::toString(getRecvFromTag())+")"; }

protected:
   ToCbVectorRcfPool<T>* rcfVectorPool;
};

#endif //CAB_RCF

#endif //TOCBVECTORRCFPOOL_H
