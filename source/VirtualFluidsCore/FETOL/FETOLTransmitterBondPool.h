//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifdef VF_FETOL

#ifndef FETOLTRANSMITTERBONDPOOL_H
#define FETOLTRANSMITTERBONDPOOL_H

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>

#include <mpi.h>

#include "fbond.h"
#include "Version.h"

#include <basics/transmitter/TbTransmitter.h>
#include <basics/container/CbVector.h>
#include <basics/container/CbVectorPool.h>

#include <boost/shared_ptr.hpp>

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//TbCbVectorBondPoolSender/Receiver
//diese verschicken immer einen VectorPool. Letztlich einen langen vector,
//der eigentlich aus vielen kleinen besteht
//jeder BondPoolVector hat einen pointer auf die startadresse in diesem vector
//die informationen werden im TbBondVectorPool verwaltet
//BondPoolVector verhaelt sich nach aussen hin mit einschraenkungen wie ein std::vector
//und kann somit bei den vector connector verwendet werden
//man kann die klassen theoretisch verallgemeinern.

template<typename T> class TbCbVectorSenderBondPool;
template<typename T> class TbCbVectorReceiverBondPool;

/*==================================================================*/
//! \brief The class sends always one VectorPool.
//! \details It is a long vector, which actually consists of many small.
//!         Each BondPoolVector has a pointer to the start address in this vector.
//!         The information is managed in TbCbBondVectorPool.
//!         BondPool vectors behaves outwardly with limitations such as a std::vector
//!         and can thus be used in the vector connector.
//! \author Kostyantyn Kucher

template<typename T>
class TbCbVectorBondPool : public CbVectorPool<T>
{
public:
   typedef boost::shared_ptr< TbCbVectorBondPool< T > > BondPoolPtr;

   //////////////////////////////////////////////////////////////////////////
   typedef std::map< int, BondPoolPtr >      BondPoolPtrMap;
   typedef typename BondPoolPtrMap::iterator BondPoolPtrMapIter;

   //da BasisKlasse templateKlasse ist MUSS man hier die typedefs nochmal wiederholen!
   typedef typename CbVector<T>::value_type value_type;
   typedef typename CbVector<T>::size_type  size_type;
   typedef std::vector< value_type >        Pool;

   typedef unsigned CbVectorKey;
   typedef std::map< CbVectorKey, CbVector< value_type >* /*ptrVector*/  > CbVectorMap;
   typedef typename CbVectorMap::iterator CbVectorMapIter;

   //////////////////////////////////////////////////////////////////////////
   friend class TbCbVectorSenderBondPool< T >; 
   friend class TbCbVectorReceiverBondPool< T >; 

protected:
   //////////////////////////////////////////////////////////////////////////
   static BondPoolPtrMap poolMap;
public:
   //////////////////////////////////////////////////////////////////////////
   //STATIC MEMBERS
   //////////////////////////////////////////////////////////////////////////
   //!createTbCbVectorBondPool:
   //! \param poolKey      : Schluessel fuer eindeutige Indizierung in Map
   //! \param bondRemoteRank: bond-rank of receiver/sender 
   //! \param bondTag       : with the bond-tag it is received/sent  

static BondPoolPtr createTbCbVectorBondPool(const int& poolKey, const int& bondRemoteRank, const int& bondTag, const size_type& startPoolSize = 20000 ) //startPoolSize*sizeof(T)/1024/1024 [MB]
   {
      if( poolMap.find(poolKey)!=poolMap.end() )
      {
         throw UbException(UB_EXARGS,"es ist bereits ein Pool mit dem key vorhanden!!!");
      }

      //pool erstellen
      BondPoolPtr bondPool(new TbCbVectorBondPool<T>(poolKey, bondRemoteRank, bondTag, startPoolSize) ); 

      //pool "speichern"
      TbCbVectorBondPool< value_type >::poolMap[poolKey] = bondPool;

      return bondPool; 
   }
   static void deleteTbCbVectorBondPool(const int& poolKey)
   {
      BondPoolPtrMapIter it = TbCbVectorBondPool< value_type >::poolMap.find(poolKey);
      if( it==poolMap.end() )
      {
         throw UbException(UB_EXARGS,"kein Pool mit dem key vorhanden");
      }
      TbCbVectorBondPool< value_type >::poolMap.erase(it);
   }
   /*==================================================================*/
   static BondPoolPtr getTbCbVectorBondPool(const int& poolKey)
   {
      BondPoolPtrMapIter it;
      if( (it=TbCbVectorBondPool< T >::poolMap.find(poolKey))!=TbCbVectorBondPool< T >::poolMap.end() ) 
      {
         return it->second;
      }
      return BondPoolPtr();
   }
   /*==================================================================*/
   static std::string getInfoString()
   {
      std::stringstream out;  
      out<<"TbCbVectorBondPool<"<< typeid( T ).name()  << ") - Info:"<<std::endl;
      for(BondPoolPtrMapIter it=poolMap.begin(); it!=poolMap.end(); ++it)
         out<<"pool with key("            <<std::setw(15)<<it->first<<") "
         <<"stores "                  <<std::setw(12)<<it->second->getNofStoredVectors() <<" vectors " 
         <<", elements to transfer = "<<std::setw(15)<<it->second->getPoolSize() 
         <<" ( "<< it->second->getPoolSize()*sizeof( T ) / ( 1024.0 * 1024.0 ) << " MB )" <<std::endl;
      return out.str();
   }
   /*==================================================================*/
   //! checks if all vectors have one to one pool-entries
   static bool consistencyCheck()
   {
      for(BondPoolPtrMapIter it=poolMap.begin(); it!=poolMap.end(); ++it)
      {
         if( !it->second-> CbVectorPool<T>::consistencyCheck() ) 
         {
            return false;         
         }
      }

      return true;
   }
   //////////////////////////////////////////////////////////////////////////
   static void eraseMap()
   {
      poolMap.clear();
   }
protected:
   //////////////////////////////////////////////////////////////////////////
TbCbVectorBondPool(const int& poolKey, const int& bondRemoteRank, const int& bondTag, const size_type& startPoolSize )
      :    CbVectorPool< value_type >( startPoolSize ) 
      , poolKey(poolKey)                           
      , nofStoredVectors(0) //=Anzahl an Vectoren im Pool, wird bei send/receiveDataOrder gesetzt
      , counterPrepareReceiveDataOrder(0)          
      , counterSendDataOrder(0)                    
      , counterReceiveDataOrder(0)                 
      , counterPrepareForReceive(0)                
      , counterReceive(0)                          
      , counterPrepareForSend(0)                   
      , counterSend(0)                             
      , bondRemoteRank(bondRemoteRank)               
      , bondTag(bondTag)                              
   {
      if     ( (std::string)typeid(value_type).name()==(std::string)typeid(double).name() ) mpiDataType = MPI_DOUBLE;
      else if( (std::string)typeid(value_type).name()==(std::string)typeid(float).name()  ) mpiDataType = MPI_FLOAT;
      else if( (std::string)typeid(value_type).name()==(std::string)typeid(int).name()    ) mpiDataType = MPI_INT;
      else throw UbException(UB_EXARGS,"no BondDataType for T"+(std::string)typeid(T).name());
   }
public:
   /*==================================================================*/
   //!returns key of Pool in BondPoolMap
   int  getPoolKey()    const { return  this->poolKey;       }
   /*==================================================================*/
   //!returns rank of process pool data will be send to/received from
   int  getRemoteRank() const { return  this->bondRemoteRank; }
   /*==================================================================*/
   //!returns tag of process pool data will be send to/received from
   int  getRemoteTag()  const { return  this->bondTag;        }

protected:
   /*==================================================================*/
   void sendDataOrder()
   {
      counterSendDataOrder++;
      if(counterSendDataOrder==this->cbVectorMap.size())
      {
         UBLOG(logDEBUG5, "TbCbVectorBondPool::sendDataOrder():start"<<" bondRemoteRank="<<bondRemoteRank<<" bondTag="<<bondTag);
         unsigned nofElements = (unsigned)this->cbVectorMap.size()*3+1;
         tmpSendOrderVec.resize(nofElements);//std::vector< unsigned > vec(nofElements);
         unsigned index = 0;
         tmpSendOrderVec[index++] = (unsigned)this->pool.size(); //= laenge des vectors
         if(this->nextCbVectorStartIndexInPool != this->pool.size())  throw UbException(UB_EXARGS,"an dieser Stelle sollten nextStartIndex und pool.size() identisch sein!!!");
         for(CbVectorMapIter it = this->cbVectorMap.begin(); it!=this->cbVectorMap.end(); ++it)
         {
            CbVectorKey vectorKey=0;
            size_type   dataSize=0, startIndexInPool=0;
            this->getCbVectorData(*it->second/*vec*/, vectorKey, startIndexInPool, dataSize );
            if(it->first != vectorKey) throw UbException(UB_EXARGS,"key mismatch!");

            tmpSendOrderVec[index++] = (unsigned)vectorKey;         //vectorKey == allocator.getAllocatorKey()
            tmpSendOrderVec[index++] = (unsigned)startIndexInPool;  //startIndex in poolVector
            tmpSendOrderVec[index++] = (unsigned)dataSize;          //dataSize
         }

         try
         {
            sendRequest = bond::sendFuture(&tmpSendOrderVec[0], (int)tmpSendOrderVec.size(), MPI_UNSIGNED, bondRemoteRank, bondTag);
         }
         catch (...)
         {
            std::cerr << "it is bond exception in sendDataOrder()" << std::endl;
            throw;
         }
                  
         counterSendDataOrder=0;

         nofStoredVectors = this->cbVectorMap.size();

         UBLOG(logDEBUG5, "TbCbVectorBondPool::sendDataOrder():end"<<" bondRemoteRank="<<bondRemoteRank<<" bondTag="<<bondTag);

#ifdef _DEBUG
         orgPoolVectorStartPointer = &this->pool[0];
#endif
      }
   }
   /*==================================================================*/
   void receiveDataOrder()
   {
      counterReceiveDataOrder++;
      if(counterReceiveDataOrder==this->cbVectorMap.size())
      {
         UBLOG(logDEBUG5, "TbCbVectorBondPool::receiveDataOrder():start"<<" bondRemoteRank="<<bondRemoteRank<<" bondTag="<<bondTag);

         unsigned nofElements = (unsigned)this->cbVectorMap.size()*3+1; //map MUSS auf beiden seiten gleich gross sein, sonst hat man ein grundsaetzliches problem ;-)

         std::vector< unsigned > tmpRecvOrderVec;
         tmpRecvOrderVec.resize(nofElements);
         
         try
         {
            bond::receiveComplete(&tmpRecvOrderVec[0], nofElements, MPI_UNSIGNED, bondRemoteRank, bondTag);
         }
         catch (...)
         {
            std::cerr << "it is bond exception in receiveDataOrder()" << std::endl;
            throw;
         }

         if(nofElements!=(unsigned)tmpRecvOrderVec.size())
            throw UbException(UB_EXARGS,"error... vec size stimmt nicht");

         unsigned index = 0;
         this->nextCbVectorStartIndexInPool = tmpRecvOrderVec[index++]; //= laenge des vectors
         this->pool.resize(this->nextCbVectorStartIndexInPool);
         CbVectorMapIter it = this->cbVectorMap.begin();
         for(/*index*/; index<nofElements; index+=3, ++it)
         {
            CbVectorKey vectorKey        = (CbVectorKey)tmpRecvOrderVec.at(index  );
            size_type   startIndexInPool = (size_type)tmpRecvOrderVec.at(index+1);
            size_type   dataSize         = (size_type)tmpRecvOrderVec.at(index+2);

            if(it==this->cbVectorMap.end() || it->first != vectorKey ) 
               throw UbException(UB_EXARGS,"entweder hat map nicht die gleiche reihenfolge oder vectorKey nicht vorhanden");

            this->setCbVectorData(*it->second/*vec*/, vectorKey, startIndexInPool, dataSize );
         }
         if(it!=this->cbVectorMap.end())
            throw UbException(UB_EXARGS,"error... in der map sind scheinbar noch weiter elemente vorhanden, die es auf der send seite nicht gibt...");

         counterReceiveDataOrder = 0;
         nofStoredVectors = this->cbVectorMap.size();

         UBLOG(logDEBUG5, "TbCbVectorBondPool::receiveDataOrder():end"<<" bondRemoteRank="<<bondRemoteRank<<" bondTag="<<bondTag);

#ifdef _DEBUG
         orgPoolVectorStartPointer = &this->pool[0];
#endif
      }
   }
   /*==================================================================*/
   void prepareForSendData()
   {
      //da sendDataOrder einen request verwendet muss man hier immer abfragen
      if(counterPrepareForSend==0)
      {
         UBLOG(logDEBUG5, "TbCbVectorBondPool::prepareForSendData():start"<<" bondRemoteRank="<<bondRemoteRank<<" bondTag="<<bondTag);
         if (sendRequest)
         {
            try
            {
               sendRequest->complete();
            }
            catch (...)
            {
               std::cerr << "it is bond exception in prepareForSendData()" << std::endl;
               throw;
            }
            sendRequest = std::tr1::shared_ptr<bond::FutureSend>();
         }
         UBLOG(logDEBUG5, "TbCbVectorBondPool::prepareForSendData():end"<<" bondRemoteRank="<<bondRemoteRank<<" bondTag="<<bondTag);
      }

      counterPrepareForSend++;

      if(counterPrepareForSend==nofStoredVectors)
      {
         counterPrepareForSend=0;  
      }
   }
   /*==================================================================*/
   void sendData()
   {
      counterSend++;
      if(counterSend==nofStoredVectors)
      {
#ifdef _DEBUG
         if(this->orgPoolVectorStartPointer != &this->pool[0] ) throw UbException(UB_EXARGS, "ups, pool array address changed - unknown behavior");
#endif
         UBLOG(logDEBUG5, "TbCbVectorBondPool::sendData():start"<<" bondRemoteRank="<<bondRemoteRank<<" bondTag="<<bondTag);
         try
         {
            bond::sendComplete(&this->pool[0],(int)this->nextCbVectorStartIndexInPool, mpiDataType, bondRemoteRank, bondTag);
         }
         catch(...)
         {
            std::cerr << "it is bond exception in sendData()" << std::endl;
            throw;
         }
         counterSend=0;
         UBLOG(logDEBUG5, "TbCbVectorBondPool::sendData():end"<<" bondRemoteRank="<<bondRemoteRank<<" bondTag="<<bondTag);
      }
   }
   /*==================================================================*/
   void prepareForReceiveData()
   {
      counterPrepareForReceive++;
      if(counterPrepareForReceive==this->nofStoredVectors)
      {
         UBLOG(logDEBUG5, "TbCbVectorBondPool::prepareForReceiveData():start"<<" bondRemoteRank="<<bondRemoteRank<<" bondTag="<<bondTag);
#ifdef _DEBUG
         if(this->orgPoolVectorStartPointer != &this->pool[0] ) throw UbException(UB_EXARGS, "ups, pool array adress changed - unknown behavoir");
#endif
         try
         {
            receiveRequest = bond::receiveFuture(&this->pool[0],(int)this->nextCbVectorStartIndexInPool, mpiDataType, bondRemoteRank, bondTag);
         }
         catch (...)
         {
         	std::cerr << "it is bond exception in prepareForReceiveData()" << std::endl;
            throw;
         }

         
         UBLOG(logDEBUG5, "TbCbVectorBondPool::prepareForReceiveData():end"<<" bondRemoteRank="<<bondRemoteRank<<" bondTag="<<bondTag);
         counterPrepareForReceive=0;
      }
   }
   /*==================================================================*/
   void receiveData()
   {
      if(counterReceive==0)
      {

         UBLOG(logDEBUG5, "TbCbVectorBondPool::receiveData():start"<<" bondRemoteRank="<<bondRemoteRank<<" bondTag="<<bondTag);
         
         try
         {
            receiveRequest->complete();
         }
         catch (...)
         {
            std::cerr << "it is bond exception in receiveData()" << std::endl;
            throw;
         }
         UBLOG(logDEBUG5, "TbCbVectorBondPool::receiveData():end"<<" bondRemoteRank="<<bondRemoteRank<<" bondTag="<<bondTag);
      }
      counterReceive++;
      if(counterReceive==this->nofStoredVectors) //alle receiver waren hier
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
   size_type counterPrepareForSend;
   size_type counterSend;

   std::vector< unsigned > tmpSendOrderVec; //wird zur temp speicherung der anordnung benoetigt
   MPI_Datatype mpiDataType;

   std::tr1::shared_ptr<bond::FutureSend> sendRequest;
   std::tr1::shared_ptr<bond::FutureReceive> receiveRequest;

   int bondRemoteRank, bondTag;

#ifdef _DEBUG
   T* orgPoolVectorStartPointer;
#endif
};

template<typename T>
typename TbCbVectorBondPool<T>::BondPoolPtrMap TbCbVectorBondPool<T>::poolMap;

//////////////////////////////////////////////////////////////////////////
//  TbSenderBondPool
//////////////////////////////////////////////////////////////////////////
template<typename T>
class TbCbVectorSenderBondPool : public TbTransmitter< CbVector< T >  >
{
public:
   typedef CbVector< T > value_type;   

public:
   TbCbVectorSenderBondPool(const unsigned int& cbVectorKey, TbCbVectorBondPool< T >* bondVectorPool)
      : bondVectorPool(bondVectorPool)
   { 
      this->getData().setAllocator( new CbVectorAllocatorPool<T>(cbVectorKey,this->bondVectorPool) );
   }
   ~TbCbVectorSenderBondPool()
   {
      if( this->bondVectorPool->getNofStoredVectors()==1 ) //last entry!
      {
         TbCbVectorBondPool< T >::deleteTbCbVectorBondPool(this->bondVectorPool->getPoolKey());  
      }
   }

   bool isLocalTransmitter()  const { return false;                        }
   bool isRemoteTransmitter() const { return !this->isLocalTransmitter();  }

   void sendDataSize()          { this->bondVectorPool->sendDataOrder(); }
   void receiveDataSize()       { throw UbException(UB_EXARGS,"TbBondPoolSender sends only");  }   
   CbVector< T >& receiveData() { throw UbException(UB_EXARGS,"TbBondPoolSender sends only");  }
   void prepareForSend()        { this->bondVectorPool->prepareForSendData(); }
   void sendData()              { this->bondVectorPool->sendData(); }

   //info-section (usable for remote transmitter)
   int  getSendTbRank()   const { return  this->bondVectorPool->getRemoteRank(); }
   int  getSendTbTag()    const { return  this->bondVectorPool->getRemoteTag();  }
   int  getRecvFromRank() const { throw UbException(UB_EXARGS,"TbCbVectorSenderBondPool sends only"); }
   int  getRecvFromTag()  const { throw UbException(UB_EXARGS,"TbCbVectorSenderBondPool sends only"); }

   std::string toString() const { return "TbCbVectorSenderBondPool<"+(std::string)typeid(T).name()+" to rank (tag)"+UbSystem::toString(getSendTbRank())+"("+UbSystem::toString(getSendTbTag())+")"; }

protected:
   TbCbVectorBondPool<T>* bondVectorPool;
};


/*==================================================================*/
template<typename T>
class TbCbVectorReceiverBondPool : public TbTransmitter< CbVector< T >  >
{
public:
   typedef CbVector< T > value_type;   

public:
   TbCbVectorReceiverBondPool(const unsigned int& cbVectorKey, TbCbVectorBondPool< T >* bondVectorPool)
      : bondVectorPool(bondVectorPool)
   { 
      this->getData().setAllocator( new CbVectorAllocatorPool<T>(cbVectorKey, this->bondVectorPool) );
   }
   ~TbCbVectorReceiverBondPool()
   {
      if( this->bondVectorPool->getNofStoredVectors()==1 ) //last entry!
      {
         TbCbVectorBondPool< T >::deleteTbCbVectorBondPool(this->bondVectorPool->getPoolKey());  
      }
   }
   bool isLocalTransmitter()  const { return false;                        }
   bool isRemoteTransmitter() const { return !this->isLocalTransmitter();  }

   void sendDataSize()      { throw UbException(UB_EXARGS,"TbCbVectorReceiverBondPool receives only");  }   
   void receiveDataSize()   { this->bondVectorPool->receiveDataOrder(); }  
   void sendData()          { throw UbException(UB_EXARGS,"TbCbVectorReceiverBondPool receives only"); }
   void prepareForReceive() { this->bondVectorPool->prepareForReceiveData(); }
   CbVector< T >& receiveData()
   { 
      this->bondVectorPool->receiveData();
      return this->getData();
   }

   //info-section (usable for remote transmitter)
   int  getSendTbRank()   const { throw UbException(UB_EXARGS,"TbCbVectorReceiverBondPool receives only"); }
   int  getSendTbTag()    const { throw UbException(UB_EXARGS,"TbCbVectorReceiverBondPool receives only"); }
   int  getRecvFromRank() const { return  this->bondVectorPool->getRemoteRank();  }
   int  getRecvFromTag()  const { return  this->bondVectorPool->getRemoteTag();  }

   std::string toString() const { return "TbCbVectorReceiverBondPool<"+(std::string)typeid(T).name()+" to rank (tag)"+UbSystem::toString(getRecvFromRank())+"("+UbSystem::toString(getRecvFromTag())+")"; }

protected:
   TbCbVectorBondPool<T>* bondVectorPool;
};



#endif //TBTRANSMITTERBONDPOOL_H
 
#endif //VF_FETOL
