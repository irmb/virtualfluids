//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef TBTRANSMITTERMPIPOOL_H
#define TBTRANSMITTERMPIPOOL_H

#ifdef VF_MPI

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>

#include <mpi.h>

#include <basics/transmitter/TbTransmitter.h>
#include <basics/container/CbVector.h>
#include <basics/container/CbVectorPool.h>

#include <PointerDefinitions.h>
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//TbCbVectorMpiPoolSender/Receiver
//diese verschicken immer einen VectorPool. Letztlich einen langen vector,
//der eigentlich aus vielen kleinen besteht
//jeder MpiPoolVector hat einen pointer auf die startadresse in diesem vector
//die informationen werden im TbMpiVectorPool verwaltet
//MpiPoolVector verhaelt sich nach aussen hin mit einschraenkungen wie ein std::vector
//und kann somit bei den vector connector verwendet werden
//man kann die klassen theoretisch verallgemeinern.

template<typename T> class TbCbVectorSenderMpiPool;
template<typename T> class TbCbVectorReceiverMpiPool;

/*==================================================================*/
template<typename T>
class TbCbVectorMpiPool : public CbVectorPool<T>
{
public:
   typedef SPtr< TbCbVectorMpiPool< T > > MpiPoolPtr;

   //////////////////////////////////////////////////////////////////////////
   typedef std::map<std::string, MpiPoolPtr >      MpiPoolPtrMap;
   typedef typename MpiPoolPtrMap::iterator MpiPoolPtrMapIter;

   //da BasisKlasse templateKlasse ist MUSS man hier die typedefs nochmal wiederholen!
   typedef typename CbVector<T>::value_type value_type;
   typedef typename CbVector<T>::size_type  size_type;
   typedef std::vector< value_type >        Pool;

   typedef std::string CbVectorKey;
   typedef std::map< CbVectorKey, CbVector< value_type >* /*ptrVector*/  > CbVectorMap;
   typedef typename CbVectorMap::iterator CbVectorMapIter;

   //////////////////////////////////////////////////////////////////////////
   friend class TbCbVectorSenderMpiPool< T >; 
   friend class TbCbVectorReceiverMpiPool< T >; 

protected:
   //////////////////////////////////////////////////////////////////////////
   static MpiPoolPtrMap poolMap;
public:
   //////////////////////////////////////////////////////////////////////////
   //STATIC MEMBERS
   //////////////////////////////////////////////////////////////////////////
   //createTbCbVectorMpiPool:
   // poolKey      : Schluessel fuer eindeutige Indizierung in Map
   // mpiRemoteRank: mpi-rank des Empfaengers/Senders
   // mpiTag       : mpi-tag mit dem empfangen/gesendet wird
   static MpiPoolPtr createTbCbVectorMpiPool(CbVectorKey poolKey, int mpiRemoteRank, int mpiTag, MPI_Comm comm, size_type startPoolSize = 20000 ) //startPoolSize*sizeof(T)/1024/1024 [MB]

   {
      if( poolMap.find(poolKey)!=poolMap.end() )
      {
         throw UbException(UB_EXARGS,"es ist bereits ein Pool mit dem key vorhanden!!!");
      }

      //pool erstellen
      MpiPoolPtr mpiPool(new TbCbVectorMpiPool<T>(poolKey, mpiRemoteRank, mpiTag, comm, startPoolSize) ); 

      //pool "speichern"
      TbCbVectorMpiPool< value_type >::poolMap[poolKey] = mpiPool;

      return mpiPool; 
   }
   static void deleteTbCbVectorMpiPool(CbVectorKey poolKey)
   {
      MpiPoolPtrMapIter it = TbCbVectorMpiPool< value_type >::poolMap.find(poolKey);
      if( it==poolMap.end() )
      {
         throw UbException(UB_EXARGS,"kein Pool mit dem key vorhanden");
      }
      TbCbVectorMpiPool< value_type >::poolMap.erase(it);
   }
   /*==================================================================*/
   static MpiPoolPtr getTbCbVectorMpiPool(CbVectorKey poolKey)
   {
      MpiPoolPtrMapIter it;
      if( (it=TbCbVectorMpiPool< T >::poolMap.find(poolKey))!=TbCbVectorMpiPool< T >::poolMap.end() ) 
      {
         return it->second;
      }
      return MpiPoolPtr();
   }
   /*==================================================================*/
   static std::string getInfoString()
   {
      std::stringstream out;  
      out<<"TbCbVectorMpiPool<"<< typeid( T ).name()  << ") - Info:"<<std::endl;
      for(MpiPoolPtrMapIter it=poolMap.begin(); it!=poolMap.end(); ++it)
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
      for(MpiPoolPtrMapIter it=poolMap.begin(); it!=poolMap.end(); ++it)
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
   TbCbVectorMpiPool(CbVectorKey poolKey, int mpiRemoteRank, int mpiTag, MPI_Comm comm, size_type startPoolSize )
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
      , comm(comm)                                 
      , receiveRequest(MPI_REQUEST_NULL)
      //, sendRequest(MPI_REQUEST_NULL)
      , mpiRemoteRank(mpiRemoteRank)               
      , mpiTag(mpiTag)                              
   {
      if     ( (std::string)typeid(value_type).name()==(std::string)typeid(double).name() ) mpiDataType = MPI_DOUBLE;
      else if( (std::string)typeid(value_type).name()==(std::string)typeid(float).name()  ) mpiDataType = MPI_FLOAT;
      else if( (std::string)typeid(value_type).name()==(std::string)typeid(int).name()    ) mpiDataType = MPI_INT;
      else throw UbException(UB_EXARGS,"no MpiDataType for T"+(std::string)typeid(T).name());

      for (int i = 0; i < 3; i++)
      {
         sendRequest[i] = MPI_REQUEST_NULL;
      }
   }

public:
   /*==================================================================*/
   //returns key of Pool in MpiPoolMap
   CbVectorKey  getPoolKey()    const { return  this->poolKey;       }
   /*==================================================================*/
   //returns rank of process pool data will be send to/received from
   int  getRemoteRank() const { return  this->mpiRemoteRank; }
   /*==================================================================*/
   //returns tag of process pool data will be send to/received from
   int  getRemoteTag()  const { return  this->mpiTag;        }

protected:
   /*==================================================================*/
   void sendDataOrder()
   {
      counterSendDataOrder++;
      if(counterSendDataOrder==this->cbVectorMap.size())
      {
         //allg.: bei MPI muss man darauf achten, dass bei unblocked operationen die puffer (aus dem oder in den 
         //geschrieben wird auch noch vorhanden sind!!! wuerde man hier z.B. einen lokalen vector mit Isend() los-
         //schicken, dann wurde der scope verlassen werden und der vector evtl geloescht werden, bevor mpi den
         //vorgang abgeschlossen hat!!! ->  tmpOrderVec ist class-member!!!
         unsigned nofElements = (unsigned)this->cbVectorMap.size()*3+1;
         tmpSendOrderVec.resize(nofElements);//std::vector< unsigned > vec(nofElements);
         unsigned index = 0;
         tmpSendOrderVec[index++] = (unsigned)this->pool.size(); //= laenge des vectors
         if(this->nextCbVectorStartIndexInPool != this->pool.size())  throw UbException(UB_EXARGS,"an dieser Stelle sollten nextStartIndex und pool.size() identisch sein!!!");
         for(CbVectorMapIter it = this->cbVectorMap.begin(); it!=this->cbVectorMap.end(); ++it)
         {
            CbVectorKey vectorKey;
            size_type   dataSize=0, startIndexInPool=0;
            this->getCbVectorData(*it->second/*vec*/, vectorKey, startIndexInPool, dataSize );
            if(it->first != vectorKey) throw UbException(UB_EXARGS,"key mismatch!");

            tmpSendOrderKeyVec += vectorKey;         //vectorKey == allocator.getAllocatorKey()
            tmpSendOrderVec[index++] = (unsigned)vectorKey.length();
            tmpSendOrderVec[index++] = (unsigned)startIndexInPool;  //startIndex in poolVector
            tmpSendOrderVec[index++] = (unsigned)dataSize;          //dataSize
         }
         
         MPI_Isend(&tmpSendOrderVec[0],(int)tmpSendOrderVec.size(), MPI_UNSIGNED, mpiRemoteRank, mpiTag, comm, &sendRequest[0]);
         
         tmpSendOrderKeyVecLength = (unsigned)tmpSendOrderKeyVec.length();
         MPI_Isend(&tmpSendOrderKeyVecLength,1, MPI_UNSIGNED, mpiRemoteRank, mpiTag, comm, &sendRequest[1]);
         MPI_Isend((char *)tmpSendOrderKeyVec.c_str(),tmpSendOrderKeyVecLength, MPI_CHAR, mpiRemoteRank, mpiTag, comm, &sendRequest[2]);


         counterSendDataOrder=0;

         nofStoredVectors = this->cbVectorMap.size();

         UBLOG(logDEBUG5, "TbCbVectorMpiPool::sendDataOrder()"<<" mpiRemoteRank="<<mpiRemoteRank<<" mpiTag="<<mpiTag);

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
         //receiveRequest.Wait();
         unsigned nofElements = (unsigned)this->cbVectorMap.size()*3+1; //map MUSS auf beiden seiten gleich gross sein, sonst hat man ein grundsaetzliches problem ;-)

         UBLOG(logDEBUG5, "TbCbVectorMpiPool::receiveDataOrder()"<<" mpiRemoteRank="<<mpiRemoteRank<<" mpiTag="<<mpiTag);

         std::vector< unsigned > tmpRecvOrderVec;
         tmpRecvOrderVec.resize(nofElements);

         std::vector<char> tmpRecvOrderKeyVec;
        
         //MPI_Status status;
         MPI_Recv(&tmpRecvOrderVec[0], nofElements, MPI_UNSIGNED, mpiRemoteRank, mpiTag, comm, MPI_STATUS_IGNORE);

         unsigned rcount;
         MPI_Recv(&rcount, 1, MPI_UNSIGNED, mpiRemoteRank, mpiTag, comm, MPI_STATUS_IGNORE);
         tmpRecvOrderKeyVec.resize(rcount);
         MPI_Recv(&tmpRecvOrderKeyVec[0], rcount, MPI_CHAR, mpiRemoteRank, mpiTag, comm, MPI_STATUS_IGNORE);

         if(nofElements!=(unsigned)tmpRecvOrderVec.size())
            throw UbException(UB_EXARGS,"error... vec size stimmt nicht");

         unsigned index = 0;
         size_type index2 = 0;
         this->nextCbVectorStartIndexInPool = tmpRecvOrderVec[index++]; //= laenge des vectors
         this->pool.resize(this->nextCbVectorStartIndexInPool);
         CbVectorMapIter it = this->cbVectorMap.begin();
         for(/*index*/; index<nofElements; index+=3, ++it)
         {
            size_type   vectorKeyLength  = (size_type)tmpRecvOrderVec.at(index  );
            size_type   startIndexInPool = (size_type)tmpRecvOrderVec.at(index+1);
            size_type   dataSize         = (size_type)tmpRecvOrderVec.at(index+2);
            CbVectorKey vectorKey = CbVectorKey(&tmpRecvOrderKeyVec[index2], vectorKeyLength);
            index2 += vectorKeyLength;

            //if(it==this->cbVectorMap.end() || it->first != vectorKey ) 
               //throw UbException(UB_EXARGS, "entweder hat map nicht die gleiche reihenfolge oder vectorKey = "+UbSystem::toString(vectorKey)+" nicht vorhanden");
            if (it==this->cbVectorMap.end())
               throw UbException(UB_EXARGS,"map ist leer");
            else if (it->first != vectorKey)
               throw UbException(UB_EXARGS, "vectorKey = "+UbSystem::toString(vectorKey)+" nicht vorhanden it->first =" + UbSystem::toString(it->first));

            this->setCbVectorData(*it->second/*vec*/, vectorKey, startIndexInPool, dataSize );
            
         }
         if(it!=this->cbVectorMap.end())
            throw UbException(UB_EXARGS,"error... in der map sind scheinbar noch weiter elemente vorhanden, die es auf der send seite nicht gibt...");

         counterReceiveDataOrder = 0;
         nofStoredVectors = this->cbVectorMap.size();

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
         UBLOG(logDEBUG5, "TbCbVectorMpiPool::prepareForSendData():start"<<" mpiRemoteRank="<<mpiRemoteRank<<" mpiTag="<<mpiTag);
         //if(sendRequest != MPI_REQUEST_NULL) MPI_Wait(&sendRequest, MPI_STATUS_IGNORE);
         if(sendRequest[2] != MPI_REQUEST_NULL) MPI_Waitall(3, sendRequest, MPI_STATUS_IGNORE);
         UBLOG(logDEBUG5, "TbCbVectorMpiPool::prepareForSendData():end"<<" mpiRemoteRank="<<mpiRemoteRank<<" mpiTag="<<mpiTag);
      }

      counterPrepareForSend++;

      if(counterPrepareForSend==nofStoredVectors)
      {
         counterPrepareForSend=0;  
      }


      //A - non blocking
      ////der ERSTE is entscheidend 
      ////Grund: wenn man 
      //// for(all trans) { trans->prepare(); trans->fillBuffer(); }
      //// aufruft, dann wuerde u.U. der Buffer neu beschrieben werden obwohl noch nicht versendet wurde!!!
      //counterPrepareForSend++;
      //if(counterPrepareForSend==1)
      //{
      //   if(sendRequest != MPI::REQUEST_NULL) sendRequest.Wait();
      //}
      //
      //if(counterPrepareForSend==nofStoredVectors)
      //   counterPrepareForSend=0;  
   }
   /*==================================================================*/
   void sendData()
   {
      //A - non blocking
      //der LETZTE is entscheidend 
      //counterSend++;
      //if(counterSend==nofStoredVectors)
      //{
      //   //std::cout<<"Isend von "<<(int)nextStartIndex<<"elementen von "<<mpiRemoteRank<<" mit tag="<<mpiTag<<std::endl;
      //   sendRequest = comm.Isend(&pool[0],(int)nextCbVectorStartIndexInPool, mpiDataType, mpiRemoteRank, mpiTag);
      //   counterSend=0;
      //}
      //B - blocking
      //der LETZTE is entscheidend 
      counterSend++;
      if(counterSend==nofStoredVectors)
      {
         //std::cout<<"Isend von "<<(int)nextStartIndex<<"elementen von "<<mpiRemoteRank<<" mit tag="<<mpiTag<<std::endl;
         UBLOG(logDEBUG5, "TbCbVectorMpiPool::sendData():start"<<" mpiRemoteRank="<<mpiRemoteRank<<" mpiTag="<<mpiTag);

         //synchronous send 
         //comm.Ssend(&this->pool[0],(int)this->nextCbVectorStartIndexInPool, mpiDataType, mpiRemoteRank, mpiTag);
#ifdef _DEBUG
         if(this->orgPoolVectorStartPointer != &this->pool[0] ) throw UbException(UB_EXARGS, "ups, pool array adress changed - unknown behavoir");
#endif

         //standard send
         MPI_Send(&this->pool[0],(int)this->nextCbVectorStartIndexInPool, mpiDataType, mpiRemoteRank, mpiTag, comm);
         UBLOG(logDEBUG5, "TbCbVectorMpiPool::sendData():end"<<" mpiRemoteRank="<<mpiRemoteRank<<" mpiTag="<<mpiTag);
////////////////////////////////////////////////////////////////////////////////////////////
//DEBUG///////////////////////////////////////
         //int irank;
         //MPI_Comm_rank(MPI_COMM_WORLD, &irank);
         //std::cout << "MPI_Send: " << irank <<  " "  << mpiRemoteRank << " "  <<mpiTag<<std::endl;
///////////////////////////////////////////////////
         counterSend=0;
      }                           
   }
   /*==================================================================*/
   void prepareForReceiveData()
   {
      //A - non blocking
      //sobald der Letzte kann man den den request holen.
      //andernfalls kann nicht gewaehrleistet werden, dass evtl noch mit dem buffer gearbeitet wird!!!
      counterPrepareForReceive++;
      if(counterPrepareForReceive==this->nofStoredVectors)
      {
         UBLOG(logDEBUG5, "TbCbVectorMpiPool::prepareForReceiveData():start"<<" mpiRemoteRank="<<mpiRemoteRank<<" mpiTag="<<mpiTag);
#ifdef _DEBUG
         if(this->orgPoolVectorStartPointer != &this->pool[0] ) throw UbException(UB_EXARGS, "ups, pool array adress changed - unknown behavoir");
#endif
         MPI_Irecv(&this->pool[0],(int)this->nextCbVectorStartIndexInPool, mpiDataType, mpiRemoteRank, mpiTag, comm, &receiveRequest);
         UBLOG(logDEBUG5, "TbCbVectorMpiPool::prepareForReceiveData():end"<<" mpiRemoteRank="<<mpiRemoteRank<<" mpiTag="<<mpiTag);
         counterPrepareForReceive=0;
      }
   }
   /*==================================================================*/
   void receiveData()
   {
      //A - non blocking
      //sobald der ERSTE reinkommt muss man warten, bis received wurde!!!
      //denn erst anschliessend stehen die empfangenen daten zur verfuegung
      if(counterReceive==0)
      {
         UBLOG(logDEBUG5, "TbCbVectorMpiPool::receiveData():start"<<" mpiRemoteRank="<<mpiRemoteRank<<" mpiTag="<<mpiTag);
         MPI_Wait(&receiveRequest, MPI_STATUS_IGNORE);
         UBLOG(logDEBUG5, "TbCbVectorMpiPool::receiveData():end"<<" mpiRemoteRank="<<mpiRemoteRank<<" mpiTag="<<mpiTag);
      }
      counterReceive++;
      if(counterReceive==this->nofStoredVectors) //alle receiver waren hier
      {
         counterReceive=0;
      }

      ////B - blocking
      ////sobald der ERSTE reinkommt muss man warten, bis received wurde!!!
      ////denn erst anschliessend stehen die empfangenen daten zur verfuegung
      //if(counterReceive==0)
      //{
      //   comm.Recv(&this->pool[0],(int)this->nextCbVectorStartIndexInPool, mpiDataType, mpiRemoteRank, mpiTag);
      //}
      //counterReceive++;
      //if(counterReceive==this->nofStoredVectors) //alle receiver waren hier
      //   counterReceive=0;
   }

protected:
   CbVectorKey poolKey; //eindeutiger schluessel fuer pool
   size_type nofStoredVectors;

   size_type counterPrepareReceiveDataOrder;
   size_type counterSendDataOrder;
   size_type counterReceiveDataOrder;
   size_type counterPrepareForReceive;
   size_type counterReceive;
   size_type counterPrepareForSend;
   size_type counterSend;

   std::vector< unsigned > tmpSendOrderVec; //wird zur temp speicherung der anordnung benoetigt
   std::string tmpSendOrderKeyVec;
   unsigned tmpSendOrderKeyVecLength;

   MPI_Comm     comm;
   MPI_Request  receiveRequest;
   //MPI_Request  sendRequest;
   MPI_Request  sendRequest[3];
   //MPI_Status   sendStatus;
   //MPI_Status   receiveStatus;
   MPI_Datatype mpiDataType;

   int mpiRemoteRank, mpiTag;

#ifdef _DEBUG
   T* orgPoolVectorStartPointer;
#endif
};

template<typename T>
typename TbCbVectorMpiPool<T>::MpiPoolPtrMap TbCbVectorMpiPool<T>::poolMap;

//////////////////////////////////////////////////////////////////////////
//  TbSenderMpiPool
//////////////////////////////////////////////////////////////////////////
template<typename T>
class TbCbVectorSenderMpiPool : public TbTransmitter< CbVector< T >  >
{
public:
   typedef CbVector< T > value_type;

public:
   TbCbVectorSenderMpiPool(std::string cbVectorKey, TbCbVectorMpiPool< T >* mpiVectorPool)
      : mpiVectorPool(mpiVectorPool)
   { 
      this->getData().setAllocator( new CbVectorAllocatorPool<T>(cbVectorKey,this->mpiVectorPool) );
   }
   ~TbCbVectorSenderMpiPool()
   {
      if( this->mpiVectorPool->getNofStoredVectors()==1 ) //last entry!
      {
         TbCbVectorMpiPool< T >::deleteTbCbVectorMpiPool(this->mpiVectorPool->getPoolKey());  
      }
   }

   bool isLocalTransmitter()  const { return false;                        }
   bool isRemoteTransmitter() const { return !this->isLocalTransmitter();  }

   void sendDataSize()          { this->mpiVectorPool->sendDataOrder(); }
   void receiveDataSize()       { throw UbException(UB_EXARGS,"TbMpiPoolSender sends only");  }   
   CbVector< T >& receiveData() { throw UbException(UB_EXARGS,"TbMpiPoolSender sends only");  }
   void prepareForSend()        { this->mpiVectorPool->prepareForSendData(); }
   void sendData()              { this->mpiVectorPool->sendData(); }

   //info-section (usable for remote transmitter)
   int  getSendTbRank()   const { return  this->mpiVectorPool->getRemoteRank(); }
   int  getSendTbTag()    const { return  this->mpiVectorPool->getRemoteTag();  }
   int  getRecvFromRank() const { throw UbException(UB_EXARGS,"TbCbVectorSenderMpiPool sends only"); }
   int  getRecvFromTag()  const { throw UbException(UB_EXARGS,"TbCbVectorSenderMpiPool sends only"); }

   std::string toString() const { return "TbCbVectorSenderMpiPool<"+(std::string)typeid(T).name()+" to rank (tag)"+UbSystem::toString(getSendTbRank())+"("+UbSystem::toString(getSendTbTag())+")"; }

protected:
   TbCbVectorMpiPool<T>* mpiVectorPool;
};


/*==================================================================*/
template<typename T>
class TbCbVectorReceiverMpiPool : public TbTransmitter< CbVector< T >  >
{
public:
   typedef CbVector< T > value_type;   

public:
   TbCbVectorReceiverMpiPool(std::string cbVectorKey, TbCbVectorMpiPool< T >* mpiVectorPool)
      : mpiVectorPool(mpiVectorPool)
   { 
      this->getData().setAllocator( new CbVectorAllocatorPool<T>(cbVectorKey, this->mpiVectorPool) );
   }
   ~TbCbVectorReceiverMpiPool()
   {
      if( this->mpiVectorPool->getNofStoredVectors()==1 ) //last entry!
      {
         TbCbVectorMpiPool< T >::deleteTbCbVectorMpiPool(this->mpiVectorPool->getPoolKey());  
      }
   }
   bool isLocalTransmitter()  const { return false;                        }
   bool isRemoteTransmitter() const { return !this->isLocalTransmitter();  }

   void sendDataSize()      { throw UbException(UB_EXARGS,"TbCbVectorReceiverMpiPool receives only");  }   
   void receiveDataSize()   { this->mpiVectorPool->receiveDataOrder(); }  
   void sendData()          { throw UbException(UB_EXARGS,"TbCbVectorReceiverMpiPool receives only"); }
   void prepareForReceive() { this->mpiVectorPool->prepareForReceiveData(); }
   CbVector< T >& receiveData()
   { 
      this->mpiVectorPool->receiveData();
      return this->getData();
   }

   //info-section (usable for remote transmitter)
   int  getSendTbRank()   const { throw UbException(UB_EXARGS,"TbCbVectorReceiverMpiPool receives only"); }
   int  getSendTbTag()    const { throw UbException(UB_EXARGS,"TbCbVectorReceiverMpiPool receives only"); }
   int  getRecvFromRank() const { return  this->mpiVectorPool->getRemoteRank();  }
   int  getRecvFromTag()  const { return  this->mpiVectorPool->getRemoteTag();  }

   std::string toString() const { return "TbCbVectorReceiverMpiPool<"+(std::string)typeid(T).name()+" to rank (tag)"+UbSystem::toString(getRecvFromRank())+"("+UbSystem::toString(getRecvFromTag())+")"; }

protected:
   TbCbVectorMpiPool<T>* mpiVectorPool;
};

#endif //VF_MPI

#endif //TBTRANSMITTERMPIPOOL_H
 
