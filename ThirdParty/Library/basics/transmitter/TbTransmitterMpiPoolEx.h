//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef TBTRANSMITTERMPIPOOLEX_H
#define TBTRANSMITTERMPIPOOLEX_H

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

#include <boost/shared_ptr.hpp>

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

template<typename T> class TbCbVectorSenderMpiPoolEx;
template<typename T> class TbCbVectorReceiverMpiPoolEx;

/*==================================================================*/
template<typename T>
class TbCbVectorMpiPoolEx : public CbVectorPool<T>
{
public:
   typedef boost::shared_ptr< TbCbVectorMpiPoolEx< T > > MpiPoolPtr;

   //////////////////////////////////////////////////////////////////////////
   typedef std::map< int, MpiPoolPtr >      MpiPoolPtrMap;
   typedef typename MpiPoolPtrMap::iterator MpiPoolPtrMapIter;

   //da BasisKlasse templateKlasse ist MUSS man hier die typedefs nochmal wiederholen!
   typedef typename CbVector<T>::value_type value_type;
   typedef typename CbVector<T>::size_type  size_type;
   typedef std::vector< value_type >        Pool;

   typedef unsigned CbVectorKey;
   typedef std::map< CbVectorKey, CbVector< value_type >* /*ptrVector*/  > CbVectorMap;
   typedef typename CbVectorMap::iterator CbVectorMapIter;

   //////////////////////////////////////////////////////////////////////////
   friend class TbCbVectorSenderMpiPoolEx< T >; 
   friend class TbCbVectorReceiverMpiPoolEx< T >; 

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
#ifdef USE_MPI_CXX_SYNTAX 
   static MpiPoolPtr createTbCbVectorMpiPool(const int& poolKey, const int& mpiRemoteRank, const int& mpiTag, MPI::Intracomm comm, const size_type& startPoolSize = 20000 ) //startPoolSize*sizeof(T)/1024/1024 [MB]
#else
   static MpiPoolPtr createTbCbVectorMpiPool(const int& poolKey, const int& mpiRemoteRank, const int& mpiTag, MPI_Comm comm, const size_type& startPoolSize = 20000 ) //startPoolSize*sizeof(T)/1024/1024 [MB]
#endif 
   {
      if( poolMap.find(poolKey)!=poolMap.end() )
      {
         throw UbException(UB_EXARGS,"es ist bereits ein Pool mit dem key vorhanden!!!");
      }

      //pool erstellen
      MpiPoolPtr mpiPool(new TbCbVectorMpiPoolEx<T>(poolKey, mpiRemoteRank, mpiTag, comm, startPoolSize) ); 

      //pool "speichern"
      TbCbVectorMpiPoolEx< value_type >::poolMap[poolKey] = mpiPool;

      return mpiPool; 
   }
   static void deleteTbCbVectorMpiPool(const int& poolKey)
   {
      MpiPoolPtrMapIter it = TbCbVectorMpiPoolEx< value_type >::poolMap.find(poolKey);
      if( it==poolMap.end() )
      {
         throw UbException(UB_EXARGS,"kein Pool mit dem key vorhanden");
      }
      TbCbVectorMpiPoolEx< value_type >::poolMap.erase(it);
   }
   /*==================================================================*/
   static MpiPoolPtr getTbCbVectorMpiPool(const int& poolKey)
   {
      MpiPoolPtrMapIter it;
      if( (it=TbCbVectorMpiPoolEx< T >::poolMap.find(poolKey))!=TbCbVectorMpiPoolEx< T >::poolMap.end() ) 
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
#ifdef USE_MPI_CXX_SYNTAX 
   TbCbVectorMpiPoolEx(const int& poolKey, const int& mpiRemoteRank, const int& mpiTag, MPI::Intracomm comm, const size_type& startPoolSize )
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
      , receiveRequest(MPI::REQUEST_NULL)
      , sendRequest(MPI::REQUEST_NULL)
      , mpiRemoteRank(mpiRemoteRank)               
      , mpiTag(mpiTag)                              
   {
      if     ( (std::string)typeid(value_type).name()==(std::string)typeid(double).name() ) mpiDataType = MPI::DOUBLE;
      else if( (std::string)typeid(value_type).name()==(std::string)typeid(float).name()  ) mpiDataType = MPI::FLOAT;
      else if( (std::string)typeid(value_type).name()==(std::string)typeid(int).name()    ) mpiDataType = MPI::INT;
      else throw UbException(UB_EXARGS,"no MpiDataType for T"+(std::string)typeid(T).name());
   }
#else
   TbCbVectorMpiPoolEx(const int& poolKey, const int& mpiRemoteRank, const int& mpiTag, MPI_Comm comm, const size_type& startPoolSize )
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
      , sendRequest(MPI_REQUEST_NULL)
      , mpiRemoteRank(mpiRemoteRank)               
      , mpiTag(mpiTag)                              
   {
      if     ( (std::string)typeid(value_type).name()==(std::string)typeid(double).name() ) mpiDataType = MPI_DOUBLE;
      else if( (std::string)typeid(value_type).name()==(std::string)typeid(float).name()  ) mpiDataType = MPI_FLOAT;
      else if( (std::string)typeid(value_type).name()==(std::string)typeid(int).name()    ) mpiDataType = MPI_INT;
      else throw UbException(UB_EXARGS,"no MpiDataType for T"+(std::string)typeid(T).name());
   }
#endif


public:
   /*==================================================================*/
   //returns key of Pool in MpiPoolMap
   int  getPoolKey()    const { return  this->poolKey;       }
   /*==================================================================*/
   //returns rank of process pool data will be send to/received from
   int  getRemoteRank() const { return  this->mpiRemoteRank; }
   /*==================================================================*/
   //returns tag of process pool data will be send to/received from
   int  getRemoteTag()  const { return  this->mpiTag;        }

protected:
   /*==================================================================*/
   /*==================================================================*/
   /*==================================================================*/
   /*==================================================================*/
   /*==================================================================*/
   //void prepareTbReceiveDataOrder() 
   //{
   //counterPrepareReceiveDataOrder++;
   //if(counterPrepareReceiveDataOrder==nofStoredVectors)
   //{
   //   unsigned nofElements = relationMap.size()*3+1; //map MUSS auf beiden seiten gleich gross sein, sonst hat man ein grundsaetzliches problem ;-)
   //   tmpOrderVec.resize(nofElements);
   //   std::cout<<RcfSystem::getRank()<<" prepForRecv from rank="<<mpiRemoteRank<<" with tag="<<mpiTag<<"e="<<nofElements<<std::endl;
   //   receiveRequest = comm.Irecv(&tmpOrderVec[0], nofElements, MPI::UNSIGNED, mpiRemoteRank, mpiTag);
   //   counterPrepareReceiveDataOrder = 0;
   //}
   //}
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
            CbVectorKey vectorKey=0;
            size_type   dataSize=0, startIndexInPool=0;
            this->getCbVectorData(*it->second/*vec*/, vectorKey, startIndexInPool, dataSize );
            if(it->first != vectorKey) throw UbException(UB_EXARGS,"key mismatch!");

            tmpSendOrderVec[index++] = (unsigned)vectorKey;         //vectorKey == allocator.getAllocatorKey()
            tmpSendOrderVec[index++] = (unsigned)startIndexInPool;  //startIndex in poolVector
            tmpSendOrderVec[index++] = (unsigned)dataSize;          //dataSize
         }
         //std::cout<<RcfSystem::getRank()<<" send to rank="<<mpiRemoteRank<<" with tag="<<mpiTag<<" e="<<nofElements<<std::endl;
         //comm.Send(&tmpOrderVec[0],nofElements, MPI::UNSIGNED, mpiRemoteRank, mpiTag);

         ////////////////////////////
         //int rank;
         //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
         //std::cout<<"rank = " << rank <<" sendDataOrder() nofElements = "<<nofElements<<std::endl;
         ////////////////////////////

#ifdef USE_MPI_CXX_SYNTAX 
         sendRequest = comm.Isend(&tmpSendOrderVec[0],(int)tmpSendOrderVec.size(), MPI::UNSIGNED, mpiRemoteRank, mpiTag);
#else
         MPI_Isend(&tmpSendOrderVec[0],(int)tmpSendOrderVec.size(), MPI_UNSIGNED, mpiRemoteRank, mpiTag, comm, &sendRequest);
#endif
         counterSendDataOrder=0;

         nofStoredVectors = this->cbVectorMap.size();

         UBLOG(logDEBUG5, "TbCbVectorMpiPool::sendDataOrder()"<<" mpiRemoteRank="<<mpiRemoteRank<<" mpiTag="<<mpiTag);

#ifdef DEBUG
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

         //////////////TODO////////////DEBUG
         //int rank;
         //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
         //std::cout<<"rank = " << rank <<" receiveDataOrder() nofElements = "<<nofElements << " from " << mpiRemoteRank <<std::endl;
         UBLOG(logDEBUG5, "TbCbVectorMpiPool::receiveDataOrder()"<<" mpiRemoteRank="<<mpiRemoteRank<<" mpiTag="<<mpiTag);
         //////////////////////////

         std::vector< unsigned > tmpRecvOrderVec;
         tmpRecvOrderVec.resize(nofElements);

#ifdef USE_MPI_CXX_SYNTAX 
         comm.Recv(&tmpRecvOrderVec[0], nofElements, MPI::UNSIGNED, mpiRemoteRank, mpiTag);
#else
         //MPI_Status status;
         MPI_Recv(&tmpRecvOrderVec[0], nofElements, MPI_UNSIGNED, mpiRemoteRank, mpiTag, comm, MPI_STATUS_IGNORE);
#endif

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

#ifdef DEBUG
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
#ifdef USE_MPI_CXX_SYNTAX 
         if(sendRequest != MPI::REQUEST_NULL) sendRequest.Wait();
#else
         //if(sendRequest != MPI_REQUEST_NULL) MPI_Wait(&sendRequest, MPI_STATUS_IGNORE);
#endif
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
#ifdef DEBUG
         if(this->orgPoolVectorStartPointer != &this->pool[0] ) throw UbException(UB_EXARGS, "ups, pool array adress changed - unknown behavoir");
#endif

         //standard send
#ifdef USE_MPI_CXX_SYNTAX 
         comm.Send(&this->pool[0],(int)this->nextCbVectorStartIndexInPool, mpiDataType, mpiRemoteRank, mpiTag);
#else
         //MPI_Send(&this->pool[0],(int)this->nextCbVectorStartIndexInPool, mpiDataType, mpiRemoteRank, mpiTag, comm);
         MPI_Isend(&this->pool[0],(int)this->nextCbVectorStartIndexInPool, mpiDataType, mpiRemoteRank, mpiTag, comm, &sendRequest);
#endif
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
#ifdef DEBUG
         if(this->orgPoolVectorStartPointer != &this->pool[0] ) throw UbException(UB_EXARGS, "ups, pool array adress changed - unknown behavoir");
#endif
         //std::cout<<"Irecv von "<<(int)nextStartIndex<<"elementen von "<<mpiRemoteRank<<" mit tag="<<mpiTag<<std::endl;
#ifdef USE_MPI_CXX_SYNTAX 
         receiveRequest = comm.Irecv(&this->pool[0],(int)this->nextCbVectorStartIndexInPool, mpiDataType, mpiRemoteRank, mpiTag);
#else
         //MPI_Irecv(&this->pool[0],(int)this->nextCbVectorStartIndexInPool, mpiDataType, mpiRemoteRank, mpiTag, comm, &receiveRequest);
#endif
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
#ifdef USE_MPI_CXX_SYNTAX 
         receiveRequest.Wait();
#else
         //MPI_Wait(&receiveRequest, MPI_STATUS_IGNORE);
         MPI_Recv(&this->pool[0],(int)this->nextCbVectorStartIndexInPool, mpiDataType, mpiRemoteRank, mpiTag, comm, MPI_STATUS_IGNORE);
#endif
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

#ifdef USE_MPI_CXX_SYNTAX 
   MPI::Intracomm comm;
   MPI::Request   receiveRequest;
   MPI::Request   sendRequest;
   MPI::Datatype  mpiDataType;
#else
   MPI_Comm     comm;
   MPI_Request  receiveRequest;
   MPI_Request  sendRequest;
   //MPI_Status   sendStatus;
   //MPI_Status   receiveStatus;
   MPI_Datatype mpiDataType;
#endif

   int mpiRemoteRank, mpiTag;

#ifdef DEBUG
   T* orgPoolVectorStartPointer;
#endif
};

template<typename T>
typename TbCbVectorMpiPoolEx<T>::MpiPoolPtrMap TbCbVectorMpiPoolEx<T>::poolMap;

//   static MpiPoolPtrMap poolMap;


//////////////////////////////////////////////////////////////////////////
//  TbSenderMpiPool
//////////////////////////////////////////////////////////////////////////
template<typename T>
class TbCbVectorSenderMpiPoolEx : public TbTransmitter< CbVector< T >  >
{
public:
   typedef CbVector< T > value_type;   

public:
   TbCbVectorSenderMpiPoolEx(const unsigned int& cbVectorKey, TbCbVectorMpiPoolEx< T >* mpiVectorPool)
      : mpiVectorPool(mpiVectorPool)
   { 
      this->getData().setAllocator( new CbVectorAllocatorPool<T>(cbVectorKey,this->mpiVectorPool) );
   }
   ~TbCbVectorSenderMpiPoolEx()
   {
      if( this->mpiVectorPool->getNofStoredVectors()==1 ) //last entry!
      {
         TbCbVectorMpiPoolEx< T >::deleteTbCbVectorMpiPool(this->mpiVectorPool->getPoolKey());  
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
   TbCbVectorMpiPoolEx<T>* mpiVectorPool;
};


/*==================================================================*/
template<typename T>
class TbCbVectorReceiverMpiPoolEx : public TbTransmitter< CbVector< T >  >
{
public:
   typedef CbVector< T > value_type;   

public:
   TbCbVectorReceiverMpiPoolEx(const unsigned int& cbVectorKey, TbCbVectorMpiPoolEx< T >* mpiVectorPool)
      : mpiVectorPool(mpiVectorPool)
   { 
      this->getData().setAllocator( new CbVectorAllocatorPool<T>(cbVectorKey, this->mpiVectorPool) );
   }
   ~TbCbVectorReceiverMpiPoolEx()
   {
      if( this->mpiVectorPool->getNofStoredVectors()==1 ) //last entry!
      {
         TbCbVectorMpiPoolEx< T >::deleteTbCbVectorMpiPool(this->mpiVectorPool->getPoolKey());  
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
   TbCbVectorMpiPoolEx<T>* mpiVectorPool;
};

#endif //VF_MPI

#endif //TBTRANSMITTERMPIPOOL_H
 
