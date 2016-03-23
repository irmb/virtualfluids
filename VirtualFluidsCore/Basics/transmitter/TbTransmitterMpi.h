//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef TBTRANSMITTERMPI_H
#define TBTRANSMITTERMPI_H

#ifdef VF_MPI

/*=========================================================================*/
/*  MPI Transmitter                                                        */
/*                                                                         */
/**
This Class provides the base for exception handling.
Old TbTransmitter was renamed in TbTransmitterCPPB (C++ binding in MPI is deprecated)
Rewrite from K. Kucher with C binding
<BR><BR>
@author <A HREF="mailto:kucher@irmb.tu-bs.de">K. Kucher</A>
@version 1.0 - 21.11.11
*/ 

/*
usage: ...
*/

#include <iostream>
#include <mpi.h>
#include <basics/transmitter/TbTransmitter.h>


//////////////////////////////////////////////////////////////////////////
// TbVectorSenderMpiUnblocked
template< typename Vector  >
class TbVectorSenderMpiUnblocked : public TbTransmitter< Vector >
{
public:
   typedef Vector value_type;

public:
   TbVectorSenderMpiUnblocked(const int& sendTbRank, const int& sendTag, MPI_Comm comm) 
      : comm(comm), request(MPI_REQUEST_NULL), sendTbRank(sendTbRank), sendTag(sendTag), dataSize(0)    
   { 
      //temporaeren vector erstellen um den datentyp zu ermittlen
      if     ( (std::string)typeid(typename Vector::value_type).name()==(std::string)typeid(double).name() ) mpiDataType = MPI_DOUBLE;
      else if( (std::string)typeid(typename Vector::value_type).name()==(std::string)typeid(float).name()  ) mpiDataType = MPI_FLOAT;
      else if( (std::string)typeid(typename Vector::value_type).name()==(std::string)typeid(int).name()    ) mpiDataType = MPI_INT;
      else UB_THROW( UbException(UB_EXARGS,"no MpiDataType for type="+(std::string)typeid(typename Vector::value_type).name()) );
   }

   bool isLocalTransmitter()  const { return false;                        }
   bool isRemoteTransmitter() const { return !this->isLocalTransmitter();  }

   void sendDataSize()   
   {
      dataSize = (unsigned int)this->getData().size();
      //MPI_Isend(&dataSize, 1, MPI_UNSIGNED, sendTbRank, sendTag, comm, &request);
      MPI_Send(&dataSize, 1, MPI_UNSIGNED, sendTbRank, sendTag, comm);
   }  
   void receiveDataSize()   { UB_THROW( UbException(UB_EXARGS,"MPIVectorSender sends only") ); }
   void prepareForSend()    { if(request!=MPI_REQUEST_NULL) MPI_Wait(&request, &status); }
   void sendData()          { MPI_Isend(&this->getData()[0],(int)this->getData().size(), mpiDataType, sendTbRank, sendTag, comm, &request);  }

   void prepareForReceive() { UB_THROW( UbException(UB_EXARGS,"MPIVectorSender sends only") ); }
   Vector& receiveData()    { UB_THROW( UbException(UB_EXARGS,"MPIVectorSender sends only") ); }

   //info-section (usable for remote transmitter)
   int  getSendTbRank()   const { return  sendTbRank; }
   int  getSendTbTag()    const { return  sendTag;    }
   int  getRecvFromRank() const { UB_THROW( UbException(UB_EXARGS,"MPIVectorSender sends only") ); }
   int  getRecvFromTag()  const { UB_THROW( UbException(UB_EXARGS,"MPIVectorSender sends only") ); }

   std::string toString()  const { return "TbVectorSenderMpiUnblocked<"+(std::string)typeid(Vector).name()+"<"+(std::string)typeid(typename Vector::value_type).name()+"> > to rank (tag)"+UbSystem::toString(sendTbRank)+"("+UbSystem::toString(sendTag)+")"; }

protected:
   MPI_Comm comm;
   MPI_Request   request;
   MPI_Datatype  mpiDataType;
   MPI_Status    status;
   int sendTbRank, sendTag;
   unsigned dataSize;
};


//////////////////////////////////////////////////////////////////////////
// TbVectorSenderMpiBlocked
template< typename Vector  >
class TbVectorSenderMpiBlocked : public TbTransmitter< Vector >
{
public:
   typedef Vector value_type;

public:
   TbVectorSenderMpiBlocked(const int& sendTbRank, const int& sendTag, MPI_Comm comm) 
      : comm(comm), request(MPI_REQUEST_NULL), sendTbRank(sendTbRank), sendTag(sendTag), dataSize(0)    
   { 
      if     ( (std::string)typeid(typename Vector::value_type).name()==(std::string)typeid(double).name() ) mpiDataType = MPI_DOUBLE;
      else if( (std::string)typeid(typename Vector::value_type).name()==(std::string)typeid(float).name()  ) mpiDataType = MPI_FLOAT;
      else if( (std::string)typeid(typename Vector::value_type).name()==(std::string)typeid(int).name()    ) mpiDataType = MPI_INT;
      else UB_THROW( UbException(UB_EXARGS,"no MpiDataType for Vector"+(std::string)typeid(Vector).name()) );
   }

   bool isLocalTransmitter()  const { return false;                        }
   bool isRemoteTransmitter() const { return !this->isLocalTransmitter();  }

   void sendDataSize()   
   {
      dataSize = (unsigned int)this->getData().size();
      MPI_Isend(&dataSize, 1, MPI_UNSIGNED, sendTbRank, sendTag, comm, &request);
   }  
   void receiveDataSize()   { UB_THROW( UbException(UB_EXARGS,"MPIVectorSender sends only") ); }
   void sendData()          { MPI_Wait(&request, &status); MPI_Send(&this->getData()[0],(int)this->getData().size(), mpiDataType, sendTbRank, sendTag, comm); }

   void prepareForReceive()  { UB_THROW( UbException(UB_EXARGS,"MPIVectorSender sends only") ); }
   value_type& receiveData() { UB_THROW( UbException(UB_EXARGS,"MPIVectorSender sends only") ); }

   //info-section (usable for remote transmitter)
   int  getSendTbRank()   const { return  sendTbRank; }
   int  getSendTbTag()    const { return  sendTag;    }
   int  getRecvFromRank() const { UB_THROW( UbException(UB_EXARGS,"MPIVectorSender sends only") ); }
   int  getRecvFromTag()  const { UB_THROW( UbException(UB_EXARGS,"MPIVectorSender sends only") ); }

   std::string toString() const { return "TbVectorSenderMpiBlocked<"+(std::string)typeid(Vector).name()+"<"+(std::string)typeid(typename Vector::value_type).name()+"> > to rank (tag)"+UbSystem::toString(sendTbRank)+"("+UbSystem::toString(sendTag)+")"; }

protected:
   MPI_Comm comm;
   MPI_Request   request;
   MPI_Datatype  mpiDataType;
   MPI_Status    status;
   int sendTbRank, sendTag;
   unsigned dataSize;
};

//////////////////////////////////////////////////////////////////////////
// TbVectorReceiverMpiUnblocked
template<typename Vector  >
class TbVectorReceiverMpiUnblocked : public TbTransmitter< Vector >
{
public:
   typedef Vector value_type;

public:
   TbVectorReceiverMpiUnblocked(const int& receiveFromRank, const int& receiveTag, MPI_Comm comm) 
      : comm(comm), request(MPI_REQUEST_NULL), receiveFromRank(receiveFromRank), receiveTag(receiveTag)       
   { 
      if     ( (std::string)typeid(typename Vector::value_type).name()==(std::string)typeid(double).name() ) mpiDataType = MPI_DOUBLE;
      else if( (std::string)typeid(typename Vector::value_type).name()==(std::string)typeid(float).name()  ) mpiDataType = MPI_FLOAT;
      else if( (std::string)typeid(typename Vector::value_type).name()==(std::string)typeid(int).name()    ) mpiDataType = MPI_INT;
      else UB_THROW( UbException(UB_EXARGS,"no MpiDataType for Vector"+(std::string)typeid(Vector).name()) );
   }

   bool isLocalTransmitter()  const { return false;                        }
   bool isRemoteTransmitter() const { return !this->isLocalTransmitter();  }

   void sendDataSize()     { UB_THROW( UbException(UB_EXARGS,"MPIVectorReceiver receives only") ); }
   void receiveDataSize()   
   {
      unsigned dataSize;
      MPI_Recv(&dataSize, 1, MPI_UNSIGNED, receiveFromRank, receiveTag, comm, &status);
      this->getData().resize(dataSize,0.0);
   }  
   void sendData()          { UB_THROW( UbException(UB_EXARGS,"MPIVectorReceiver receives only") ); }

   void prepareForReceive() { MPI_Irecv(&this->getData()[0],(int)this->getData().size(), mpiDataType, receiveFromRank, receiveTag, comm, &request);  }
   Vector& receiveData()    { MPI_Wait(&request, &status); return this->getData(); }

   //info-section (usable for remote transmitter)
   int  getSendTbRank()   const { UB_THROW( UbException(UB_EXARGS,"MPIVectorReceiver receives only") ); }
   int  getSendTbTag()    const { UB_THROW( UbException(UB_EXARGS,"MPIVectorReceiver receives only") ); }
   int  getRecvFromRank() const { return  receiveFromRank; }
   int  getRecvFromTag()  const { return  receiveTag;      }

   std::string toString() const { return "TbVectorReceiverMpiUnblocked<"+(std::string)typeid(Vector).name()+"<"+(std::string)typeid(typename Vector::value_type).name()+"> > to rank (tag)"+UbSystem::toString(receiveFromRank)+"("+UbSystem::toString(receiveTag)+")"; }

protected:
   MPI_Comm comm;
   MPI_Request   request;
   MPI_Datatype  mpiDataType;
   MPI_Status    status;
   int receiveFromRank, receiveTag;
};


//////////////////////////////////////////////////////////////////////////
template<typename Vector>
class TbVectorReceiverMpiBlocked : public TbTransmitter< Vector >
{
public:
   typedef Vector value_type;

public:
   TbVectorReceiverMpiBlocked(const int& receiveFromRank, const int& receiveTag, MPI_Comm comm) 
      : comm(comm), request(MPI_REQUEST_NULL), receiveFromRank(receiveFromRank), receiveTag(receiveTag)
   { 
      if     ( (std::string)typeid(typename Vector::value_type).name()==(std::string)typeid(double).name() ) mpiDataType = MPI_DOUBLE;
      else if( (std::string)typeid(typename Vector::value_type).name()==(std::string)typeid(float).name()  ) mpiDataType = MPI_FLOAT;
      else if( (std::string)typeid(typename Vector::value_type).name()==(std::string)typeid(int).name()    ) mpiDataType = MPI_INT;
      else UB_THROW( UbException(UB_EXARGS,"no MpiDataType for Vector+(std::string)typeid(Vector).name()") );
   }

   bool isLocalTransmitter()  const { return false;                        }
   bool isRemoteTransmitter() const { return !this->isLocalTransmitter();  }

   void sendDataSize()     { UB_THROW( UbException(UB_EXARGS,"MPIVectorReceiver receives only") ); }
   void receiveDataSize()   
   {
      unsigned dataSize;
      MPI_Recv(&dataSize, 1, MPI_UNSIGNED, receiveFromRank, receiveTag, comm, &status);
      this->getData().resize(dataSize,0.0);
   }  
   void sendData()         { UB_THROW( UbException(UB_EXARGS,"MPIVectorReceiver receives only") ); }
   Vector& receiveData()
   {
      MPI_Recv(&this->getData()[0],(int)this->getData().size(), mpiDataType, receiveFromRank, receiveTag, comm, &status);
      return this->getData();
   }

   //info-section (usable for remote transmitter)
   int  getSendTbRank()   const { UB_THROW( UbException(UB_EXARGS,"MPIVectorReceiver receives only") ); }
   int  getSendTbTag()    const { UB_THROW( UbException(UB_EXARGS,"MPIVectorReceiver receives only") ); }
   int  getRecvFromRank() const { return  receiveFromRank; }
   int  getRecvFromTag()  const { return  receiveTag;      }

   std::string toString() const { return "TbVectorReceiverMpiBlocked<"+(std::string)typeid(Vector).name()+"<"+(std::string)typeid(typename Vector::value_type).name()+"> > to rank (tag)"+UbSystem::toString(receiveFromRank)+"("+UbSystem::toString(receiveTag)+")"; }

protected:
   MPI_Comm comm;
   MPI_Request   request;
   MPI_Datatype  mpiDataType;
   MPI_Status    status;
   int receiveFromRank, receiveTag;
};

#endif //VF_MPI

#endif //TBTRANSMITTERMPI_H
