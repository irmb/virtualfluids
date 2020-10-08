//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef PbMpi_H
#define PbMpi_H

#include <vector>
#include <sstream>

#ifndef VF_MPI
#  error VF_MPI has to be defined
#endif

// As we doing a lot of const-cast here we define PbMpi.h to system_header to mute clang-tidy
#pragma clang system_header

//#undef SEEK_SET
//#undef SEEK_CUR
//#undef SEEK_END
#include <mpi.h>

#include <basics/utilities/UbException.h>

#ifdef USE_MPI_CXX_SYNTAX
   #define PbMpi_COMM_WORLD MPI::COMM_WORLD
   #define PbMpi_INT        MPI::INT  
   #define PbMpi_CHAR       MPI::CHAR 
   #define PbMpi_SHORT      MPI::SHORT
   #define PbMpi_FLOAT      MPI::FLOAT
   #define PbMpi_DOUBLE     MPI::DOUBLE
   #define PbMpi_COMM_NULL  MPI::COMM_NULL   


namespace PbMpi
{
   typedef MPI::Intracomm Comm;
   typedef MPI::Group     Group;
   typedef MPI::Request   Request;
   typedef MPI::Status    Status;

   inline void Init( )  
   {
      MPI::Init(); 
      MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS); 
   }
   inline void Init(int& argc, char** argv )  
   {
      MPI::Init(argc, argv); 
      MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS); 
   }
   inline void Finalize()                      { MPI::Finalize();        }

   inline int  GetCommSize( const Comm& comm ) { return comm.Get_size(); }
   inline int  GetCommRank( const Comm& comm ) { return comm.Get_rank(); }
   inline void Barrier( const Comm& comm  )    { comm.Barrier();         }

   inline double Wtime()                       { return MPI::Wtime();    }
   inline double Wtick()                       { return MPI::Wtick();    }

   inline void Wait( Request& request, Status* outStatus=NULL) 
   { 
      if(outStatus) request.Wait(*outStatus); 
      else request.Wait();   
   }

   inline Group GetCommGroup(Comm& comm)                               { return comm.Get_group();     }
   inline Group GetGroupIncl( Group& group, const int& n, int* ranks ) { return group.Incl(n, ranks); }
   inline Comm  CommCreateComm( Comm& comm, Group& group )             { return comm.Create(group);   }

   inline void Alltoall( Comm& comm, void* sendBuffer, const int& sn, const MPI_Datatype& sdatatype, void* recvBuffer, const int& rn, const MPI_Datatype& rdatatype)
   {
      comm.Alltoall(sendBuffer, sn, sdatatype, recvBuffer, rn, rdatatype);
   }
   inline void Bcast(Comm& comm, void* data, const int& n, const MPI_Datatype& datatype , const int& srcRank )
   {
      comm.Bcast(data, n, datatype, srcRank);
   }
   inline void Send(Comm& comm,  const void* data, const int& length, const MPI_Datatype& dataType, const int& destRank, const int& tag)
   {
      try
      { 
         comm.Send(data, length, dataType, destRank, tag);
      }
      catch(MPI::Exception& e)
      {
         std::stringstream ss; 
         ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
         throw UbException(UB_EXARGS,"MPI:Exception catched\n"+ss.str());
      }
      catch(...)  
      {
         throw UbException(UB_EXARGS,"unknown exception"); 
      }
   }
   inline void Recv(Comm& comm,  const void* data, const int& length, const MPI_Datatype& dataType, const int& srcRank, const int& tag)
   {
      try
      { 
         comm.Recv(const_cast<void*>(data), length, dataType, srcRank, tag);
      }
      catch(MPI::Exception& e)
      {
         std::stringstream ss; 
         ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
         throw UbException(UB_EXARGS,"MPI:Exception catched \n"+ss.str());
      }
      catch(...)  
      {
         throw UbException(UB_EXARGS,"unknown exception"); 
      }
   }

   inline void Irecv(Comm comm,  const void* data, const int& length, const MPI_Datatype& dataType, const int& srcRank, const int& tag, Request& outRequest)
   {
      outRequest = comm.Irecv(const_cast<void*>(data), length, dataType, srcRank, tag);
   }
   inline void Ssend(Comm& comm,  const void* data, const int& length, const MPI_Datatype& dataType, const int& destRank, const int& tag)
   {
      try
      { 
         comm.Ssend(data, length, dataType, destRank, tag);
      }
      catch(MPI::Exception& e)
      {
         std::stringstream ss; 
         ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
         throw UbException(UB_EXARGS,"MPI:Exception catched\n"+ss.str());
      }
      catch(...)  
      {
         throw UbException(UB_EXARGS,"unknown exception"); 
      }
   }

}
#else //////////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////
   // C-Syntax 
   //////////////////////////////////////////////////////////////////////////
   namespace PbMpi
   {
      using Comm = MPI_Comm;
      using Group = MPI_Group;
      using Request = MPI_Request;
      using Status = MPI_Status;
   }

   #define PbMpi_COMM_WORLD ((PbMpi::Comm)MPI_COMM_WORLD)
   #define PbMpi_INT        MPI_INT
   #define PbMpi_CHAR       MPI_CHAR
   #define PbMpi_SHORT      MPI_SHORT
   #define PbMpi_FLOAT      MPI_FLOAT
   #define PbMpi_DOUBLE     MPI_DOUBLE
   #define PbMpi_COMM_NULL  MPI_COMM_NULL   

namespace PbMpi
{
   inline void Init( )  
   {
      int    argc = 1;
      char** argv = new char*[1];
      argv[0]     = new char[1];
      argv[0][0]  = 'n';
      MPI_Init(&argc, &argv); 
   }
   inline void Init(int& argc, char** argv ) { MPI_Init(&argc, &argv);                         }
   inline void Finalize( )                   { MPI_Finalize();                                 }
   inline int  GetCommSize( Comm comm )      { int tmp; MPI_Comm_size(comm, &tmp); return tmp; }
   inline int  GetCommRank( Comm comm )      { int tmp; MPI_Comm_rank(comm, &tmp); return tmp; }
   inline void Barrier(  Comm comm )         { MPI_Barrier( comm );                            }
   inline double Wtime()                     { return MPI_Wtime();                             }
   inline double Wtick()                     { return MPI_Wtick();                             }
   inline void   Wait( Request& request, Status* outStatus=NULL) { MPI_Wait( &request, outStatus);   }

   inline Group GetCommGroup(Comm comm)                               { Group out; MPI_Comm_group(comm, &out);            return out; }
   inline Group GetGroupIncl( Group group, const int& n, int* ranks ) { Group out; MPI_Group_incl(group, n, ranks, &out); return out; }
   inline Comm  CommCreateComm( Comm comm, Group& group )             { Comm out;  MPI_Comm_create(comm, group, &out);    return out; }

   inline void Alltoall( Comm comm, void* sendBuffer, const int& sn, const MPI_Datatype& sdatatype, void* recvBuffer, const int& rn, const MPI_Datatype& rdatatype)
   {
      MPI_Alltoall(sendBuffer, sn, sdatatype, recvBuffer, rn, rdatatype, comm);
   }
   inline void Bcast(Comm comm, void* data, const int& n, const MPI_Datatype& datatype , const int& srcRank )
   {
      MPI_Bcast(data, n, datatype, srcRank, comm);
   }
   inline void Send(Comm comm,  const void* data, const int& length, const MPI_Datatype& dataType, const int& destRank, const int& tag)
   {
      MPI_Send(const_cast<void*>(data), length, dataType, destRank, tag, comm);
   }
   inline void Recv(Comm comm,  const void* data, const int& length, const MPI_Datatype& dataType, const int& srcRank, const int& tag)
   {
      MPI_Recv(const_cast<void*>(data), length, dataType, srcRank, tag, comm, MPI_STATUS_IGNORE);
   }
   inline void Ssend(Comm comm,  const void* data, const int& length, const MPI_Datatype& dataType, const int& destRank, const int& tag)
   {
      MPI_Ssend(const_cast<void*>(data), length, dataType, destRank, tag, comm);
   }
   inline void Irecv(Comm comm,  const void* data, const int& length, const MPI_Datatype& dataType, const int& srcRank, const int& tag, Request& outRequest)
   {
      MPI_Irecv(const_cast<void*>(data), length, dataType, srcRank, tag, comm, &outRequest);
   }

}
#endif

namespace PbMpi
{
   /*======================================================================*/  
   // send a single value "value" of MPI_Datatype
   template <class T>
   inline void sendSingleValue(const T& value, MPI_Datatype datatype, int dest, int tag, PbMpi::Comm comm);
   
   /*======================================================================*/  
   // receives a single value "value" of MPI_Datatype
   template <class T>
   inline void receiveSingleValue(T& value, MPI_Datatype datatype, int source, int tag, PbMpi::Comm comm);
   
   /*======================================================================*/  
   // receives and returns a single value of MPI_Datatype
   // expample: int value = PbMpi::receiveSingleValue<int>(MPI::INT,0,10,comm);
   template <class T>
   inline T receiveSingleValue(MPI_Datatype datatype, int source, int tag, PbMpi::Comm comm);
   
   /*======================================================================*/  
   // sends bool value (doesn't work with template, why ever... stupid MPI)
   inline void sendBoolValue(const bool& value,int dest, int tag, PbMpi::Comm comm);

   /*======================================================================*/  
   // receives bool value (doesn't work with template, why ever... stupid MPI)
   inline bool receiveBoolValue(int source, int tag, PbMpi::Comm comm);

   /*======================================================================*/  
   // sends bool value (doesn't work with template, why ever... stupid MPI)
   inline void sendStringValue(const std::string& value,int dest, int tag, PbMpi::Comm comm);

   /*======================================================================*/  
   // receives bool value (doesn't work with template, why ever... stupid MPI)
   inline std::string receiveStringValue(int source, int tag, PbMpi::Comm comm);

   /*======================================================================*/  
   // send a vector of MPI_Datatype
   template <class T>
	inline void sendVector(const std::vector<T>& v, MPI_Datatype datatype, int dest, int tag, PbMpi::Comm comm);
	
   /*======================================================================*/  
   // receive a std::vector of MPI_Datatype
   template <class T>
   inline void receiveVector(std::vector<T>& v, MPI_Datatype datatype, int source, int tag, PbMpi::Comm comm);

   /*======================================================================*/  
   // receive a vector of MPI_Datatype and adds this vector to existing vector
   // ans returns number of received elements
   template <class T>
   inline int receiveVectorAndAddToVector(std::vector<T>& v, MPI_Datatype datatype, int source, int tag, PbMpi::Comm comm);

   /*======================================================================*/  
   // send a std::vector of strings
   inline void sendStringVector(const std::vector<std::string>& v, int dest, int tag, PbMpi::Comm comm);

   /*======================================================================*/  
   // send a vector of strings
   inline void receiveStringVector(std::vector<std::string>& v, int dest, int tag, PbMpi::Comm comm);
};

/*======================================================================*/  
// send a single value of MPI_Datatype
template <class T>
void PbMpi::sendSingleValue(const T& value, MPI_Datatype datatype, int dest, int tag, PbMpi::Comm comm)
{
   PbMpi::Send(comm, &value, 1, datatype, dest, tag); 
   //comm.Send(&value, 1, datatype, dest, tag); 
}
/*======================================================================*/  
template <class T>
void PbMpi::receiveSingleValue(T& value, MPI_Datatype datatype, int source, int tag, PbMpi::Comm comm) 
{
   PbMpi::Recv(comm, &value, 1, datatype, source, tag); 
   //comm.Recv(&value, 1, datatype, source, tag); 
}
/*======================================================================*/  
template <class T>
T PbMpi::receiveSingleValue(MPI_Datatype datatype, int source, int tag, PbMpi::Comm comm) 
{
   T value;
   PbMpi::Recv(comm, &value, 1, datatype, source, tag); 
   //comm.Recv(&value, 1, datatype, source, tag); 

   return value;
}
/*======================================================================*/  
// send a bool value (bool doesn't work with template, why ever)
void PbMpi::sendBoolValue(const bool& value,int dest, int tag, PbMpi::Comm comm)
{
   short dummy;
   if(value) dummy=1;                  
   else      dummy=0;

   PbMpi::Send(comm, &dummy, 1, PbMpi_SHORT, dest, tag); 
   //comm.Send(&dummy, 1, MPI::SHORT, dest, tag); 
}
/*======================================================================*/  
bool PbMpi::receiveBoolValue(int source, int tag, PbMpi::Comm comm) 
{
   short dummy {0};
   PbMpi::Recv(comm, &dummy, 1, PbMpi_SHORT, source, tag); 
   //comm.Recv(&dummy, 1, MPI::SHORT, source, tag);
 
   return (dummy==1);
}
/*======================================================================*/  
// sends bool value (doesn't work with template, why ever... stupid MPI)
void PbMpi::sendStringValue(const std::string& value,int dest, int tag, PbMpi::Comm comm)
{
   std::vector<char> vec;
   for(char i : value)
      vec.push_back(i);
 
   PbMpi::sendVector(vec,PbMpi_CHAR,dest,tag,comm);
}

/*======================================================================*/  
// receives bool value (doesn't work with template, why ever... stupid MPI)
std::string PbMpi::receiveStringValue(int source, int tag, PbMpi::Comm comm)
{
   std::vector<char> vec;
   PbMpi::receiveVector(vec,PbMpi_CHAR,source,tag,comm);
   
   std::string str;
   for(char i : vec)
      str+=i;

   return str;
}
/*======================================================================*/  
// send a vector of MPI_Datatype
template <class T>
void PbMpi::sendVector(const std::vector<T>& v, MPI_Datatype datatype, int dest, int tag, PbMpi::Comm comm)
{
   // send size
   int size = (int)v.size();
   
   PbMpi::Send(comm, &size, 1, PbMpi_INT, dest, tag);
   //comm.Send(&size, 1, MPI::INT, dest, tag); 
   
   if(size>0)
	{
      PbMpi::Send(comm, &v[0], size, datatype, dest, tag);
      //comm.Send(&v[0], size, datatype, dest, tag);
   }
}
/*======================================================================*/  
// receive a vector of MPI_Datatype
template <class T>
void PbMpi::receiveVector(std::vector<T>& v, MPI_Datatype datatype, int source, int tag, PbMpi::Comm comm) 
{
   int size {0};

   PbMpi::Recv(comm, &size, 1, PbMpi_INT, source, tag);
   //comm.Recv(&size, 1, MPI::INT, source, tag); 

   v.resize(size);

   if( size>0 )
   {
      PbMpi::Recv(comm, &v[0], size, datatype, source, tag);
      //comm.Recv(&v[0], size, datatype, source, tag); 
   }
}
/*======================================================================*/  
// receive a vector of MPI_Datatype and adds this vector to existing vector
// return value is size of received elements
template <class T>
int PbMpi::receiveVectorAndAddToVector(std::vector<T>& v, MPI_Datatype datatype, int source, int tag, PbMpi::Comm comm) 
{
   int incommingSize;

   PbMpi::Recv(comm, &incommingSize, 1, PbMpi_INT, source, tag);
   //comm.Recv(&incommingSize, 1, MPI::INT, source, tag);

   int oldSize = (int)v.size();
   v.resize(oldSize+incommingSize);

   if( incommingSize>0 )
   {
      PbMpi::Recv(comm, &v[oldSize], incommingSize, datatype, source, tag);
      //comm.Recv(&v[oldSize], incommingSize, datatype, source, tag);
   }

   return incommingSize;
}
/*======================================================================*/  
// send a vector of strings
void PbMpi::sendStringVector(const std::vector<std::string>& v, int dest, int tag, PbMpi::Comm comm)
{
   // send size
   int stringVectorSize = (int)v.size();

   PbMpi::Send(comm, &stringVectorSize, 1, PbMpi_INT, dest, tag); 
   //comm.Send(&stringVectorSize, 1, MPI::INT, dest, tag); 

   if(stringVectorSize>0)
   {
      std::vector<int> singleStringSizes(stringVectorSize+1);
      int nofChars = 0;
      for(int i=0; i<stringVectorSize; i++)
         nofChars += singleStringSizes[i] = (int)v[i].length();
      singleStringSizes[stringVectorSize] = nofChars;

      PbMpi::Send(comm, &singleStringSizes[0], stringVectorSize+1, PbMpi_INT, dest, tag); 

      std::vector<char> charVector(nofChars);
      int pos = 0;
      for(int i=0; i<stringVectorSize; i++)
         for(int j=0; j<singleStringSizes[i]; j++)
            charVector[pos++] = v[i][j];      

      PbMpi::Send(comm, &charVector[0], nofChars, PbMpi_CHAR, dest, tag); 
      //comm.Send(&charVector[0], nofChars, MPI::CHAR, dest, tag); 
   }
}
/*======================================================================*/  
// send a vector of strings
void PbMpi::receiveStringVector(std::vector<std::string>& v, int source, int tag, PbMpi::Comm comm)
{
   // send size
   int stringVectorSize {0};
   PbMpi::Recv(comm, &stringVectorSize, 1, PbMpi_INT, source, tag);
   //comm.Recv(&stringVectorSize, 1, MPI::INT, source, tag); 

   v.clear();
   v.resize(stringVectorSize);

   if(stringVectorSize>0)
   {
      std::vector<int> singleStringSizes(stringVectorSize+1);

      PbMpi::Recv(comm, &singleStringSizes[0], stringVectorSize+1, PbMpi_INT, source, tag); 
      //comm.Recv(&singleStringSizes[0], stringVectorSize+1, MPI::INT, source, tag); 

      int nofChars = singleStringSizes[stringVectorSize];
      std::vector<char> charVector(nofChars);

      PbMpi::Recv(comm, &charVector[0], nofChars, PbMpi_CHAR, source, tag); 
      //comm.Recv(&charVector[0], nofChars, MPI::CHAR, source, tag); 

      int pos=0;
      for(int i=0; i<stringVectorSize; i++)
         for(int j=0; j<singleStringSizes[i]; j++)
            v[i].push_back(charVector[pos++]);      
   }
}

#endif //PbMpi_H
