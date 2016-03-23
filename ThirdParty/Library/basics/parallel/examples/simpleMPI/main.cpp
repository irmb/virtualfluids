#include <iostream>
#include <vector>
#include <algorithm>
#include <mpi.h>

#include <basics/utilities/UbSystem.h>
#include <basics/utilities/UbException.h>
#include <basics/utilities/UbLogger.h>

using namespace std;

int randomNumber () { return (rand()%100); }

struct RankSetter{
   RankSetter(int rank) : rank(rank) {}
   
   int operator()() 
   {
      return rank;
   }
  
   int rank;
} /*rankSetter*/;


//////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
   MPI::Init(argc, argv);
   MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS); 

   try
   {  
      MPI::Intracomm comm = MPI::COMM_WORLD;
      
      int rank = comm.Get_rank();
      
      vector<int> sendData(1000,0);
      generate(sendData.begin(), sendData.end(), RankSetter(rank+1) );

      vector<int> recvData(1000,0);

      if(rank==0)
      {
         UBLOG(logINFO,"rank="<<rank<<" - recv request");
         MPI::Request request = comm.Irecv(&recvData[0], (int)recvData.size(), MPI::INT, 1, 100);
         UBLOG(logINFO,"rank="<<rank<<" - sendData");
         comm.Ssend(&sendData[0],(int)sendData.size(), MPI::INT, 1, 100);
         sendData.back() = 999;

         UBLOG(logINFO,"rank="<<rank<<" - Wait");
         request.Wait();
         UBLOG(logINFO,"rank="<<rank<<" - all data received, last = "<<recvData.back());
      }
      else if(rank == 1)
      {
         UbSystem::sleepS(5);
         UBLOG(logINFO,"rank="<<rank<<" - recv request");
         MPI::Request request = comm.Irecv(&recvData[0],(int)recvData.size(), MPI::INT, 0, 100);
         
         request.Wait();
         UBLOG(logINFO,"rank="<<rank<<" - all data received, last = "<<recvData.back());

         UbSystem::sleepS(5);
         UBLOG(logINFO,"rank="<<rank<<" - sendData");
         comm.Ssend(&sendData[0],(int)sendData.size(), MPI::INT, 0, 100);
         sendData.back() = 999;
         UBLOG(logINFO,"rank="<<rank<<" - data sent");
      }
      else 
      {
         throw UB_THROW( UbException(UB_EXARGS,"only two ranks allwoed") );
      }

      UBLOG(logINFO,"rank="<<rank<<" barrier start");
      MPI::COMM_WORLD.Barrier();
      UBLOG(logINFO,"rank="<<rank<<" barrier done ");

   }
   catch(const std::exception& e)
   {
      UBLOG2(  logERROR,std::cerr, "caught exception:" );
      UBLOG2(  logERROR,std::cerr, "type: " << typeid(e).name() );
      UBLOG2ML(logERROR,std::cerr, "what: " << e.what() );
   }
   catch(MPI::Exception e)
   { 
      UBLOG2ML(logERROR,std::cerr, "caught exception:" << e.Get_error_string());

      MPI::COMM_WORLD.Abort(99); 
   } 
   catch(...)
   {
      UBLOG2(logERROR,std::cerr,"Verdammte Scheisse - mal wieder Mist gebaut!");
   }

   MPI::Finalize();

   return 0;
}

