#include <iostream>
#include <string>

#include <vfluids.h>

#include "fbond.h"
#include "Version.h"


using namespace std;


int agent_main();

int main(int argc, char* argv[])
{
   int returnval = 0;
   try
   {
      bond::init();
      returnval = agent_main();
      bond::finalize();
   }
   catch(std::runtime_error& e)
   {
      std::cerr<<"\nRUNTIME ERROR: "<<e.what()<<"\n"<<std::endl;
   }
   catch(...)
   {
      std::cerr<<"unknown error"<<std::endl;
   }
   return returnval;
}


int agent_main()
{
   cout<<"\n=== bond lib info:\n"<<bond::Version::info()<<"\n===\n\n";

   string pathname = "/work/koskuche/scratch/bond_benchmark";

   // try to work around a bug in mpich (at least mpich2-1.4.1p1 and mpich2-1.5a1)
   int _mpiInitialized = (int)false;
   MPI_Initialized(&_mpiInitialized);
   if(!_mpiInitialized)
   {
      MPI_Init(0, 0);    
      _mpiInitialized = true;
   }

   int mpi_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
   int mpi_size;
   MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
   cout<<"I am process "<<bond::processID()<<" of "<<bond::processCount()
      <<", bundle ID "<<bond::bundleID()<<" of "<<bond::bundleCount()
      <<", MPI rank "<<mpi_rank<<" of "<<mpi_size<<"\n";

   stringstream logFilename;
   logFilename <<  pathname + "/logfile"+UbSystem::toString(UbSystem::getTimeStamp())+"_"+UbSystem::toString(bond::bundleID())+".txt";
   UbLog::output_policy::setStream(logFilename.str());

   vector<double> data(42);
   data[0] = 123.1;
   data[data.size()-1] = -999.1;  

   vector<double> data2(42);
   data2[0] = 123.2;
   data2[data2.size()-1] = -999.2;     

   std::tr1::shared_ptr<bond::FutureReceive> receiveRequest[8];
   int id = bond::processID();

   UBLOG(logINFO, "I am process "<<bond::processID()<<" of "<<bond::processCount()
      <<", bundle ID "<<bond::bundleID()<<" of "<<bond::bundleCount()
      <<", MPI rank "<<mpi_rank<<" of "<<mpi_size);

   for(int i = 0; i < 8; i++)
   {
      UBLOG(logINFO, "bond::receiveFuture:start "<< i);
      if(i == id) continue;
      receiveRequest[i] = bond::receiveFuture(&data[0], data.size(), MPI_DOUBLE, i, 0);  
      UBLOG(logINFO, "bond::receiveFuture:end "<< i);
   }

   for(int i = 0; i < 8; i++)
   {
      UBLOG(logINFO, "bond::sendComplete:start "<< i);
      if(i == id) continue;
      bond::sendComplete(&data2[0], data2.size(), MPI_DOUBLE, i, 0); 
      UBLOG(logINFO, "bond::sendComplete:end "<<i);
   }

   for(int i = 0; i < 8; i++)
   {
      UBLOG(logINFO, "receiveRequest->complete:start "<<i);
      if(i == id) continue;
      receiveRequest[i]->complete();  
      UBLOG(logINFO, "receiveRequest->complete:end "<<i);
   }

   return 0;
}
