#include <iostream>
#include <string>

#include "VirtualFluids.h"

using namespace std;


void run(string configname)
{
   try
   {
      vf::basics::ConfigurationFile   config;
      config.load(configname);

      string path = config.getString("path");
      int    step = config.getValue<int>("step");
      int    numberOfProcesses = config.getValue<int>("numberOfProcesses");
      
      SPtr<vf::mpi::Communicator> comm = vf::mpi::MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      SPtr<Grid3D> grid(new Grid3D(comm));

      CheckpointConverter converter(grid, path, comm);
      converter.convert(step, numberOfProcesses);
   }
   catch (std::exception& e)
   {
      cerr << e.what() << endl << flush;
   }
   catch (std::string& s)
   {
      cerr << s << endl;
   }
   catch (...)
   {
      cerr << "unknown exception" << endl;
   }

}
int main(int argc, char* argv[])
{
   if (argv != NULL)
   {
      if (argv[1] != NULL)
      {
         run(string(argv[1]));
      }
      else
      {
         cout << "Configuration file is missing!" << endl;
      }
   }

}

