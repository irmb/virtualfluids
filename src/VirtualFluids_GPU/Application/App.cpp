#include "LBM/LB.h"
#include "App.h"
#include "Communication/Communicator.h"
#include "LBM/Simulation.h"

App* App::instanz = 0;
App* App::getInstanz()
{
   if( instanz == 0 )
      instanz = new App();
   return instanz;
}
void App::run(std::string &cstr)
{
	try
	{
      Simulation sim;
      sim.init(cstr);
      sim.run();
	}
   catch (std::string e)
	{
      std::cout << e << std::flush;
      //MPI_Abort(MPI_COMM_WORLD, -1);
	}
}
