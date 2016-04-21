#include "NUPSCounterCoProcessor.h"


NUPSCounterCoProcessor::NUPSCounterCoProcessor(Grid3DPtr grid, UbSchedulerPtr s, int numOfThreads, CommunicatorPtr comm)
                                                   : CoProcessor(grid, s),
                                                     numOfThreads(numOfThreads),
                                                     comm(comm),
                                                     nup(0),
                                                     nup_t(0),
                                                     nupsStep(0.0)
{
   if (comm->getProcessID() == comm->getRoot())
   {
      timer.resetAndStart();

      double nop = comm->getNumberOfProcesses();
      int minInitLevel = grid->getCoarsestInitializedLevel();
      int maxInitLevel = grid->getFinestInitializedLevel();
      int gl = 2;
      UbTupleInt3 blocknx = grid->getBlockNX();
      //double nod = (val<1>(blocknx)+gl) * (val<2>(blocknx)+gl) * (val<3>(blocknx)+gl);
      double nod = (double)(val<1>(blocknx)) * (double)(val<2>(blocknx)) * (double)(val<3>(blocknx));
      nup = 0;

      for(int level = minInitLevel; level<=maxInitLevel; level++)
      {
         int nob = grid->getNumberOfBlocks(level);
         nup_t += (double)(1<<level) * nob * nod;
      }
      nup = nup_t / nop;
   }
}
//////////////////////////////////////////////////////////////////////////
NUPSCounterCoProcessor::~NUPSCounterCoProcessor() 
{
}
//////////////////////////////////////////////////////////////////////////
void NUPSCounterCoProcessor::process(double step)
{
   if(scheduler->isDue(step) )
      collectData(step);
}
//////////////////////////////////////////////////////////////////////////
void NUPSCounterCoProcessor::collectData(double step)
{
   if (comm->getProcessID() == comm->getRoot())
   {
      double time = timer.stop();
      //double time = timer.elapsed();
      //std::ofstream ostr;
      //std::string fname = path;
      //ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
      //if(!ostr)
      //{ 
      //   ostr.clear(); 
      //   std::string path = UbSystem::getPathFromString(fname);
      //   if(path.size()>0){ UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);}
      //   if(!ostr) throw UbException(UB_EXARGS,"couldn't open file "+fname);
      //}
      double nups_t = nup_t*(step-nupsStep)/time;
      double nups = nup*(step-nupsStep)/time;//timer.getTotalTime();
      double tnups = nups/(double)numOfThreads;
      //ostr << nups << std::endl;
      //ostr.close();
      UBLOG(logINFO, "Calculation step = "<<step);
      UBLOG(logINFO, "Total performance = "<<nups_t<<" NUPS");
      UBLOG(logINFO, "Performance per process = "<<nups<<" NUPS");
      UBLOG(logINFO, "Performance per thread = "<<tnups<<" NUPS");
      UBLOG(logINFO, "Time for " << step-nupsStep <<" steps = "<< time <<" s");
      //timer.restart();
      nupsStep = step;
      timer.resetAndStart();
   }
}
