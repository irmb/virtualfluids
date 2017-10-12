/*
*  TimeseriesWriterCoProcessor.h
*
*  Created on: 08.05.2013
*  Author: uphoff
*/

#include "TimeseriesCoProcessor.h"


#include <iostream>
#include <fstream>

using namespace std;

TimeseriesCoProcessor::TimeseriesCoProcessor(Grid3DPtr grid, UbSchedulerPtr s,
                                                             IntegrateValuesHelperPtr h1,
                                                             const std::string& path, CommunicatorPtr comm)
                                                             : CoProcessor(grid, s),                                                
                                                               h1(h1),
                                                               path(path),
                                                               comm(comm)
{
   if (comm->getProcessID() == comm->getRoot())
   {
      std::ofstream ostr;
      //fname = path+"/timeseries/timeseries"+UbSystem::toString(grid->getTimeStep())+".csv";
      fname = path+".csv";
      UBLOG(logINFO, "TimeseriesWriterCoProcessor::fname:" << fname);
      ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
      if(!ostr)
      { 
         ostr.clear();
         string path = UbSystem::getPathFromString(fname);
         if (path.size()>0) { UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app); }
         if(!ostr) throw UbException(UB_EXARGS,"couldn't open file "+fname);
      }
      ostr << "step;rho;vx;vy;vz;volume\n";
      ostr.close();
      UBLOG(logINFO, "TimeseriesWriterCoProcessor::Constructor:end");
   }
}
//////////////////////////////////////////////////////////////////////////
TimeseriesCoProcessor::~TimeseriesCoProcessor() 
{
}
//////////////////////////////////////////////////////////////////////////
void TimeseriesCoProcessor::process(double step)
{
   if(scheduler->isDue(step) )
      collectData(step);
}
//////////////////////////////////////////////////////////////////////////
void TimeseriesCoProcessor::collectData(double step)
{
   h1->calculateMQ();

   UBLOG(logDEBUG3, "TimeseriesWriterCoProcessor::update:" << step);

   if (comm->getProcessID() == comm->getRoot())
   {
      int istep = static_cast<int>(step);
      std::ofstream ostr;
      double cellsVolume = h1->getCellsVolume();

      double rho=(h1->getRho())/cellsVolume;
      double vx= (h1->getVx1())/cellsVolume;
      double vy= (h1->getVx2())/cellsVolume;
      double vz= (h1->getVx3())/cellsVolume;
      double volume = cellsVolume;

      ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
      if(!ostr)
      { 
         ostr.clear();
         string path = UbSystem::getPathFromString(fname);
         if(path.size()>0){ UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);}
         if(!ostr) throw UbException(UB_EXARGS,"couldn't open file "+fname);
      }

      ostr << istep << ";" << rho <<";" << vx << ";" << vy << ";" << vz << ";" << volume << "\n";
      ostr.close();
   }
}
