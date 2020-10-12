//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __         
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |        
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |        
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |        
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____    
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|   
//      \    \  |    |   ________________________________________________________________    
//       \    \ |    |  |  ______________________________________________________________|   
//        \    \|    |  |  |         __          __     __     __     ______      _______    
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)   
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______    
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/   
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file NUPSCounterCoProcessor.cpp
//! \ingroup CoProcessors
//! \author Konstantin Kutscher
//=======================================================================================

#include "NUPSCounterCoProcessor.h"

#include "Communicator.h"
#include "UbScheduler.h"
#include "Grid3D.h"

NUPSCounterCoProcessor::NUPSCounterCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, int numOfThreads, SPtr<Communicator> comm)
                                                   : CoProcessor(grid, s),
                                                     numOfThreads(numOfThreads),
                                                     nup(0),
                                                     nup_t(0),
                                                     nupsStep(0.0),
                                                     comm(comm)
{
   if (comm->getProcessID() == comm->getRoot())
   {
      timer.resetAndStart();

      double nop = comm->getNumberOfProcesses();
      int minInitLevel = grid->getCoarsestInitializedLevel();
      int maxInitLevel = grid->getFinestInitializedLevel();
      UbTupleInt3 blocknx = grid->getBlockNX();
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
= default;
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
      double nups_t = nup_t*(step-nupsStep)/time;
      double nups = nup*(step-nupsStep)/time;
      double tnups = nups/(double)numOfThreads;
      UBLOG(logINFO, "Calculation step = "<<step);
      UBLOG(logINFO, "Total performance = "<<nups_t<<" NUPS");
      UBLOG(logINFO, "Performance per process = "<<nups<<" NUPS");
      UBLOG(logINFO, "Performance per thread = "<<tnups<<" NUPS");
      UBLOG(logINFO, "Time for " << step-nupsStep <<" steps = "<< time <<" s");
      nupsStep = step;
      timer.resetAndStart();
   }
}
