/*
*  NUPSCounterPostprocessor.h
*
*  Created on: 03.05.2011
*  Author: K. Kucher
*/

#ifndef NUPSCOUNTERPOSTPROCESSOR_H_
#define NUPSCOUNTERPOSTPROCESSOR_H_

#include "Postprocessor.h"
#include "Communicator.h"
#include "basics/utilities/UbTiming.h"
#include <boost/timer.hpp>

class NUPSCounterPostprocessor: public Postprocessor {
public:
   NUPSCounterPostprocessor(Grid3DPtr grid, UbSchedulerPtr s, int numOfThreads, CommunicatorPtr comm);
   virtual ~NUPSCounterPostprocessor();
   void update(double step);
protected:
   void collectPostprocessData(double step);
   UbTimer timer;
   //boost::timer timer;
   int numOfThreads;
   double numberOfNodes;
   double numberOfBlocks;
   double nup;
   double nup_t;
   double nupsStep;
   CommunicatorPtr comm;
};


#endif 
