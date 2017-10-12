/*
*  NUPSCounterCoProcessor.h
*
*  Created on: 03.05.2011
*  Author: K. Kucher
*/

#ifndef NUPSCOUNTERCoProcessor_H_
#define NUPSCOUNTERCoProcessor_H_

#include "CoProcessor.h"
#include "Communicator.h"
#include "basics/utilities/UbTiming.h"
#include <boost/timer.hpp>

class NUPSCounterCoProcessor: public CoProcessor {
public:
   NUPSCounterCoProcessor(Grid3DPtr grid, UbSchedulerPtr s, int numOfThreads, CommunicatorPtr comm);
   virtual ~NUPSCounterCoProcessor();
   void process(double step);
protected:
   void collectData(double step);
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
