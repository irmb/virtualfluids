/*
*  NUPSCounterCoProcessor.h
*
*  Created on: 03.05.2011
*  Author: K. Kucher
*/

#ifndef NUPSCOUNTERCoProcessor_H_
#define NUPSCOUNTERCoProcessor_H_

#include <PointerDefinitions.h>

#include "CoProcessor.h"
#include "basics/utilities/UbTiming.h"

class Communicator;
class Grid3D;
class UbScheduler;

class NUPSCounterCoProcessor: public CoProcessor
{
public:
   NUPSCounterCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, int numOfThreads, SPtr<Communicator> comm);
   virtual ~NUPSCounterCoProcessor();

   void process(double step)override;

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
   SPtr<Communicator> comm;
};


#endif 
