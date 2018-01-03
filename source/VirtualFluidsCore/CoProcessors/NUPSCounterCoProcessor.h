/*
*  NUPSCounterCoProcessor.h
*
*  Created on: 03.05.2011
*  Author: K. Kucher
*/

#ifndef NUPSCOUNTERCoProcessor_H_
#define NUPSCOUNTERCoProcessor_H_

#include <memory>

#include "CoProcessor.h"
#include "basics/utilities/UbTiming.h"

class Communicator;
class Grid3D;
class UbScheduler;

class NUPSCounterCoProcessor: public CoProcessor
{
public:
   NUPSCounterCoProcessor(std::shared_ptr<Grid3D> grid, std::shared_ptr<UbScheduler> s, int numOfThreads, std::shared_ptr<Communicator> comm);
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
   std::shared_ptr<Communicator> comm;
};


#endif 
